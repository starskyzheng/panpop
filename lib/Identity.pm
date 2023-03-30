#!/usr/bin/perl

package Identity;

require strict;
require warnings;
use v5.10;
no warnings qw( experimental::smartmatch );
use Carp qw/confess carp/; # carp=warn;confess=die

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Data::Dumper;
use MCE::Mutex;
use Carp;

use zzIO;

our $use_mcl = 1;
if($use_mcl==1) {
    #use Algorithm::MCL;
    use zzmcl;
}

sub new {
    my ($class, $args) = @_;
    my %args = %$args;
    confess 'Identity::new usage : $args{nseqs}' unless
            exists $args{seqs_maxi};
    my $self={
        _seqs_maxi => $args{seqs_maxi},
        groups => [],
        pair2indentities => {},
        parallel => 0,
        mcl_group_threshold_diff => $args{mcl_group_threshold_diff} // 20,
        mcl_group_threshold_identity => $args{mcl_group_threshold_identity} // 0.8,
        seqs => $args{seqs} // undef,
        lock => undef,# => MCE::Mutex->new,
        use_mcl => $args{use_mcl} // 0,
    };
    bless($self, $class || ref($class));
    $self->{mcl_group_threshold_diff_percent} = 1 - $self->{mcl_group_threshold_identity};
    $self->cal_seq_lens() if defined $self->{seqs};
    return $self;
}

sub cal_seq_lens {
    my ($self)  = @_;
    my $seqs = $self->{seqs};
    my $seqs_maxi = $self->{_seqs_maxi};
    my $seq_lens = [];
    for(my $i=0; $i <= $seqs_maxi; $i++) {
        my $seq = $$seqs[$i];
        my $len = length($seq);
        $$seq_lens[$i] = $len;
    }
    $self->{seq_lens} = $seq_lens;
}

sub addpair {
    my ($self, $id1, $id2, $indentity) = @_;
    if($use_mcl==1) {
        $self->addpair_mcl($id1, $id2, $indentity);
    } else {
        $self->addpair_simple($id1, $id2, $indentity) if $indentity > 0;
    }
}

sub init_pair2indentities {
    # set all pair to 0
    my ($self) = @_;
    my $nseqs = $self->{_seqs_maxi};
    my $pair2indentities = $self->{pair2indentities};
    for(my $i=0; $i < $nseqs; $i++) {
        for(my $j=$i+1; $j < $nseqs; $j++) {
            next if $i>=$j;
            my $pair = "$i:$j";
            $$pair2indentities{$pair} = 0;
        }
    }
}

sub addpair_mcl {
    my ($self, $id1, $id2, $indentity) = @_;
    $indentity //= 1;
    # say STDERR "add pair: $id1 $id2 $indentity";
    # must: $id1 < $id2
    ($id1, $id2) = ($id2, $id1) if ($id1 > $id2);
    my $pair = "$id1:$id2";
    my $pair2indentities = $self->{pair2indentities};
    $$pair2indentities{$pair} = $indentity;
    # say STDERR Dumper $pair2indentities;
}

sub addpair_simple {
    my ($self, $id1, $id2, $indentity) = @_;
    my $groups = $self->{groups};
    ($id1, $id2) = ($id2, $id1) if ($id1 > $id2);
    my $pair = "$id1:$id2";
    my $pair2indentities = $self->{pair2indentities};
    $$pair2indentities{$pair} = $indentity;
    return if $indentity <= 0;
    ADDPAIRREDO:
    my $id1g = $self->is_in_group($id1);
    my $id2g = $self->is_in_group($id2);
    foreach my $idg ($id1g, $id2g) {
        if ( ref($idg) ne '' ) { # id in multiple groups
            $self->merge_groups($idg);
            redo ADDPAIRREDO;
        }
    }
    if ($id1g == -1 and $id2g == -1) {
        push @$groups, [$id1, $id2];
    } elsif ( $id1g != -1 and $id2g == -1 ) {
        push $$groups[$id1g]->@*, $id2;
    } elsif ( $id1g == -1 and $id2g != -1 ) {
        push $$groups[$id2g]->@*, $id1;
    } elsif ( $id1g != -1 and $id2g != -1 and $id1g != $id2g ) {
        # merge two groups
       $self->merge_groups( [$id1g, $id2g] );
    }
}

sub merge_groups {
    my ($self, $tomerge) = @_;
    my $groups = $self->{groups};
    my @tomerge = sort{$b<=>$a} @$tomerge;
    my $dst = pop @tomerge; # min
    foreach my $ig (@tomerge) {
        push $$groups[$dst]->@*, $$groups[$ig]->@*;
        splice @$groups, $ig, 1;
    }
}

sub needs_run {
    my ($self, $id1, $id2, $check_len) = @_;
    if($check_len>0) {
        my $seq_lens = $self->{seq_lens};
        my $len1 = $$seq_lens[$id1] // die;
        my $len2 = $$seq_lens[$id2] // die;
        my $diff = abs($len1-$len2);
        my $min_len = $len1 < $len2 ? $len1 : $len2;
        if ($diff > $self->{mcl_group_threshold_diff}) {
            $self->addpair($id1, $id2, 0);
            return 0;
        } elsif($diff > $self->{mcl_group_threshold_diff_percent} * $min_len) {
            $self->addpair($id1, $id2, 0);
            return 0;
        } elsif($len1 == 0 and $len2  == 0) {
            $self->addpair($id1, $id2, 1);
            return 0;
        } elsif($len1 == 0 or $len2  == 0) {
            $self->addpair($id1, $id2, 0);
           return 0;
        }
    }
    my $id = "$id1:$id2";
    my $pair2indentities = $self->{pair2indentities};
    if($use_mcl==1) { # mcl mode
        return 0 if exists $$pair2indentities{$id};
        return 1; # force need run
    } else { # simple mode
        return 0 if exists $$pair2indentities{$id};
        my $id1g = $self->is_in_group($id1);
        my $id2g = $self->is_in_group($id2);
        if ( ref($id1g) eq '' and ref($id2g) eq '' and 
                $id1g eq $id2g and 
                $id1g != -1 and $id2g != -1 ) {
            return 0; # no need run, in same group
        } else {
            return 1; # need run
        }
    }
}



sub is_in_group {
    my ($self, $id) = @_;
    my $groups = $self->{groups};
    my $max_groupsi = scalar(@$groups)-1;
    my @ingroups;
    foreach my $i (0..$max_groupsi) {
        if( $id ~~ $$groups[$i]->@* ) {
            push @ingroups, $i;
        }
    }
    if (! @ingroups) {
        return -1;
    }elsif (scalar(@ingroups) >= 2) {
        return(\@ingroups);
    } else {
        return($ingroups[0])
    }
}

sub get_final_groups {
    my ($self) = @_;
    if($use_mcl==1) {
        return $self->get_final_groups_mcl();
    } else {
        return $self->get_final_groups_simple();
    }
}

sub get_final_groups_simple {
    my ($self) = @_;
    my $groups = $self->{groups};
    my $nseqs = $self->{_seqs_maxi};
    my $max_groupsi = scalar(@$groups)-1;
    my %ids_ingroup;
    foreach my $i (0..$max_groupsi) {
        foreach my $id ($$groups[$i]->@*) {
            $ids_ingroup{$id}++;
        }
    }
    # single id with no group
    foreach my $id (0..$nseqs) {
        if (! exists $ids_ingroup{$id}) {
            push @$groups, [$id];
        }
    }
    return $groups;
}


sub get_final_groups_mcl {
    my ($self) = @_;
    my $pair2indentities = $self->{pair2indentities};
    # say STDERR Dumper $pair2indentities;
    if($self->{parallel}==1) {
        my $pair2indentities2 = $pair2indentities->export({ unbless => 1 });
        $self->{pair2indentities} = $pair2indentities2;
        $pair2indentities = $pair2indentities2;
        undef $pair2indentities2;
    }

    # del 0 value from pair2indentities
    foreach my $id (keys %$pair2indentities) {
        delete $$pair2indentities{$id} if $$pair2indentities{$id} == 0;
    }
    my $min_group_count = 2;
    my @identites = sort {$a<=>$b} grep {defined} values %$pair2indentities;
    my $maxi = scalar(@identites)-1;
    my $threshold_i = 0;
    my $groups;
    while(1) {
        my $threshold = $identites[$threshold_i]; # start from min
        #say STDERR "$threshold_i / $maxi : $threshold";
        if($threshold_i == $maxi) {
            if(defined $groups) {
                return $group;
            } else {
                return $self->get_final_groups_mcl_run($threshold);
            }
        }
        my $groups = $self->get_final_groups_mcl_run($threshold);
        my $group_count = scalar(@$groups);
        if($group_count >= $min_group_count) {
            return $groups;
        }
        my $add = int($maxi/10);
        $add = 1 if $add < 1;
        $threshold_i += $add;
        if($threshold_i > $maxi) {
            $threshold_i = $maxi;
        }
    }

}

sub get_final_groups_mcl_run {
    my ($self, $mcl_group_threshold_identity) = @_;
    $mcl_group_threshold_identity = $self->{mcl_group_threshold_identity} unless defined $mcl_group_threshold_identity;
    my $nseqs = $self->{_seqs_maxi};
    my $max_refaltsi = $nseqs;
    my $items = $self->{pair2indentities};
    #die Dumper $items;
    #my ($items, $max_refaltsi) = @_;
    #say STDERR '$items: ' . Dumper $items;
    #my $mcl1 = Algorithm::MCL->new();
    my $mcl1 = zzmcl->new({I=>2});
    my $need_cal_mcl=0;
    my (%ids, %ids_mcl);
    foreach my $pair (keys %$items) {
        $pair =~ /^([^>]+):([^:]+)$/ or next;
        my $identity = $$items{$pair};
        $identity = 0 if $identity < 0;
        $identity = 1 if $identity > 1;
        my ( $id1, $id2) = ( int($1), int($2) );
        next if $identity < $mcl_group_threshold_identity;
        #say STDERR join "\t", "add edge: ", $id1, $id2, $identity;
        $mcl1->addEdge($id1, $id2, $identity);
        $need_cal_mcl++;
        #$ids_mcl{$id1}++; $ids_mcl{$id2}++;
    }
    my $mcl_clusters=[];
    $mcl_clusters = $mcl1->run() if $need_cal_mcl>0;
    #say STDERR "mcl_clusters: " . Dumper $mcl_clusters;
    my @groups;
    foreach my $mcl_cluster ( @$mcl_clusters ) {
        #print "Cluster size: ". scalar @$cluster. "\n";
        my @ids_now;
        push @ids_now, $_ foreach @$mcl_cluster;
        push @groups, \@ids_now;
        $ids_mcl{$_}++ foreach @ids_now;
    }
    # single id with no group
    foreach my $id (0..$max_refaltsi) {
        unless (exists $ids_mcl{$id}) {
            push @groups, [int($id)];
            #say STDERR $id;
        }
    }
    #say STDERR '@groups: ' . Dumper \@groups;
    $self->{groups} = \@groups;
    #die;
    return \@groups;
}



1;
