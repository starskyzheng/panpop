#!/usr/bin/perl

package Identity;

require strict;
require warnings;
use v5.10;
no warnings qw( experimental::smartmatch );

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Data::Dumper;
use Carp;

use zzIO;

sub new {
    my ($class, $args) = @_;
    my %args = %$args;
    confess 'Identity::new usage : $args{nseqs}' unless
            exists $args{nseqs};
    my $nseqs = $args{nseqs};
    my $self={
        _nseqs => $args{nseqs},
        groups => [],
        lock=>undef,
    };
    bless($self, $class || ref($class));
    return $self;
}

sub addpair() {
    my ($self, $id1, $id2) = @_;
    my $lock = $self->{lock};
    if (defined $lock) {
        $lock->lock();
    }
    my $groups = $self->{groups};
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
    if (defined $lock) {
        $lock->unlock();
    }
}

sub needs_run {
    my ($self, $id1, $id2) = @_;
    my $id1g = $self->is_in_group($id1);
    my $id2g = $self->is_in_group($id2);
    if ( ref($id1g) eq '' and ref($id2g) eq '' and 
            $id1g eq $id2g and 
            $id1g != -1 and $id2g != -1 ) {
        return 0; # no need run
    } else {
        return 1; # need run
    }
}

sub get_final_groups {
    my ($self) = @_;
    my $groups = $self->{groups};
    my $nseqs = $self->{_nseqs};
    my $max_groupsi = scalar(@$groups)-1;
    my %ids_ingroup;
    foreach my $i (0..$max_groupsi) {
        foreach my $id ($$groups[$i]->@*) {
            $ids_ingroup{$id}++;
        }
    }
    foreach my $id (0..$nseqs) {
        if (! exists $ids_ingroup{$id}) {
            push @$groups, [$id];
        }
    }
    return $groups;
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


=old_mcl
use Algorithm::MCL;
sub cal_group_by_mcl {
    my ($items, $max_refaltsi) = @_;
    say STDERR '$items: ' . Dumper $items;
    my $mcl1 = Algorithm::MCL->new();
    my $need_cal=0;
    my (%ids, %ids_mcl);
    my @ref;
    $ref[$_] = $_ foreach (0..$max_refaltsi);
    foreach my $pair (keys %$items) {
        my $identity = $$items{$pair};
        $pair =~ /^([^>]+):([^:]+)$/ or die $pair;
        my ( $id1, $id2) = ( int($1), int($2) );
        $ids{$id1}++; $ids{$id2}++;
        next if $identity < $mcl_group_threshold_identity;
        say STDERR "add edge: $id1 $id2";
        $mcl1->addEdge(\$ref[$id1], \$ref[$id2]);
        $need_cal++;
        $ids_mcl{$id1}++; $ids_mcl{$id2}++;
    }
    my $mcl_clusters=[];
    $mcl_clusters = $mcl1->run() if $need_cal>0;
    #$mcl_clusters = [[sort keys %ids_mcl]];
    my @groups;
    say STDERR '$mcl_clusters: ' . Dumper $mcl_clusters;
    foreach my $id (keys %ids) {
        unless ($id ~~ [keys %ids_mcl]) {
            push @groups, [int($id)];
            #say STDERR $id;
        }
    }
    foreach my $mcl_cluster ( @$mcl_clusters ) {
        #print "Cluster size: ". scalar @$cluster. "\n";
        my @ids;
        push @ids, $$_ foreach @$mcl_cluster;
        push @groups, \@ids;
    }
    say STDERR '@groups: ' . Dumper \@groups;
    return \@groups;
}
=cut



1;
