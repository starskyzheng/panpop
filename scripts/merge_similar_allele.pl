#!/usr/bin/env perl
#===============================================================================
#
#         FILE: sv2snp_indel.pl
#
#        USAGE: ./sv2snp_indel.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zeyu Zheng (LZU), zhengzy2014@lzu.edu.cn
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 2/17/2022 11:54:32 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
no warnings qw( experimental::smartmatch );
use v5.24;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Carp qw/confess carp/; # carp=warn;confess=die
use zzIO;
use Data::Dumper;
use Tie::CharArray;
use MCE::Flow;
use MCE::Candy;
use MCE::Shared;
use MCE::Mutex;
use Getopt::Long;
use MCE::Channel;
#use Coro::Generator;
use File::Temp;
use List::Util qw/ max min shuffle /;
use File::Basename;
use Tie::CharArray;

use realign_alts qw/aln_halign process_alts_aln_only/;
use read_config qw/read_config_yaml/;
use Identity;

our $config = read_config_yaml("$Bin/../config.yaml");

my $stretcher = $$config{'stretcher'};

sub help {
    die <<EOF;
usage: perl $0 [parameters]
--invcf <file>
--outvcf <file>
--tmpdir <dir>
--sv2pav_merge_identity_threshold <int>
--sv2pav_merge_diff_threshold <int>
--threads <int>
--type <int>
type bitcode:
    0 for simple mode
    1 for mcl mode
    2 for fast and ambigous
EOF
}
my ($opt_help, $invcf, $outvcf);
my $mcl_group_threshold_identity = 0.8;
my $mcl_group_threshold_diff = 20;
my $tmp_dir_def = '/run/user/' . `id -u`; chomp $tmp_dir_def;
our $tmp_dir = $tmp_dir_def;
my $threads = 32;
our $verb = 0;
our $debug = 0;
our $type_bitcode = 0;
my $max_refalts_threshold = 16; # too big may slow down speed

my $ARGVs = join " ", $0, @ARGV; 
GetOptions (
        'help|h!' => \$opt_help,
        't|threads=i' => \$threads,
        'sv2pav_merge_identity_threshold=s' => \$mcl_group_threshold_identity,
        'sv2pav_merge_diff_threshold=i' => \$mcl_group_threshold_diff,
        'i|invcf=s' => \$invcf,
        'o|outvcf=s' => \$outvcf,
        'T|tmpdir=s' => \$tmp_dir,
        'max_refalts_threshold=i' => \$max_refalts_threshold,
        'debug!' => \$debug,
        'v|verb!' => \$verb,
        'type=i' => \$type_bitcode,
);

my $use_mcl = $type_bitcode & 1;
&help() unless $invcf and $outvcf;

$tmp_dir = $tmp_dir_def if $tmp_dir eq 'Default';
realign_alts::init();

my $I = open_in_fh($invcf);
my $O = open_out_fh($outvcf);



while ( <$I> ){
    if ( /^#/ ) {
        if (/^#CHROM/) {
            say $O qq(##CommandLine="$ARGVs");
            print $O $_;
            last;
        } else {
            print $O $_;
        }
        next;
    }
}


my $lines_cpx;

my $not_use_cpx = $type_bitcode & 2 ? 1 : 0; # 2 for fast and not use cpx; not_use_cpx==1 means force_cpx==1
if ($threads==1) {
    $lines_cpx = zzarray->new();
    while(my $line = <$I>) {
        my ($result, $is_cpx) = &prase_line($line, 1, $not_use_cpx);
        if ($is_cpx==1) {
            $lines_cpx->push($result);
        } else {
            say $O $result if $result;
        }
    }
} else {
    $lines_cpx = MCE::Shared->array();
    MCE::Flow->init(
        chunk_size => 1,
        max_workers => $threads,
        gather => MCE::Candy::out_iter_fh($O)
    );
    mce_flow_f \&flt_mce, $I;
    MCE::Flow->finish();
    my $lines_cpx2 = $lines_cpx->export({ unbless => 1 });
    undef $lines_cpx;
    $lines_cpx = $lines_cpx2;
}

#close $O;
#close $I; # cannot close here, or MCE::Flow will not work


my $lines_cpx_len = scalar(@$lines_cpx);
if ($lines_cpx_len>0) {
    say STDERR "Now processing: lines_cpx_len lines: $lines_cpx_len !\n";
    foreach my $line (@$lines_cpx) {
        my ($result, $is_cpx) = &prase_line($line, $threads, 1);
        say $O $result if $result;
    }
}

#close $O;

say STDERR "Done.";

exit;

sub flt_mce {
    my ( $mce, $chunk_ref, $chunk_id ) = @_;
    my ($result, $is_cpx) = &prase_line($$chunk_ref[0], 1, $not_use_cpx);
    if ($is_cpx==1) {
        $lines_cpx->push($result);
        $mce->gather($chunk_id);
    } else {
        if (defined $result) {
            $mce->gather($chunk_id, $result, "\n");
        } else {
            $mce->gather($chunk_id);
        }
    }
}


sub prase_line {
    my ($line, $threads, $force_cpx) = @_;
    $threads//=1;
    $force_cpx//=0;
    chomp $line;
    return (undef, 0) if $line eq '' or $line =~ /^#/;
    #say STDERR "line::: ", $line;
    my @F = split(/\t/, $line);
    say STDERR "$F[0] $F[1]" if $debug; 
    my $ref = $F[3];
    my @alts = split(/,/, $F[4]);
    my @ref_alts = ($ref, @alts);
    # * to -
    @ref_alts = map { s/^\*$//g; $_ } map { s/^-$//g; $_ } @ref_alts;
    my $max_refaltsi = scalar(@ref_alts)-1;
    if ($force_cpx==0 and $max_refaltsi >= $max_refalts_threshold) {
        say STDERR "Too much alleles at $F[0]:$F[1] ($max_refaltsi)! push to lines_cpx" if $debug==1;
        return($line, 1);
    }
    #my $random = time() . "_" . rand();
    #my $temp_prefile = "$tmp_dir/$random";
    my $groups = &get_identity_halign(\@ref_alts, $threads);
    #die Dumper $groups;
    my ($replace, $newseqs) = &cal_maxlen_allele(\@ref_alts, $groups);
    my $newline = &buildline(\@F, $replace, $newseqs);
    return($newline, 0);
}

sub buildline {
    my ($F, $replace, $newseqs) = @_;
    my @newline = ( $F->@[0,1,2,3,4,5,6,7,8] );
    my $max_allei = scalar(@$newseqs)-1;
    foreach my $allei (1..$max_allei) {
        my $seq = $$newseqs[$allei];
        $seq = '*' if $seq eq '' or $seq eq '-';
        $$newseqs[$allei] = $seq;
    }
    my $alts = join ",", $newseqs->@[1..$max_allei];
    $newline[4] = $alts;
    $newline[4] = "*" if $newline[4] eq '';
    my $max_idi = scalar(@$F) - 1;
    foreach my $idi (9..$max_idi) {
        my $infos = $$F[$idi];
        if ($infos=~s#^(\d+)([/|])(\d+)##) {
            my $a1 = $1;
            my $sep = $2;
            my $a2 = $3;
            if ( ! exists $$replace{$a1} ) {
                say STDERR "!!!!!!";
                say STDERR Dumper $replace;
                say STDERR Dumper $F;
                die "$a1";
            }
            my $a1n = $$replace{$a1} // die "!!!!!!!!$a1 @$F";
            my $a2n = $$replace{$a2} // die;
            $infos = "$a1n$sep$a2n" . $infos;
            push @newline, $infos;
        } else {
            push @newline, $infos;
        }
    }
    my $ret = join "\t", @newline;
    return($ret);
}

sub cal_maxlen_allele {
    my ($seqs, $groups) = @_;
    my %ret;
    foreach my $ids (@$groups) {
        my $maxlen = 0;
        my $maxlen_id;
        die if @$ids==0;
        foreach my $id (@$ids) {
            my $seq = $$seqs[$id];
            my $len = length($seq);
            if ($id==0) { # ref must keep
                $maxlen = $len;
                $maxlen_id = $id;
                last;
            }
            if ($len>=$maxlen) { # keep max len allele
                $maxlen = $len;
                $maxlen_id = $id;
            }
        }
        $ret{$maxlen_id} = $ids;
    }
    die unless exists $ret{0};
    my (%replace, @newseqs);
    my $newid=0;
    foreach my $maxlen_id (sort {$a<=>$b} keys %ret) {
        push @newseqs, $$seqs[$maxlen_id];
        $replace{$maxlen_id} = $newid;
        foreach my $id ($ret{$maxlen_id}->@*) {
            $replace{$id} = $newid;
        }
        $newid++;
    }
    return(\%replace, \@newseqs);
}



sub get_identity_halign {
    my ($seqs, $threads) = @_;
    $threads//=1;
    my $len = scalar(@$seqs) - 1;
    my $obj = new Identity({seqs_maxi=>$len, seqs=>$seqs, 
                    mcl_group_threshold_diff=>$mcl_group_threshold_diff, 
                    mcl_group_threshold_identity=>$mcl_group_threshold_identity,
                    use_mcl => $use_mcl, });
    #$obj->init_pair2indentities() if($threads>1);
    &run_get_identity_halign_all($seqs, $obj, $type_bitcode & 2);
    #die Dumper $obj->{pair2indentities};
    my @pairs;
    foreach my $i1 (0..$len) {
        foreach my $i2 ($i1..$len) {
            next if $i2<=$i1;
            next if $obj->needs_run($i1, $i2, 1)==0; # no needs to run
            push @pairs, "$i1:$i2";
        }
    }
    #say STDERR "paris length: " . scalar(@pairs);
    # die Dumper \@pairs;
    if(@pairs) {
        if($threads==1) {
            foreach my $pair (shuffle @pairs) {
                &run_get_identity_halign_pair($pair, $seqs, $obj);
            }
        } else {
            my $mutex = MCE::Mutex->new;
            $obj->{lock} = $mutex;
            my $pair2indentities_ori = $obj->{pair2indentities};
            my $group_ori = $obj->{groups};
            #$obj->{pair2indentities} = MCE::Shared->hash({_DEEPLY_ => 1}, %$pair2indentities_ori);
            #$obj->{groups} = MCE::Shared->array({_DEEPLY_ => 1}, @$group_ori);
            $obj->{pair2indentities} = MCE::Shared->share({_DEEPLY_ => 1}, $pair2indentities_ori);
            $obj->{groups} = MCE::Shared->share({_DEEPLY_ => 1}, $group_ori);
            $obj->{parallel} = 1;
            say STDERR "start MCE::Flow, count: " . scalar(@pairs);
            MCE::Flow->init(
               chunk_size => 1,
               max_workers => $threads, 
            );
            mce_flow sub{ &run_get_identity_halign_pair($_, $seqs, $obj) }, \@pairs;
            MCE::Flow->finish;
        }
    }
    # die Dumper $obj->{groups};
    my $groups = $obj->get_final_groups();
    #die Dumper $groups;
    return($groups);
}

sub run_get_identity_halign_all {
    my ($alts, $obj, $fastmode) = @_;
    $fastmode //= 0;
    my $max_alts = scalar(@$alts)-1;
    my ($aln_alts, $reflen) = &process_alts_aln_only($alts, -1);
    #die Dumper $aln_alts;
    foreach my $id1 (0..$max_alts) {
        foreach my $id2 ($id1..$max_alts) {
            next if $id2<=$id1;
            next if $obj->needs_run($id1, $id2, 1)==0; # no needs to run
            my ($identity, $diff, $gap) = &seq2identity([$$aln_alts[$id1], $$aln_alts[$id2]]);
            say STDERR "get_identity_halign : $id1 $id2 $identity $diff $gap" if $verb;
            if ($identity > $mcl_group_threshold_identity and $diff < $mcl_group_threshold_diff)  {
                #$seq2identity{"$id1:$id2"} = $identity;
                $obj->addpair($id1, $id2, $identity);
            } elsif($fastmode>0) {
                $obj->addpair($id1, $id2, 0);
            }
        }
    }
}

sub run_get_identity_halign_pair {
    my ($pair, $seqs, $obj, $mutex) = @_;
    #my $mutex = $obj->{lock};
    my ($id1, $id2) = split(/:/, $pair);
    my $seq1 = $$seqs[$id1];
    my $len1 = length($seq1);
    return if $obj->needs_run($id1, $id2) == 0;
    say STDERR "get_identity_halign : $id1 $id2" if $debug;
    my $seq2 = $$seqs[$id2];
    my $len2 = length($seq2);
    my $lendiff = abs($len1-$len2);
    #my ($aln_exit_status, $aln_seqs) = aln_halign('HAlignC', [$seq1, $seq2], 1);
    my ($aln_seqs) = process_alts_aln_only([$seq1, $seq2], -1, -1);
    my ($identity, $diff, $gap) = &seq2identity($aln_seqs);
    #say STDERR "$identity";
    if ($identity > $mcl_group_threshold_identity and $diff < $mcl_group_threshold_diff) {
        $mutex->lock() if $mutex;
        $obj->addpair($id1, $id2, $identity); # 多线程会出错
        $mutex->unlock() if $mutex;
    }
    return;
}

sub seq2identity {
    my ($seqs) = @_;
    tie my @s1, 'Tie::CharArray', $$seqs[0];
    tie my @s2, 'Tie::CharArray', $$seqs[1];
    my $len1 = scalar(@s1);
    my $len2 = scalar(@s2);
    my $lendiff = abs($len1-$len2);
    if ($len1 != $len2) {
        confess "???? len not equal: $len1 $len2 @$seqs";
    }
    my $same=0;
    my $score = 0;
    my $lenr1 = 0;
    my $lenr2 = 0;
    my $gap = 0;
    my $diff_now = 0;
    for(my $i=0; $i < $len1; $i++) {
        my $b1 = $s1[$i];
        my $b2 = $s2[$i];
        if ($b1 eq $b2 and $b1 ne '-') { # A/A
            $score++;
            $same++;
        } elsif ($b1 ne $b2 and ($b1 eq '-' or $b2 eq '-')) { # A/- or -/A
            $diff_now++;
            $score -= 0.5;
            $gap++;
        } elsif ($b1 eq '-' and $b2 eq '-') { # -/-
            # nothing
        } elsif ($b1 ne $b2 and ($b1 ne '-' and $b2 ne '-')) { # A/C
            $diff_now++;
            # nothing
        } else {
            confess "???? $b1 $b2";
        }
        $lenr1++ if $b1 ne '-';
        $lenr2++ if $b2 ne '-';
        return(0, $diff_now, $gap) if $diff_now > $mcl_group_threshold_diff;
    }
    #my $meanlenr = ($lenr1+$lenr2)/2;
    my $mexlenr = $lenr1 > $lenr2 ? $lenr1 : $lenr2; # max
    my $diff = $mexlenr - $same;
    my $identity = $same/$mexlenr;
    return($identity, $diff, $gap);
}





{
    package zzarray;
    sub new {
        my $class = shift;
        my $self = [];
        bless $self, $class;
        return $self;
    }
    sub push {
        my $self = shift;
        my $item = shift;
        push $self->@*, $item;
    }
    sub len {
        my $self = shift;
        return scalar $self->@*;
    }
}
