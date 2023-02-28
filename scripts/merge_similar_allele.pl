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
use IPC::Open2;
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

use realign_alts qw/aln_halign/;
use read_config qw/read_config_yaml/;
use Identity;

our $config = read_config_yaml("$Bin/../config.yaml");

my $seqs2identity = $$config{'seqs2identity'};
my $stretcher = $$config{'stretcher'};

sub help {
    die <<EOF;
usage: perl $0 [parameters]
--invcf <file>
--outvcf <file>
--tmpdir <dir>
--sv2pav_merge_identity_threshold <int>
--threads <int>
EOF
}
my ($opt_help, $invcf, $outvcf);
my $mcl_group_threshold_identity = 0.8;
my $tmp_dir_def = '/run/user/' . `id -u`; chomp $tmp_dir_def;
our $tmp_dir = $tmp_dir_def;
my $threads = 32;
our $verb = 0;
our $debug = 0;
my $max_refalts_threshold = 16;

GetOptions (
        'help|h!' => \$opt_help,
        't|threads=i' => \$threads,
        'sv2pav_merge_identity_threshold=s' => \$mcl_group_threshold_identity,
        'i|invcf=s' => \$invcf,
        'o|outvcf=s' => \$outvcf,
        'T|tmpdir=s' => \$tmp_dir,
        'max_refalts_threshold=i' => \$max_refalts_threshold,
        'debug!' => \$debug,
);


&help() unless $invcf and $outvcf;

$tmp_dir = $tmp_dir_def if $tmp_dir eq 'Default';
realign_alts::init();

my $I = open_in_fh($invcf);
my $O = open_out_fh($outvcf);



while ( <$I> ){
    if ( /^#/ ){
        print $O $_;
        last if /^#CHROM/;
        next;
    }
}

my $lines_cpx = MCE::Shared->array();

if ($threads==1) {
    while(my $line = <$I>) {
        my ($result) = &prase_line($line);
        say $O $result if $result;
    }
} else {
    mce_flow_f {
            chunk_size => 1,
            max_workers => $threads,  
            gather => MCE::Candy::out_iter_fh($O)
        },
        \&flt_mce,
        $I;
}

close $I;

my $lines_cpx_len = $lines_cpx->len();
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
    my ($result, $is_cpx) = &prase_line($$chunk_ref[0], 1, 0);
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
    #say STDERR "line::: ", $line;
    my @F = split(/\t/, $line);
    say STDERR "$F[0] $F[1]" if $debug; 
    my $ref = $F[3];
    my @alts = split(/,/, $F[4]);
    my @ref_alts = ($ref, @alts);
    my $max_refaltsi = scalar(@ref_alts)-1;
    if ($force_cpx==0 and $max_refaltsi >= $max_refalts_threshold) {
        say STDERR "Too much alleles at $F[0]:$F[1] ($max_refaltsi)! Keep this loc unchanged";
        return($line, 1);
    }
    #my $random = time() . "_" . rand();
    #my $temp_prefile = "$tmp_dir/$random";
    my $groups = &get_identity_halign(\@ref_alts, $threads);
    my ($replace, $newseqs) = &cal_maxlen_allele(\@ref_alts, $groups);
    my $newline = &buildline(\@F, $replace, $newseqs);
    return($newline, 0);
}

sub buildline {
    my ($F, $replace, $newseqs) = @_;
    my @newline = ( $F->@[0,1,2,3,4,5,6,7,8] );
    my $max_allei = scalar(@$newseqs)-1;
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
    my $obj = new Identity({nseqs=>$len});
    #return {"@ids"=>-1} if length( $$seqs{$ids[0]} ) > 1000; ######## DEBUG
    my @pairs;
    foreach my $i1 (0..$len) {
        foreach my $i2 ($i1..$len) {
            next if $i2<=$i1;
            push @pairs, "$i1:$i2";
        }
    }
    if($threads==1) {
        foreach my $pair (shuffle @pairs) {
            &run_get_identity_halign_pair($pair, $seqs, $obj);
        }
    } else {
        my $mutex = MCE::Mutex->new;
        $obj->{lock} = $mutex;
        mce_flow {
                chunk_size => 1,
                max_workers => $threads, 
            },
            sub{ &run_get_identity_halign_pair($_, $seqs, $obj, $mutex) },
            \@pairs;
    }
    my $groups = $obj->get_final_groups();
    return($groups);
}

sub run_get_identity_halign_pair {
    my ($pair, $seqs, $obj) = @_;
    my ($id1, $id2) = split(/:/, $pair);
    my $seq1 = $$seqs[$id1];
    my $len1 = length($seq1);
    return if $obj->needs_run($id1, $id2) == 0;
    say STDERR "get_identity_halign : $id1 $id2" if $debug;
    my $seq2 = $$seqs[$id2];
    my $len2 = length($seq2);
    my $lendiff = abs($len1-$len2);
    if ($seq1 eq '' and $seq2 eq '') {
        say STDERR "WARN: seq for $id1 and $id2 Both empty, identity set to 1" if $verb;
        $obj->addpair($id1, $id2);
        return;
    } elsif($seq1 eq '' or $seq2 eq '') {
        say STDERR "WARN: seq for $id1 or $id2 is empty! identity set to 0" if $verb;
        return;
    } elsif ($lendiff>0.5*$len1 or $lendiff>0.5*$len2) {
        say STDERR "WARN: lendiff for $id1 $id2 to large: $lendiff! identity set to 0" if $verb;
        return;
    } else {
        my ($aln_exit_status, $aln_seqs) = aln_halign('HAlignC', [$seq1, $seq2], 1);
        if ($aln_exit_status != 0) {
            next;
        }
        my $identity = &seq2same($aln_seqs);
        if ($identity > $mcl_group_threshold_identity) {
            $obj->addpair($id1, $id2);
        }
        return;
    }
}

sub seq2same {
    my ($seqs) = @_;
    tie my @s1, 'Tie::CharArray', $$seqs[0];
    tie my @s2, 'Tie::CharArray', $$seqs[1];
    my $len1 = scalar(@s1);
    my $len2 = scalar(@s2);
    if ($len1 != $len2) {
        confess "????";
    }
    my $same=0;
    my $lenr1 = 0;
    my $lenr2 = 0;
    for(my $i=0; $i < $len1; $i++) {
        my $b1 = $s1[$i];
        my $b2 = $s2[$i];
        if ($b1 eq $b2 and $b1 ne '-') {
            $same++;
        }
        $lenr1++ if $b1 ne '-';
        $lenr2++ if $b2 ne '-';
    }
    my $meanlenr = ($lenr1+$lenr2)/2;
    my $identity = $same/$meanlenr;
    return($identity);
}






