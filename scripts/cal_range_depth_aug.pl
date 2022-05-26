#!/usr/bin/env perl
#===============================================================================
#
#         FILE: cal_sv_depth_by_bam.pl
#
#        USAGE: ./cal_sv_depth_by_bam.pl  
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
#      CREATED: 10/21/2021 02:16:06 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
no warnings qw( experimental::smartmatch );
use v5.24;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Carp qw/confess/;
use zzIO;
use Data::Dumper;
use MCE::Loop;
use Getopt::Long;
use read_config qw/read_config_yaml/;




sub help() {
        die "usage: perl x.pl [options]
        -v|--vcfposs <str>     input vcfposs file
        -l|--list_packpg <str> input pack_pg file
        -o|--outdir <str>      output dir
        -t                     threads number
        -h|--help              help
        --chr                  process this chromsome only
        ";
}

my $config = read_config_yaml("$Bin/../config.yaml");
my $vg = $$config{vg};

my ($vcfposs_file, $packpg_list_file, $outdir, $need_chr);
my ($opt_help);
my $threads = 10;
GetOptions (
        'help|h!' => \$opt_help,
        'v|vcfposs=s' => \$vcfposs_file,
        't=i' => \$threads,
        'l|list_packpg=s' => \$packpg_list_file,
        'o|outdir=s' => \$outdir,
        'chr=s' => \$need_chr,
);
&help() if $opt_help;
&help() unless defined $vcfposs_file and defined $packpg_list_file and defined $outdir;

$need_chr = undef if $need_chr eq 'all';

my $dp_around_bp=100;

my ($gams) = &read_packpg_list_file($packpg_list_file);
my $need_ids = [ sort keys %$gams ]; ##########
$threads = scalar(@$need_ids) if $threads > scalar(@$need_ids);
say STDERR "Start to calculate depth for: " . join(' ', @$need_ids);



my $lastpos = -1;
say STDERR "Now read vcf_poss_file: $vcfposs_file";
my ($vcfposs) = &read_vcf_poss_file($vcfposs_file);
say STDERR "Done read vcf_poss_file";

if ($threads==1) {
    foreach my $sid (@$need_ids) {
        die unless exists $$gams{$sid};
        my ($pack, $pg) = $$gams{$sid}->@*;
        my $outfile = "$outdir/$sid.depth.txt.gz";
        my $O = open_out_fh($outfile);
        say STDERR "$sid start";
        &process_dp($pack, $pg, $sid, $O);
        say STDERR "$sid Done";
    }
} else { # multi-threads
    MCE::Loop->init(
        max_workers => $threads, chunk_size => 1
    );
    mce_loop {
        my $sid = $_;
        die unless exists $$gams{$sid};
        my ($pack, $pg) = $$gams{$sid}->@*;
        my $outfile = "$outdir/$sid.depth.txt.gz";
        my $O = open_out_fh($outfile);
        say STDERR "$sid start";
        &process_dp($pack, $pg, $sid, $O);
        say STDERR "$sid Done";
    } (@$need_ids);
    MCE::Loop->finish;
}


exit;

sub read_vcf_poss_file {
    my ($infile) = @_;
    my $I = open_in_fh($infile);
    my %ret;
    #my $ret = MCE::Shared->hash;
    while(<$I>) {
        chomp;
        next unless $_;
        my ($chr, $pos, @lens) = split(/\t/, $_);
        if ( scalar(@lens)>1 ) {
            $ret{$chr}{$pos} = \@lens;
        } else { # only 1 item
            $ret{$chr}{$pos} = $lens[0];
        }
        #last if $.>100; # for debug
    }
    return(\%ret);
}

sub process_dp {
    my ($pack, $pg, $sid, $O) = @_;
    my (%counts, %sum_dps);
    die unless -e $pack;
    die unless -e $pg;
    my $cmd = "$vg depth --threads 4 --pack $pack $pg";
    if (defined $need_chr) {
        $cmd .= " --ref-path $need_chr";
    }
    #$cmd = 'zcat ~/a.gz'; ############# debug
    open(my $P, "$cmd |");
    my %buf;
    my $chr_old;
    LINEDP:while(<$P>) {
        # chr pos dp
        #  0   1   2
        chomp;
        next unless $_;
        my ($chr, $pos, $dp) = split(/\t/, $_); # 1 based
        $chr_old //= $chr;
        if ($chr_old ne $chr) {
            foreach my $end (sort {$a<=>$b} keys %buf) {
                &print1($buf{$end}, $chr_old, $end, $O);
            }
            %buf = ();
            $chr_old = $chr;
            redo LINEDP;
        }
        $lastpos = -1 if ($lastpos > $pos); # change chr
        for(;$lastpos<=$pos; $lastpos++) {
            #say("$lastpos $pos $chr"); # debug
            if (exists $$vcfposs{$chr}{$lastpos}) {
                my $lens = $$vcfposs{$chr}{$lastpos};
                if (ref($lens) eq '') {
                    $lens = [$lens];
                }
                foreach my $len (@$lens) {
                    my $end = $lastpos + $len - 1;
                    $buf{$end}{$lastpos} = [0, 0]; # count, sum_dp
                }
            }
        }
        foreach my $end(keys %buf) {
            if ($pos <= $end) {
                foreach my $start (sort {$a<=>$b} keys %{$buf{$end}}) {
                    die unless $pos >= $start;
                    $buf{$end}{$start}[0] += 1 if $dp > 3;
                    $buf{$end}{$start}[1] += $dp;
                }
            } else {
                my $toprint = delete $buf{$end};
                &print1($toprint, $chr_old, $end, $O);
            }
        }
    }
}

sub print1 {
    my ($toprint, $chr, $end, $O) = @_;
    foreach my $start (sort {$a<=>$b} keys %$toprint) {
        my $len = $end - $start + 1;
        my $sumdp = $$toprint{$start}[1];
        my $count = $$toprint{$start}[0];
        next if $sumdp==0 and $count==0;
        say $O join "\t", $chr, $start, $len, $sumdp, $count;
    }
}

sub read_packpg_list_file {
    my ($file) = @_;
    my $I = open_in_fh($file);
    my %ret;
    while(<$I>) {
        # SampleID xx.pack xx.pg
        chomp;
        next unless $_;
        next if /^#/;
        my ($sid, $pack, $pg) = split(/\s+/, $_);
        die "$pack not exists!" unless -s $pack;
        die "$pg not exists!" unless -s $pg;
        $ret{$sid} = [$pack, $pg];
    }
    return \%ret;
}

