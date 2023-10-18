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
no warnings qw( experimental::smartmatch qw );
use v5.24;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Carp qw/confess/;
use zzIO;
use Data::Dumper;
use MCE::Loop;
use Getopt::Long;
use read_config qw/read_config_yaml/;


my $config = read_config_yaml("$Bin/../config.yaml");
my $vg = $$config{vg};
my $samtools = $$config{samtools};

my ($vcfposs_file, $packpg_list_file, $outdir, $need_chr);
my ($output_dpinfo, $output_dpinfo_fh);
my ($opt_help);
my $threads = 1;
my $min_dp = 3;
my $input_format = 'vg';
my $ref_only=0;

GetOptions (
        'help|h!' => \$opt_help,
        'v|vcfposs=s' => \$vcfposs_file,
        't|threads=i' => \$threads,
        'l|list_packpg=s' => \$packpg_list_file,
        'o|outdir=s' => \$outdir,
        'chr=s' => \$need_chr,
        'input_format=s' => \$input_format,
        'output_avgdpinfo=s' => \$output_dpinfo,
        'min_dp=i' => \$min_dp,
        'ref_only!' => \$ref_only,
);
&help() if $opt_help;
&help() unless defined $vcfposs_file and defined $packpg_list_file and defined $outdir;
&help("input_format error!") unless $input_format eq 'vg' or $input_format eq 'bam';

if($input_format eq 'vg') {
    die "no vg path in config" unless defined $vg;
} elsif ($input_format eq 'bam') {
    die "no samtools path in config" unless defined $samtools;
}

$need_chr = undef if defined $need_chr and $need_chr eq 'all';
if (defined $output_dpinfo) {
    $output_dpinfo_fh = open_out_fh($output_dpinfo);
    say $output_dpinfo_fh join "\t", qw/#sid  ref_count  ref_dp  nr_count  nr_dp/;
}

my $dp_around_bp=100;

my $infiless = $input_format eq 'vg' ? &read_packpg_list_file($packpg_list_file) : &read_bam_list_file($packpg_list_file);

my $need_ids = [ sort keys %$infiless ]; ##########
$threads = scalar(@$need_ids) if $threads > scalar(@$need_ids);
say STDERR "Start to calculate depth for: " . join(' ', @$need_ids);



my $lastpos = -1;
say STDERR "Now read vcf_poss_file: $vcfposs_file";
my ($vcfposs) = &read_vcf_poss_file($vcfposs_file);
say STDERR "Done read vcf_poss_file";

if ($threads==1) {
    foreach my $sid (@$need_ids) {
        die unless exists $$infiless{$sid};
        my $infiles = $$infiless{$sid};
        my $outfile = "$outdir/$sid.depth.txt.gz";
        my $O = open_out_fh($outfile);
        say STDERR "$sid start";
        &process_dp($infiles, $sid, $O);
        say STDERR "$sid Done";
    }
} else { # multi-threads
    MCE::Loop->init(
        max_workers => $threads, chunk_size => 1
    );
    mce_loop {
        my $sid = $_;
        die unless exists $$infiless{$sid};
        my $infiles = $$infiless{$sid};
        my $outfile = "$outdir/$sid.depth.txt.gz";
        my $O = open_out_fh($outfile);
        say STDERR "$sid start";
        &process_dp($infiles, $sid, $O);
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
    my ($infiles, $sid, $O) = @_;
    my (%counts, %sum_dps);
    my $cmd;
    if($input_format eq 'vg') {
        my ($pack, $pg) = @$infiles;
        die unless -e $pack;
        die unless -e $pg;
        $cmd = "$vg depth --threads 4 --pack $pack $pg";
        if (defined $need_chr) {
            $cmd .= " --ref-path $need_chr";
        }
    } elsif ($input_format eq 'bam') {
        my ($bam) = @$infiles;
        die unless -e $bam;
        $cmd = "$samtools depth -a $bam";
        if (defined $need_chr) {
            $cmd .= " -r $need_chr";
        }
    } else {
        die "input_format error!";
    }
    #$cmd = 'zcat ~/a.gz'; ############# debug
    open(my $P, "$cmd |");
    my $is_first_line = 1;
    my %buf;
    my @buf_all_tmp;
    my %buf_all_all;
    my $chr_old = '???';
    LINEDP:while(<$P>) {
        # chr pos dp
        #  0   1   2
        chomp;
        next unless $_;
        my ($chr, $pos, $dp) = split(/\t/, $_); # 1 based
        next if $dp==0;
        if ($chr_old ne $chr) {
            if ($is_first_line==1) {
                $chr_old = $chr;
                $is_first_line=0;
                redo;
            }
            &add_dp_summary($chr_old, \%buf_all_all, \@buf_all_tmp);
            @buf_all_tmp=();
            foreach my $end (sort {$a<=>$b} keys %buf) {
                &print1($buf{$end}, $chr_old, $end, $O);
            }
            %buf = ();
            $chr_old = $chr;
            redo LINEDP;
        }
        $buf_all_tmp[0]++;
        $buf_all_tmp[1] += $dp;
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
                    $buf{$end}{$start}[0] += 1 if $dp > $min_dp;
                    $buf{$end}{$start}[1] += $dp;
                }
            } else {
                my $toprint = delete $buf{$end};
                &print1($toprint, $chr_old, $end, $O);
            }
        }
    }
    &add_dp_summary($chr_old, \%buf_all_all, \@buf_all_tmp);
    if(defined $output_dpinfo) {
        my $ref_count = $buf_all_all{ref}[0] // 0;
        my $ref_sumdp = $buf_all_all{ref}[1] // 0;
        my $nonref_count = $buf_all_all{nonref}[0] // 0;
        my $nonref_sumdp = $buf_all_all{nonref}[1] // 0;
        say $output_dpinfo_fh join "\t", $sid, $ref_count, $ref_sumdp, $nonref_count, $nonref_sumdp;
    }
}

sub add_dp_summary {
    my ($chr, $buf_all_all, $buf_all_tmp) = @_;
    if($ref_only==1 or $chr=~/^\d+$/) { # ref
        $$buf_all_all{ref}[0] += $$buf_all_tmp[0];
        $$buf_all_all{ref}[1] += $$buf_all_tmp[1];
    } else { # non-ref
        $$buf_all_all{nonref}[0] += $$buf_all_tmp[0];
        $$buf_all_all{nonref}[1] += $$buf_all_tmp[1];
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

sub read_bam_list_file {
    my ($file) = @_;
    my $I = open_in_fh($file);
    my %ret;
    while(<$I>) {
        # SampleID xx.pack xx.pg
        chomp;
        next unless $_;
        next if /^#/;
        my ($sid, $bam) = split(/\s+/, $_);
        die "$bam not exists!" unless -s $bam;
        $ret{$sid} = [$bam];
    }
    return \%ret;
}


sub help() {
    my ($info) = @_;
    if (defined $info) {
        say STDERR $info;
        say STDERR "----------------------------------------";
    }
    die "usage: perl x.pl [options]
        -v|--vcfposs <file>     input vcfposs file
        -l|--list_packpg <file> input pack_pg file
        -o|--outdir <dir>       output dir
        -t|--threads| <int>     threads number. Default: $threads
        -h|--help               help
        --chr <str>             process this chromsome only. Default: all
        --input_format <str>    vg or bam. Default: vg
        --output_avgdpinfo <file>  output avg dpinfo file. Default: not output
        --min_dp <int>          min depth to count. Default: $min_dp
        --ref_only              Treated all chromsomes as ref, not non-ref
        ";
}
