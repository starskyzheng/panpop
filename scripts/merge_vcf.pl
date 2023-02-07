#!/usr/bin/env perl
#===============================================================================
#
#         FILE: merge_vcf.pl
#
#        USAGE: ./merge_vcf.pl
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
#      CREATED: 12/31/2021 05:41:54 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Long;
use v5.10;


use FindBin qw($Bin);
use lib "$Bin/../lib";

use Carp qw/confess carp/; # carp=warn;confess=die
use Data::Dumper;
use zzIO;
use MCE::Loop;
use File::Temp;
use Getopt::Long;
use read_config qw/read_config_yaml/;


sub help {
    print STDERR <<EOF;
usage: perl xx.pl 
    --inlist FILE.vcfs.list 
    --out FILE.vcf.gz
    --tmp_dir PATH 
    [--chrome STRING] Chromosome name
    [--threads INT (default: 1)]  Number of bcftools running together
    [--vcfs_per_run INT (default: 999999)]
    [--bcftools_threads INT (default: 4)]  Number of threads of 'bcftools merge'
    [--noclean]       Do not clean tmp files
    [--invcf_nocheck] Do not check input vcf.gz files wather they are exists and has .tbi files.
    [--help]          Print this help message
EOF
    say STDERR @_ if @_;
    exit(-1);
}

my ($inlist, $out, $tmp_dir_top, $chrome);

my $opt_help = 0;
my $invcf_nocheck = 0;
my $vcfs_per_run = 999999;
my $threads = 1;
my $noclean = 0;
my $bcftools_threads = 4;
my $m_none=1;
my $bcftools;

GetOptions (
        'help|h!' => \$opt_help,
        'inlist|i=s' => \$inlist,
        'out|o=s' => \$out,
        'threads|t=i' => \$threads,
        'tmp_dir|T=s' => \$tmp_dir_top,
        'vcfs_per_run|V=i' => \$vcfs_per_run,
        'invcf_nocheck!' => \$invcf_nocheck,
        'noclean!' => \$noclean,
        'bcftools_threads=i' => \$bcftools_threads,
        'chrome|r=s' => \$chrome,
        'bcftools_bin=s' => \$bcftools,
        'm_none=i' => \$m_none,
);

&help() if $opt_help;
&help() unless($inlist && $out && $tmp_dir_top);

my $bcftools_appends="--threads $bcftools_threads";
$bcftools_appends .= " -m none --non_normalize_alleles" if $m_none==1;
$bcftools_appends .= " --regions $chrome" if defined $chrome;

my $config = read_config_yaml("$Bin/../config.yaml");
my $tabix = $$config{tabix};
$bcftools //= $$config{bcftools};

my $tmp_dir_ = &get_tmp_dir($tmp_dir_top, $noclean);
my $tmp_dir = $$tmp_dir_;
my $vcfs = &read_list($inlist);
my $lists = &split_filelist($vcfs, $vcfs_per_run);
my $lists_files = &write_splitlist($lists, $tmp_dir);
say STDERR "Total ", scalar(@$lists), " jobs for Step 1.";
say STDERR "Now start to merge raw vcf files... (Step 1/2)";
my $vcfs1 = &merge_vcfs1($lists_files);
say STDERR "Now start to merge merged files... (Step 2/2)";
&merge_vcfs2($vcfs1, $out);
say STDERR "Done!";

exit;


sub merge_vcfs2 {
    my ($vcfs, $final_vcf) = @_;
    my $vcfs_count = scalar(@$vcfs);
    if($vcfs_count == 1) {
        my $vcf = $$vcfs[0];
        my $cmd = "cp $vcf $final_vcf";
        system($cmd);
        say STDERR "Only one vcf file, so just copy it to $final_vcf";
        return $final_vcf;
    } else {
        my $vcfs_str = join(" ", @$vcfs);
        my $log = "$final_vcf.log";
        my $cmd = "$bcftools merge $bcftools_appends --output-type z --output $final_vcf $vcfs_str";
        say STDERR $cmd;
        system($cmd);
        return $final_vcf;
    }
}



sub merge_vcfs1 {
    my ($lists_files) = @_;
    my @ret_vcfs;
    my @cmds;
    foreach my $list_file (@$lists_files) {
        my $outvcf = "$list_file.vcf.gz";
        my $log = "$list_file.vcf.gz.log";
        my $cmd = "$bcftools merge $bcftools_appends --output-type z --output $outvcf --file-list $list_file; $tabix $outvcf";
        push @ret_vcfs, $outvcf;
        push @cmds, $cmd;
    }
    &run_cmd(\@cmds, $threads);
    return \@ret_vcfs;
}

sub run_cmd {
    my ($cmds, $threads) = @_;
    MCE::Loop->init( max_workers => $threads, use_slurpio => 1, chunk_size => 1 );
    mce_loop {
        my $cmd = $_;
        system($cmd);
    } @$cmds;
    MCE::Loop->finish;
}



sub get_tmp_dir {
    my ($dir, $noclean) = @_;
    -e $dir or mkdir $dir or confess "cannot mkdir $dir: $!";
    my $unlink=1;
    if($noclean==1) {
        $unlink = 0;
        #$File::Temp::KEEP_ALL = 1;
    }
    my $ret = File::Temp->newdir( DIR => $dir, UNLINK => $unlink );
    $ret->unlink_on_destroy(0) if $noclean==1;
    say STDERR "tmp_dir: $ret";
    return \$ret;
}



sub write_splitlist {
    my ($lists, $outdir) = @_;
    my @listfiles;
    my $i = 0;
    for my $list (@$lists) {
        my $out = "$outdir/1.$i.list";
        open my $fh, ">$out" or confess $!;
        for my $file (@$list) {
            say $fh $file;
        }
        push @listfiles, $out;
        close $fh;
        $i++;
    }
    return \@listfiles;
}

sub split_filelist {
    my ($filelist, $num) = @_;
    my @filelist = @$filelist;
    my $max_i = scalar(@filelist) - 1;
    my @split_filelist;
    my $i = 0;
    while ($i < $max_i) {
        my $to = $i + $num;
        $to = $max_i+1 if $to > $max_i;
        my @tmp = @filelist[$i .. $to-1];
        push @split_filelist, \@tmp;
        #say STDERR join " ", $i, $to-1;
        $i += $num;
    }
    return \@split_filelist;
}


sub read_list {
    my $list = shift;
    my @files;
    open my $fh, '<', $list or confess $!;
    while (<$fh>) {
        chomp;
        next unless $_;
        if ($invcf_nocheck!=1) {
            confess "file not exists: $_" unless -e $_;
            my $tbi = "$_.tbi";
            confess "tbi not exists: $tbi" unless -e "$tbi";
        }
        push @files, $_;
    }
    close $fh;
    return \@files;
}



