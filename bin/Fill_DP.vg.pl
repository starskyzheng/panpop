#!/usr/bin/env perl
#===============================================================================
#
#         FILE: Fill_DP.pl
#
#        USAGE: ./Fill_DP.pl  
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
#      CREATED: 10/17/2023 11:29:22 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
no warnings qw( experimental::smartmatch );
use v5.24;

use FindBin qw($Bin);
use lib "$Bin/../lib";

#use Data::Dumper;
use Getopt::Long;
use MCE::Loop;
use zzIO;

use read_config qw/read_config_yaml/;

my ($invcf_ori, $outdir, $ref_fasta_file, $sid2files_file);
my ($force_run, $opt_help);
my $threads = 32;

sub usage {
    my $msg = shift;
    print STDERR "$msg\n" if $msg;
    print STDERR <<"EOF";
USAGE: perl $0 [options]
options:
    -i, --in_vcf <FILE>    Input vcf file
    -o, --outdir <DIR>     output dir, should not exists or empty
    -r, --ref_fasta_file <FILE>      reference fasta file
    -t, --threads <int>    threads, default: $threads
    --force                force run, skip check outdir is empty.
    --sid2files <FILE>     sid2files file (sid pack vg vcf)
    -h, --help             help, print this message
EOF
    die;
}

&usage() if $opt_help;


GetOptions (
    'h|help!' => \$opt_help,
    'i|in_vcf=s' => \$invcf_ori,
    'o|outdir=s' => \$outdir,
    'b|sid2files=s' => \$sid2files_file,
    'r|ref_fasta_file=s' => \$ref_fasta_file,
    'force' => \$force_run,
    't|threads=i' => \$threads,
);

&usage() if !$invcf_ori or !$outdir or !$ref_fasta_file;

&chech_outdir_is_empty($outdir, $force_run);

our $config = read_config_yaml("$Bin/../config.yaml");
my $scripts_dir = "$Bin/../scripts";
my $sid2bam = read_sid2file_list($sid2files_file);

my $sids = read_vcf_header($invcf_ori, $sid2bam);

### 1.vcf2poss
{
    my $invcf = $invcf_ori;
    my $outpos = "$outdir/1.poss";
    my $log = "$outpos.log";
    my $cmd = "perl $scripts_dir/vcf2pos_range.pl $invcf $outpos 2>&1 | tee $log";
    &zzsystem($cmd);
}

### 2.cal_aug_dp
{
    my $outdir_now = "$outdir/2.cal_aug_dp";
    my $invcf = $invcf_ori;
    my $ref = $ref_fasta_file;
    my $outpos = "$outdir/1.poss";
    mkdir $outdir_now unless -e $outdir_now;
    MCE::Loop->init( max_workers => $threads, chunk_size => 1 );
    #foreach (@$sids) {
    mce_loop {
        my $sid = $_;
        my $pack = $$sid2bam{$sid}[0];
	my $vg = $$sid2bam{$sid}[1];
        my $dpfile = "$outdir_now/$sid.depth.txt.gz";
        my $list_packpg = "$outdir_now/$sid.list_packpg";
        my $log = "$dpfile.log";
        `echo '$sid $pack $vg' > $list_packpg`;
        my $cmd = "perl $scripts_dir/cal_range_depth_aug.pl --input_format vg --vcfposs $outpos --list_packpg $list_packpg --outdir $outdir_now -t 1 2>&1 | tee $log";
        &zzsystem($cmd);
    } (@$sids);
    MCE::Loop->finish;
}


### 3.fill_aug_dp
{
    my $outdir_now = "$outdir/3.fill_dp";
    my $outpos = "$outdir/1.poss";
    my $ref = $ref_fasta_file;
    mkdir $outdir_now unless -e $outdir_now;
    MCE::Loop->init( max_workers => $threads, chunk_size => 1 );
    #foreach (@$sids) {
    mce_loop {
        my $sid = $_;
        my $vcf_sample_in = $$sid2bam{$sid}[2];
        my $vcf_sample_out = "$outdir_now/$sid.vcf.gz";
        my $vcf_sample_out_sorted = "$outdir_now/$sid.sort.vcf.gz";
        my $dpfile = "$outdir/2.cal_aug_dp/$sid.depth.txt.gz";
        my $list_packpg = "$outdir_now/$sid.list_packpg";
        my $min_dp = 0;
        my $min_cov = $$config{'aug_nonmut_min_cov'};
        my $log = "$dpfile.log";
        my $cmd = "perl $scripts_dir/cal_range_depth_aug.fillvcf.pl --in_vcf $vcf_sample_in --out_vcf $vcf_sample_out --ref $ref --min_cov $min_cov --min_dp $min_dp --dp_file $dpfile 2>&1 | tee $log";
        &zzsystem($cmd);
        my $cmd2 = "bcftools sort $vcf_sample_out -o $vcf_sample_out_sorted; tabix $vcf_sample_out_sorted";
        &zzsystem($cmd2);
    }(@$sids);
    MCE::Loop->finish;
}


### 4.bcftools merge
{
    my $outdir_now = "$outdir/4.bcftools_merge";
    my $listfile = "$outdir_now/1.list.txt";
    my $outvcf = "$outdir_now/2.merged.vcf.gz";
    my $log = "$outvcf.log";
    mkdir $outdir_now unless -e $outdir_now;
    my $O = open_out_fh($listfile);
    foreach my $sid (@$sids) {
        my $vcf_sample_in = "$outdir/3.fill_dp/$sid.sort.vcf.gz";
        say $O $vcf_sample_in;
    }
    close $O; undef $O;
    my $cmd = "bcftools merge -l $listfile -m none 2>> $log | bcftools sort 2>>$log | bgzip -c > $outvcf 2>>$log";
    &zzsystem($cmd);
}


exit;


sub zzsystem {
    my ($cmd) = @_;
    say STDERR $cmd;
    system($cmd);
}

sub read_vcf_header {
    my ($invcf, $sid2bam) = @_;
    my $I = open_in_fh($invcf);
    my @header;
    while(<$I>) {
        chomp;
        next if /^##/;
        if(/^#/) {
            @header = split /\t/, $_;
            last;
        } else {
            die;
        }
    }
    my @sids = @header[9..$#header];
    foreach my $sid (@sids) {
        if(!exists $$sid2bam{$sid}) {
            die "sid $sid not exists in bamlist file";
        }
    }

    return \@sids;
}

sub read_sid2file_list {
    my ($infile) = @_;
    my %bams;
    my $I = open_in_fh($infile);
    while(<$I>) {
        chomp;
        my ($sid, $pack, $vg, $vcf) = split /\s+/;
        $bams{$sid} = [$pack, $vg, $vcf];
    }
    return \%bams;
}


sub chech_outdir_is_empty {
    my ($dir, $force_run) = @_;
    if(! -e $dir) {
        mkdir $dir;
        return;
    } else {
        my @files = <$dir/*>;
        if(@files) {
            say STDERR "outdir ` $dir ` is not empty, please check it\n";
            if($force_run) {
                #say STDERR "force run, delete all files in $dir\n";
                #system("rm -rf $dir/*");
            } else {
                exit(1);
            }
        }
    }
}


