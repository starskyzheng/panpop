#!/usr/bin/env perl
#===============================================================================
#
#         FILE: PART_run.pl
#
#        USAGE: ./PART_run.pl  
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
#      CREATED: 10/17/2023 09:54:28 PM
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

use read_config qw/read_config_yaml/;

my ($invcf_ori, $outdir, $ref_fasta_file);
my ($force_run, $opt_help);
my $threads = 32;
my $not_first_merge = 0;

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
    --not_first_merge      Not first merge.
    --force                force run, skip check outdir is empty.
    -h, --help             help, print this message
EOF
    die;
}

&usage() if $opt_help;


GetOptions (
    'h|help!' => \$opt_help,
    'i|in_vcf=s' => \$invcf_ori,
    'o|outdir=s' => \$outdir,
    'r|ref_fasta_file=s' => \$ref_fasta_file,
    'force' => \$force_run,
    't|threads=i' => \$threads,
    'not_first_merge' => \$not_first_merge,
);

&usage() if !$invcf_ori or !$outdir or !$ref_fasta_file;

&chech_outdir_is_empty($outdir, $force_run);

our $config = read_config_yaml("$Bin/../config.yaml");
my $scripts_dir = "$Bin/../scripts";


#### 1. realign0
{
    my $invcf = $invcf_ori;
    my $outvcf = "$outdir/1.realign0.unsorted.vcf.gz";
    my $outvcf_sorted = "$outdir/1.realign0.sorted.vcf.gz";
    my $log = "$outvcf.log";
    my $append = '';
    $append .= '--first_merge' if $not_first_merge==0;
    my $cmd = "perl $scripts_dir/realign.pl --chr_tolerance --in_vcf $invcf --out_vcf $outvcf --ref_fasta_file $ref_fasta_file --threads $threads --ext_bp_max 500 --ext_bp_min 50 --skip_mut_at_same_pos 2 --level 1 $append 2>&1 | tee $log";
    zzsystem($cmd);
    &sort_vcf($outvcf, $outvcf_sorted);
}

#### 2.thin1
{
    my $invcf = "$outdir/1.realign0.sorted.vcf.gz";
    my $outvcf = "$outdir/2.thin1.unsorted.vcf.gz";
    my $outvcf_sorted = "$outdir/2.thin1.sorted.vcf.gz";
    my $log = "$outvcf.log";
    my $sv2pav_merge_identity_threshold = $config->{sv2pav_merge_identity_threshold};
    my $sv2pav_merge_diff_threshold = $config->{sv2pav_merge_diff_threshold};
    my $cmd = "perl $scripts_dir/merge_similar_allele.pl --type 3 --invcf $invcf --outvcf $outvcf --sv2pav_merge_identity_threshold $sv2pav_merge_identity_threshold --sv2pav_merge_diff_threshold $sv2pav_merge_diff_threshold --threads $threads 2>&1 | tee $log";
    zzsystem($cmd);
    &sort_vcf($outvcf, $outvcf_sorted);
}


#### 2.thin2
{
    my $invcf = "$outdir/2.thin1.sorted.vcf.gz";
    my $outvcf = "$outdir/2.thin2.unsorted.vcf.gz";
    my $outvcf_sorted = "$outdir/2.thin2.sorted.vcf.gz";
    my $log = "$outvcf.log";
    my $cmd = "perl $scripts_dir/sv2pav.pl --invcf $invcf --outvcf $outvcf --enable_norm_alle 1 --max_len_tomerge 20 --sv_min_dp 40 --threads $threads 2>&1 | tee $log";
    zzsystem($cmd);
    &sort_vcf($outvcf, $outvcf_sorted);
}

#### final
{
    my $invcf = "2.thin2.sorted.vcf.gz";
    my $outvcf = "$outdir/3.final.vcf.gz";
    my $cmd = "ln -s $invcf $outvcf";
    zzsystem($cmd);
}

exit;

sub zzsystem {
    my ($cmd) = @_;
    say STDERR $cmd;
    system($cmd);
}

sub sort_vcf {
    my ($invcf, $outvcf) = @_;
    my $cmd = "bcftools sort $invcf -o $outvcf";
    zzsystem($cmd);
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
