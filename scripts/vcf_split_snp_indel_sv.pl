#!/usr/bin/env perl
#===============================================================================
#
#         FILE: vcf2snp_indel_sv.pl
#
#        USAGE: ./vcf2snp_indel_sv.pl  
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
#      CREATED: 11/21/2021 11:01:43 PM
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

sub help {
    die <<EOF;
usage: perl $0 IN.vcf [OUT_PREFIX] [SV min length (default 50)]
EOF
}

my $ARGVs = join(' ',@ARGV);
my ($in, $out_prefix, $min_sv) = @ARGV;

&help() unless $in;

my $threads = 12;
$out_prefix //= $in;
$min_sv //= 50;


my $I = open_in_fh($in);
my $O_SNP = open_out_fh("$out_prefix.snp.vcf.gz", $threads);
my $O_INDEL = open_out_fh("$out_prefix.indel.vcf.gz", $threads);
my $O_SV = open_out_fh("$out_prefix.sv.vcf.gz", $threads);

while(my $line = <$I>) {
    chomp $line;
    next unless $line;
    if ($line=~/^#/) {
        if ($line=~/^#CHROM/) {
            say $O_SNP qq(##CommandLine="$0 $ARGVs");
            say $O_INDEL qq(##CommandLine="$0 $ARGVs");
            say $O_SV qq(##CommandLine="$0 $ARGVs");
        }
        say $O_SNP $line;
        say $O_INDEL $line;
        say $O_SV $line;
        next;
    }
    my $max_len = -1;
    $line=~m/^\S+\t\S+\t\S+\t(\S+)\t(\S+)\t/ or die $line;
    my ($ref, $alts) = ($1, $2);
    $max_len = length($ref);
    foreach my $alt (split(/,/, $alts)) {
        my $len = length($alt);
        $max_len = $len if $max_len < $len;
    }
    if ($max_len < 0) {
        die;
    } elsif ($max_len == 1) {
        say $O_SNP $line;
    } elsif ($max_len >= $min_sv) {
        say $O_SV $line;
    } else {
        say $O_INDEL $line;
    }
}

