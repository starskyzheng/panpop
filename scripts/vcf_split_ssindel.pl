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
use Getopt::Long;

sub help {
    die <<EOF;
usage: perl $0 [parameters]
--invcf <file>
--outvcf_short <file>
--outvcf_long <file>
--ssindel_maxlen <int>
--ssindel_maxallele <int>
--threads <int>
EOF
}


my ($invcf, $outvcfs, $outvcfl, $opt_help);

my $SV_min_length = 50;
my $ssindel_maxallele = 3;
my $ssindel_maxlen = 5;
my $threads = 6;

GetOptions (
        'help|h!' => \$opt_help,
        't|threads=i' => \$threads,
        'ssindel_maxlen=s' => \$ssindel_maxlen,
        'ssindel_maxallele=s' => \$ssindel_maxallele,
        'SV_min_length' => \$SV_min_length,
        'i|invcf=s' => \$invcf,
        'o|outvcf_short=s' => \$outvcfs, # small indel/snp
        'O|outvcf_long=s' => \$outvcfl, # long indel/sv
);
&help() if $opt_help;
&help() unless defined $invcf;



my $I = open_in_fh($invcf);
my $OS = open_out_fh($outvcfs, $threads); # small indel/snp
my $OL = open_out_fh($outvcfl, $threads); # long indel/sv

while(my $line = <$I>) {
    chomp $line;
    next unless $line;
    if ($line=~/^#/) {
        say $OS $line;
        say $OL $line;
        next;
    }
    my $max_len = -1;
    $line=~m/^\S+\t\S+\t\S+\t(\S+)\t(\S+)\t/ or die $line;
    my ($ref, $alts) = ($1, $2);
    $max_len = length($ref);
    my $allele_count=0;
    foreach my $alt (split(/,/, $alts)) {
        $allele_count++;
        my $len = length($alt);
        $max_len = $len if $max_len < $len;
    }
    if ($max_len < 0) {
        die;
    } elsif ( $max_len <= $ssindel_maxlen ) {
        say $OS $line;
    } elsif ( $max_len < $SV_min_length and 
                $allele_count <= $ssindel_maxallele )  {
        say $OS $line;
    } else {
        say $OL $line;
    }
}





