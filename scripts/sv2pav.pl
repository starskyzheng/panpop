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
#      CREATED: 11/24/2021 11:54:32 AM
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
use Getopt::Long;
use List::Util qw/max min/;
use MCE::Flow;
use MCE::Candy;

my ($in, $out, $opt_help);# = @ARGV;
my $thread = 1;

sub usage {
    print <<EOF;
Usage: $0 [options]
Options:
    -i | --invcf	    input vcf file
    -o | --outvcf	    output file, default: -
    -p | --threads  	thread number, default: 1
EOF
    exit;
}

GetOptions (
    'help|h!' => \$opt_help,
    'i|invcf=s' => \$in,
    'o|outvcf=s' => \$out,
    't|threads=i' => \$thread,
);

&usage() if $opt_help or not $in or not $out;

my $I = open_in_fh($in);
my $O = open_out_fh($out);


if ($thread eq 1) {
    while(<$I>) {
        my $result = &process_line($_);
        say $O $result if defined $result;
    }
    exit();
}

sub flt_mce {
    my ( $mce, $chunk_ref, $chunk_id ) = @_;
    my $result = &process_line($$chunk_ref[0]);
    if (defined $result) {
        $mce->gather($chunk_id, $result, "\n") 
    } else {
        $mce->gather($chunk_id);
    }
}

mce_flow_f {
    chunk_size => 1,
    max_workers => $thread,  
    gather => MCE::Candy::out_iter_fh($O)
    },
    \&flt_mce,
    $I;
close $I;
exit 0;

sub process_line {
    my ($line) = @_;
    chomp $line;
    return undef unless $line;
    if ($line=~/^#/) {
        return $line;
    }
    my @F = split /\t/, $line;
    my $ref = $F[3];
    my @alts = split /,/, $F[4];
    my @ref_alts = ($ref, @alts);
    my $reflen = length $ref;
    my $alt_max_len = max map {length} @alts;
    my $alt_min_len = min map {length} @alts;
    my %replase;
    if ($reflen <= 2 and $alt_max_len > 10) { # ins
        foreach my $ialt (1..$#ref_alts) {
            my $len = length $ref_alts[$ialt];
            $replase{$ialt} = 0 if $len <= 2;
        }
    } elsif ($reflen > 10 and $alt_min_len <= 2) { # del
        my $shortest_ialt = 0;
        my $shortest_len = length $alts[0];
        foreach my $ialt (1..$#ref_alts) {
            my $len = length $ref_alts[$ialt];
            if ($len < $shortest_len) {
                $shortest_ialt = $ialt;
                $shortest_len = $len;
            }
        }
        foreach my $ialt (0..$#ref_alts) {
            my $len = length $ref_alts[$ialt];
            $replase{$ialt} = $shortest_ialt if $len <= 2;
        }
    }
    #say STDERR Dumper \%replase;
    foreach my $isid (9..$#F) {
        my $alle_infos = $F[$isid];
        $alle_infos =~ /^([^:]+)(:|$)/ or die;
        my $GT = $1;
        my @alleles = split /[\|\/]/, $GT;
        die unless scalar(@alleles)==2;
        foreach my $ialt (0..$#alleles) {
            my $alle = $alleles[$ialt];
            if (exists $replase{$alle}) {
                $alleles[$ialt] = $replase{$alle};
            }
        }
        $F[$isid] = join "/", @alleles;
    }
    return join("\t", @F);
}


