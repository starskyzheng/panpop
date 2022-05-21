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
use MCE::Flow;
use MCE::Candy;
use Getopt::Long;
use File::Basename;


#our $config = read_config_yaml("$Bin/../config.yaml");

my ($in1, $in2) = @ARGV;
#    pan   gatk

my $outfile = '-';
my $O = open_out_fh($outfile);


my ($h2, $v2) = &read_vcf(open_in_fh($in2), 0.01, 1e6);
my $I = open_in_fh($in1);


my $all=0;
my $haspos = 0;


say $O join "\t", qw/ chr pos haspos totalid sameid notsameid in2miss1 in1miss2 missboth/;

my $totalid;
my $ids_union;

my @header;
while(<$I>) {
    chomp;
    if (/^#/) {
        next if /^##/;
        @header = split(/\t/);
        next;
    }
    /^(\S+)\t(\S+)/ or die;
    my ($chr, $pos) = ($1, $2);
    if (! exists $$v2{$chr}{$pos}) {
        next;
    }
    my @F = split(/\t/);
    my %vcf;
    my @ref_alts = ($F[3], split(/,/, $F[4]));
    #@ref_alts = (0, (1) x split(/,/, $F[4]) ); ############
    foreach my $i (9..$#F) {
        my $infos = $F[$i];
        my $id = $header[$i];
        if ($infos=~m#^(\d)[/|](\d+)#) {
            my $a1 = $ref_alts[$1] // die "$1 @ref_alts";
            my $a2 = $ref_alts[$2] // die "$1 @ref_alts";
            $vcf{$id} = "$a1:$a2";
        } else {
            $vcf{$id} = '-';
        }
    }
    my $sameid = 0;
    my $notsameid = 0;
    my $in2miss1 = 0;
    my $in1miss2 = 0;
    my $missboth = 0;
    $ids_union //= &find_union_ids(\@header, $h2);
    $totalid //= scalar(keys %$ids_union);
    #foreach my $id (sort keys $$v2{$chr}{$pos}->%*) {
    foreach my $id (keys %$ids_union) {
        next unless exists $$ids_union{$id}; # pan must have this sample
        my $genotype1 = $vcf{$id} // die;
        my $genotype2 = $$v2{$chr}{$pos}{$id} // die;
        if ($genotype1 eq '-' and $genotype2 ne '-' ) {
            $in2miss1++;
        } elsif ($genotype1 ne '-' and $genotype2 eq '-' ) {
            $in1miss2++;
        } elsif ($genotype1 eq '-' and $genotype2 eq '-' ) {
            $missboth++;
        } else {
            $genotype2=~m/^(.+):(.+)$/ or die $genotype2;
            my ($a1_v2, $a2_v2) = sort($1, $2);
            $genotype1=~m/^(.+):(.+)$/ or die;
            my ($a1_v1, $a2_v1) = sort($1, $2);
            if ($a1_v1 eq $a1_v2 and $a2_v1 eq $a2_v2) {
                $sameid++;
            } else {
                $notsameid++;
            }
        }
    }
    delete $$v2{$chr}{$pos};
    say $O join "\t", $chr, $pos, 1, $totalid, $sameid, $notsameid, $in2miss1, $in1miss2, $missboth;
}

foreach my $chr (sort keys %$v2) {
    foreach my $pos (sort {$a<=>$b} keys $$v2{$chr}->%*) {
        my $miss=0;
        foreach my $id (keys %$ids_union) {
            if ($$v2{$chr}{$pos}{$id} eq '-') {
                $miss++;
            }
        }
        say $O join "\t", $chr, $pos, 0, $totalid, -1, -1, -1, -1, $miss;
    }
}


sub find_union_ids {
    my ($h1, $h2) = @_;
    my %ret;
    my $lh1 = scalar(@$h1)-1;
    my $lh2 = scalar(@$h2)-1;
    foreach my $item ($h1->@[9..$lh1]) {
        $ret{$item}++ if $item~~[$h2->@[9..$lh2]];
    }
    return \%ret;
}


sub read_vcf {
    my ($I, $samplerate, $max_sample) = @_;
    my @header;
    my %vcf;
    my $npos=0;
    while(<$I>) {
        chomp;
        if (/^#/) {
            next if /^##/;
            @header = split(/\t/);
            next;
        }
        if (defined $samplerate) {
            last if $npos > $max_sample;
            my $int=rand(1);
            next unless $int >= $samplerate;
        }
        my @F = split(/\t/);
        my ($chr, $pos) = ($F[0], $F[1]);
        my @ref_alts = ($F[3], split(/,/, $F[4]));
        #@ref_alts = (0, (1) x split(/,/, $F[4]) ); ############
        foreach my $i (9..$#F) {
            my $infos = $F[$i];
            my $id = $header[$i];
            if ($infos=~m#^(\d)[/|](\d+)#) {
                my $a1 = $ref_alts[$1] // die "$1 @ref_alts";
                my $a2 = $ref_alts[$2] // die "$1 @ref_alts";
                $vcf{$chr}{$pos}{$id} = "$a1:$a2";
            } else {
                $vcf{$chr}{$pos}{$id} = '-';
            }
        }
        $npos++;
    }
    return(\@header, \%vcf);
}

