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
use List::Util qw/min max/;

#our $config = read_config_yaml("$Bin/../config.yaml");

my $rangebp = 5;
my ($in1, $in2) = @ARGV;
#   long  short

my $outfile = '-';
my $O = open_out_fh($outfile);

say STDERR "Now read sort vcf: $in2";
my ($h2, $v2, $ranges) = &read_vcf(open_in_fh($in2), 0.01, 1e6);

my $I = open_in_fh($in1);


my $all=0;
my $haspos = 0;


say $O join "\t", qw/ chr pos haspos totalid sameid notsameid 
                      in2miss1 in1miss2 missboth dist
                      reflen altslenmin altslenmax is_snp/;

my $totalid;
my $ids_union;
my %ranget;

my @header;
say STDERR "Now process long vcf: $in1";
while(my $line = <$I>) {
    #last if $.>10000;
    chomp $line;
    if ($line=~/^#/) {
        next if $line=~/^##/;
        @header = split(/\t/, $line);
        next;
    }
    $line=~s/^pchr(0)?//;
    $line=~/^(\S+)\t(\S+)/ or die;
    my ($chr, $pos) = ($1, $2);
    if ( ! exists $$v2{$chr}{$pos} and ! exists $$ranges{$chr}{$pos} ) {
        next;
    } elsif ( ! exists $$v2{$chr}{$pos} and exists $$ranges{$chr}{$pos} ) {
        my @F = split(/\t/, $line);
        my $reflen = length($F[3]);
        my ($altslenmin, $altslenmax, $is_snp) = &cal_alts_len($F[4]);
        $is_snp=0 if $reflen!=1;
        my $svpos = $$ranges{$chr}{$pos};
        if ($is_snp==0) {
            $ranget{$chr}{$svpos} = \@F;
        }
        next;
    } elsif( exists $$v2{$chr}{$pos} ) {
        my @F = split(/\t/, $line);
        my $print = &cal($v2, \@F, \%ranget, $ranges);
        say $O $print;
    } else {
        die;
    }
}

# snp with sv/indel around
foreach my $chr (sort keys %ranget) {
    foreach my $pos (sort {$a<=>$b} keys $ranget{$chr}->%*) {
        my $F = $ranget{$chr}{$pos};
        #my $pos = $$ranges{$chr}{$i} // die;
        #my $print = &cal( $v2, $F );
        #say $O $print;
        my $reflen = length($$F[3]);
        my ($altslenmin, $altslenmax, $is_snp) = &cal_alts_len($$F[4]);
        $is_snp=0 if $reflen!=1;
        my $miss=0;
        die "$chr $pos" unless exists $$v2{$chr}{$pos};
        foreach my $id (keys %$ids_union) {
            if ($$v2{$chr}{$pos}{$id} eq '-') {
                $miss++;
            }
        }
        delete $$v2{$chr}{$pos};
        say $O join "\t", $chr, $pos, 2, $totalid, -1, -1,
                            -1, -1, $miss, abs($$F[1]-$pos), 
                            $reflen, $altslenmin, $altslenmax, $is_snp;
    }
}

# print miss
foreach my $chr (sort keys %$v2) {
    foreach my $pos (sort {$a<=>$b} keys $$v2{$chr}->%*) {
        my $miss=0;
        die unless exists $$v2{$chr}{$pos};
        foreach my $id (keys %$ids_union) {
            if ($$v2{$chr}{$pos}{$id} eq '-') {
                $miss++;
            }
        }
        say $O join "\t", $chr, $pos, 0, $totalid, -1, -1,
                                -1, -1, $miss, -1,
                                -1, -1, -1, -1;
    }
}

exit;

sub cal {
    my ($v2, $F, $ranget, $ranges) = @_;
    my $chr = $$F[0];
    my $pos = $$F[1];
    # remove ranges and ranget
    for(my $i = $pos-$rangebp; $i<=$rangebp+$pos; $i++) {
        delete $$ranget{$chr}{$i} if $ranget;
        delete $$ranges{$chr}{$i} if $ranges;
    }
    my %vcf;
    my @ref_alts = ($$F[3], split(/,/, $$F[4]));
    #@ref_alts = (0, (1) x split(/,/, $F[4]) ); ############
    my $maxi = scalar(@$F)-1;
    foreach my $i (9..$maxi) {
        my $infos = $$F[$i];
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
    my $reflen = length($$F[3]);
    my ($altslenmin, $altslenmax, $is_snp) = &cal_alts_len($$F[4]);
    $is_snp=0 if $reflen!=1;
    my $ret = join "\t", $chr, $pos, 1, $totalid, $sameid, $notsameid, 
                            $in2miss1, $in1miss2, $missboth, 0, 
                            $reflen, $altslenmin, $altslenmax, $is_snp;
    return($ret);
}


sub cal_alts_len {
    my ($alts) = @_;
    my @altslen;
    my $is_snp=1;
    foreach my $alt (split(/,/, $alts)) {
        push @altslen, length($alt);
        $is_snp=0 if $alt eq '*';
        $is_snp=0 if length($alt)>=2;
    }
    my $altslenmin = min(@altslen);
    my $altslenmax = max(@altslen);
    return($altslenmin, $altslenmax, $is_snp);
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
    my %ranges;
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
        for(my $i = $pos-$rangebp; $i<=$rangebp+$pos; $i++) {
            if ( exists $ranges{$chr}{$i} ) { # remove oldpos
                my $oldpos = $ranges{$chr}{$i};
                delete $ranges{$chr}{$_} foreach $oldpos-$rangebp .. $oldpos+$rangebp;
                delete $vcf{$chr}{$i};
                #die "$chr $pos: $i exists!" 
            }
            $ranges{$chr}{$i} = $pos;
        }
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
    return(\@header, \%vcf, \%ranges);
}

