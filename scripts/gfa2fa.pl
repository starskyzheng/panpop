#!/usr/bin/env perl
#===============================================================================
#
#         FILE: gfa2fa.pl
#
#        USAGE: ./gfa2fa.pl in.gfa out.fa
#
#  DESCRIPTION: 根据path信息，将gfa转换为fasta
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zeyu Zheng (LZU), zhengzy2014@lzu.edu.cn
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 11/12/2021 02:01:55 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
no warnings qw( experimental::smartmatch );
use v5.10;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use zzIO;

sub help {
    die "usage: perl xx.pl in.gfa out.fa";
}

my ($in, $out) = @ARGV;

$out//='-';
&help() unless $in;

my $O = open_out_fh($out);
my ($segs2seq, $path2links) = &read_gfa($in);

foreach my $path (sort keys %$path2links) {
    my $final_seq;
    my $links = $$path2links{$path};
    my @links = split(/,/, $links);
    foreach my $seg_phase (@links) {
        $seg_phase=~/^([a-zA-Z0-9]+)([+-])$/ or die "$seg_phase error!";
        my ($seg, $phase) = ($1, $2);
        if ($phase eq '+') {
            $final_seq .= $$segs2seq{$seg};
        } else {
            die;
        }
    }
    say $O ">$path";
    say $O $final_seq;
}

sub read_gfa {
    my ($in) = @_;
    my $I = open_in_fh($in);
    my (%segs, %paths);
    while(<$I>) {
        chomp;
        next unless $_;
        my @F = split(/\t/);
        if($F[0] eq 'S') {
            $segs{$F[1]} = $F[2];
        } elsif ($F[0] eq 'P') {
            $paths{$F[1]} = $F[2];
        } elsif ($F[0] eq 'L') {
            #
        } else {
            # headers
        }
    }
    return(\%segs, \%paths);
}

