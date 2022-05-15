#!/usr/bin/env perl
use warnings;
use v5.24;
use strict;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use zzIO;
use Carp qw/confess/;
no warnings qw( experimental::smartmatch );
#use Data::Dumper;
#use Getopt::Long;

my ($vcf, $out) = @ARGV;

my $I = open_in_fh($vcf);
my $O = open_out_fh($out);

my $chr_old;
my $pos_old;
my %lens;
while (<$I>) {
    chomp;
    next if /^#/;
    #my @F = split(/\t/);
    /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t/ or die;
    my ($chr, $pos, undef, $ref) = ($1, $2, $3, $4);
    my $ref_len = length($ref);
    $chr_old //= $chr; # init
    $pos_old //= $pos; # init
    if ($chr ne $chr_old or $pos ne $pos_old) {
        say $O join("\t", $chr_old, $pos_old, sort {$a<=>$b} keys %lens);
        $chr_old = $chr;
        $pos_old = $pos;
        %lens=();
        redo;
    }
    else {
        $lens{$ref_len}++;
    }
}

