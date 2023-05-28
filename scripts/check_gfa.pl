#!/usr/bin/env perl
use warnings;
use strict;
use v5.10;


use FindBin qw($Bin);
use lib "$Bin/../lib";
use zzIO;

my ($in, $out) = @ARGV;

$in //= '-';
$out //= '-';

my $I = open_in_fh($in);
my $O = open_out_fh($out);

while(my $line = <$I>) {
    # gfa1
    my @chrs;
    if($line=~/^P/) {
        $line=~/^P\t(\S+)/ or die;
        my $chr = $1;
        $chr=~/^\d$/ or die;
        push @chrs, $chr;
    }
}

die "No chromosome found in $in" unless @chrs;

my $chr_count = scalar @chrs;
say STDERR "OK: $chr_count chr found in $in";
