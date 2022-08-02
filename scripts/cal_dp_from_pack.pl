#!/usr/bin/env perl
use v5.10;
use warnings;
use strict;
use File::Basename;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use read_config qw/read_config_yaml/;
use zzIO;

my $config = read_config_yaml("$Bin/../config.yaml");
my $vg = $$config{vg};

my ($pack, $xg, $outfile, $logfile, $sid, $merge_ref_nonref) = @ARGV;
&cal($pack, $xg, $outfile, $logfile, $sid, $merge_ref_nonref);
exit;

sub cal {
    my ($pack, $xg, $outfile, $logfile, $sid, $merge_ref_nonref) = @_;
    $merge_ref_nonref //= 0;
    die "$pack not exists" unless -s $pack;
    my ($ref_count, $ref_dp, $nr_count, $nr_dp) = (0, 0, 0, 0);
    my $cmd = qq# $vg depth --threads 1 --min-coverage 1 -k $pack $xg 2>>$logfile #;
    # cost 8.6G mem for bov graph genome
    open(my $I, "$cmd |");
    while(<$I>) {
        chomp;
        next unless $_;
        my @F = split(/\t/, $_);
        # my $is_ref = $F[0]=~/^\d+$/ ? 1 : 0;
        my $dp=$F[2];
        next if ($dp == 0);
        if ($merge_ref_nonref or $F[0]=~/^\d+$/) { # ref
            $ref_count++;
            $ref_dp += $dp;
        } else { # non-ref
            $nr_count++;
            $nr_dp += $dp;
        }
    }
    my $O = open_out_fh($outfile);
    if ($merge_ref_nonref==1) {
        say $O join "\t", qw/ #sid count dp  /;
        say $O join "\t", $sid, $ref_count, $ref_dp;
    } else {
        say $O join "\t", qw/ #sid ref_count ref_dp nr_count nr_dp /;
        say $O join "\t", $sid, $ref_count, $ref_dp, $nr_count, $nr_dp;
    }
    close $O;
}




