#!/usr/bin/env perl
#===============================================================================
#
#         FILE: sv2snp_indel.pl
#
#        USAGE: ./sv2snp_indel.pl  
#
#  DESCRIPTION:  Merge small alleles for SVs
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
my $sv_min_dp = 10;
my $max_len_tomerge = 5;
my $debug = 0;
my $min_diff_fold = 5;

sub usage {
    print <<EOF;
Usage: $0 [options]
Options:
    -i | --invcf	    input vcf file
    -o | --outvcf	    output file, default: -
    -t | --threads  	thread number, default: $thread
    --sv_min_dp         sv_min_dp, default: $sv_min_dp
    --max_len_tomerge   max_len_tomerge, default: $max_len_tomerge
    --min_diff_fold     min different length (fold) to merge, default: $min_diff_fold
EOF
    exit;
}
$min_diff_fold = 1 if $min_diff_fold < 1;

GetOptions (
    'help|h!' => \$opt_help,
    'i|invcf=s' => \$in,
    'o|outvcf=s' => \$out,
    't|threads=i' => \$thread,
    'max_len_tomerge=i' => \$max_len_tomerge,
    'sv_min_dp=i' => \$sv_min_dp,
    'debug!' => \$debug,
);

&usage() if $opt_help or not $in or not $out;

my $I = open_in_fh($in);
my $O = open_out_fh($out, 8);


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
    die $line unless @alts;
    my @ref_alts = ($ref, @alts);
    #my $replace = &gen_replace_by_type($ref, \@alts, \@ref_alts);
    my $replace = &gen_replace_by_length_only($ref, \@alts, \@ref_alts);
    #say STDERR Dumper $replace;
    my ($replace2, $newalts) = &gen_replace2($replace, scalar(@ref_alts)-1);
    foreach my $isid (9..$#F) {
        my $alle_infos = $F[$isid];
        $alle_infos =~ /^([^:]+)(:|$)/ or die;
        my $GT = $1;
        if ($GT=~/\./) {
            $F[$isid] = './.';
            next;
        }
        my @alleles = split /[\|\/]/, $GT;
        die unless scalar(@alleles)==2;
        foreach my $ialt (0..$#alleles) {
            my $alle = $alleles[$ialt];
            #if (exists $replace{$alle}) {
            #    $alleles[$ialt] = $replace{$alle};
            #}
            $alleles[$ialt] = $$replace2{$alle} // die "$alle";
        }
        $F[$isid] = join "/", @alleles;
    }
    if (scalar(@$newalts)==0) {
        $F[4] = '.';
    } else {
        $F[4] = join ",", @ref_alts[@$newalts]; 
    }
    if($debug==1 and $F[1] eq '90651') { # debug
        say STDERR "replace: " . Dumper $replace;
        say STDERR "replace2: " . Dumper $replace2;
        say STDERR "newalts: " . Dumper $newalts;
        die;
    }
    return join("\t", @F);
}

sub gen_replace_by_length_only {
    my ($ref, $alts, $ref_alts) = @_;
    my $ref_alts_maxi = scalar(@$ref_alts)-1;
    my $reflen = length $ref;
    my @lengths = map {length} @$ref_alts;
    my @lengths_alts = @lengths[1..$ref_alts_maxi];
    my $alt_max_len = max @lengths_alts;
    my $alt_min_len = min @lengths_alts;
    my $ref_alt_max_len = max($reflen, $alt_max_len);
    my $ref_alt_min_len = min($reflen, $alt_min_len);
    my $diff_len = abs($reflen - $alt_max_len);
    my %replace;
    if ($ref_alt_max_len >= $sv_min_dp and $ref_alt_min_len <= $max_len_tomerge and
            $ref_alt_min_len*$min_diff_fold <= $ref_alt_max_len ) {
        # merge small alles
        my $max_len_tomerge_now = min($max_len_tomerge, int($ref_alt_max_len/$min_diff_fold));
        if($reflen < $max_len_tomerge_now) {
            # merge small alles to ref
            foreach my $ialt (1..$ref_alts_maxi) {
                $replace{$ialt} = 0 if $lengths[$ialt] <= $max_len_tomerge_now;
            }
        } else {
            my $shortest_ialt;
            my $shortest_len = 999999999;
            my @short_ialts;
            foreach my $ialt (1..$ref_alts_maxi) {
                my $len = $lengths[$ialt];
                if($len <= $max_len_tomerge_now) {
                    push @short_ialts, $ialt;
                }
                if ($len < $shortest_len) {
                    $shortest_ialt = $ialt;
                    $shortest_len = $len;
                }
            }
            if (scalar(@short_ialts) >= 2) {
                foreach my $ialt (@short_ialts) {
                    $replace{$ialt} = $shortest_ialt if $ialt != $shortest_ialt;
                }
            }
        }
    }
    return \%replace;
}

sub gen_replace_by_type {
    my ($ref, $alts, $ref_alts) = @_;
    my $ref_alts_maxi = scalar(@$ref_alts)-1;
    my $reflen = length $ref;
    my $alt_max_len = max map {length} @$alts;
    my $alt_min_len = min map {length} @$alts;
    my $ref_alt_max_len = max($reflen, $alt_max_len);
    my $ref_alt_min_len = min($reflen, $alt_min_len);
    my $diff_len = abs($reflen - $alt_max_len);
    my %replace;
    if ($reflen >= $sv_min_dp and $alt_max_len>=$sv_min_dp) {
        # cannot determine type
    } elsif ($reflen <= $sv_min_dp and $alt_max_len<=$sv_min_dp) {
        # cannot determine type
    } elsif ($reflen <= $max_len_tomerge and $alt_max_len >= $sv_min_dp and 
            $reflen*$min_diff_fold <= $alt_max_len) { # ins
        my $max_len_tomerge_now = min($max_len_tomerge, int($alt_max_len/$min_diff_fold));
        foreach my $ialt (1..$ref_alts_maxi) {
            my $len = length $$ref_alts[$ialt];
            if ($len <= $max_len_tomerge_now) {
                $replace{$ialt} = 0 if $ialt != 0;
            }
        }
    } elsif ($reflen >= $sv_min_dp and $alt_min_len <= $max_len_tomerge and
            $alt_min_len*$min_diff_fold <= $reflen) { # del
        my $max_len_tomerge_now = min($max_len_tomerge, int($reflen/$min_diff_fold));
        my $shortest_ialt = 1;
        my $shortest_len = length $$alts[0];
        foreach my $ialt (1..$ref_alts_maxi) {
            my $len = length $$ref_alts[$ialt];
            if ($len < $shortest_len) {
                $shortest_ialt = $ialt;
                $shortest_len = $len;
            }
        }
        foreach my $ialt (0..$ref_alts_maxi) {
            my $len = length $$ref_alts[$ialt];
            if ($len <= $max_len_tomerge_now) {
                $replace{$ialt} = $shortest_ialt if $ialt != $shortest_ialt;
            }
        }
    }
    return \%replace;
}

sub gen_replace2 {
    my ($replace, $nref_alts) = @_;
    my %replace2;
    # init
    foreach my $ialt (0..$nref_alts) {
        $replace2{$ialt} = $ialt;
    }
    # replace
    foreach my $ialt (keys %$replace) {
        my $dst = $$replace{$ialt};
        $replace2{$ialt} = $dst;
    }
    # fix
    my %to_subtracts;
    foreach my $ialt (sort {$a<=>$b} keys %$replace) {
        my $ialt1 = $ialt+1;
        foreach my $ii ($ialt1..$nref_alts) {
            $to_subtracts{$ii}++;
        }
    }
    foreach my $ialt (sort keys %to_subtracts) {
        next if exists $$replace{$ialt};
        my $to_subtract = $to_subtracts{$ialt};
        $replace2{$ialt} -= $to_subtract;
    }
    # fix2
    foreach my $ialt (keys %replace2) {
        my $dst_ori = $replace2{$ialt};
        if (exists $$replace{$ialt}) {
            my $dst = $replace2{$dst_ori};
            $replace2{$ialt} = $dst;
        }
    }
    my @newalts;
    foreach my $ialt (1..$nref_alts) {
        next if exists $$replace{$ialt};
        push @newalts, $ialt unless $ialt~~@newalts;
    }
    return(\%replace2, \@newalts);
}


