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
use Tie::CharArray;


my ($in, $out, $opt_help);# = @ARGV;
my $thread = 1;
my $sv_min_dp = 10;
my $max_len_tomerge = 5;
my $min_lendiff_tomerge = 100;
my $max_diff_fold = 0.6;
my $debug = 0;
my $min_diff_fold = 5;
my $enable_norm_alle = 1;
$out = '-';
my $mode = 'by_length'; # by_length by_type
my $by_type_force_pav = 0;
my $by_length_merge_big = 0; # by_length mode: merge big allele to longest allele. 0: disable, 1: enable
my $by_length_merge_small = 1; # by_length mode: merge small allele to shortest allele. 0: disable, 1: enable

sub usage {
    print <<EOF;
Usage: $0 [options]
Options:
    -i | --invcf       <FILE>  input vcf file
    -o | --outvcf      <FILE>  output file, default: $out
    -t | --threads     <INT>   thread number, default: $thread
    --sv_min_dp        <INT>   sv_min_bp, default: $sv_min_dp
    --max_len_tomerge  <INT>   max_len_tomerge SMALL allele, 0 for +∞, -1 to disable, default: $max_len_tomerge
    --min_diff_fold    <INT>   min different length (fold) to merge SMALL allele, default: $min_diff_fold
    --min_lendiff_tomerge  <INT>   min_len_diff_tomerge for BIG allele, 0 for -∞, default: $min_lendiff_tomerge
    --max_diff_fold    <FLOAT> max different length (fold) to merge BIG allele, 1 to disable, default: $max_diff_fold
    --enable_norm_alle <INT>   enable normalize alleles, default: $enable_norm_alle
    --mode             <STR>   by_length or by_type, default: $mode
    --by_type_force_pav
    --by_length_merge_big
    --by_length_merge_small
EOF
    exit;
}

$min_diff_fold = 1 if $min_diff_fold < 1;

my $ARGVs = join " ", @ARGV;

GetOptions (
    'help|h!' => \$opt_help,
    'i|invcf=s' => \$in,
    'o|outvcf=s' => \$out,
    't|threads=i' => \$thread,
    'max_len_tomerge=i' => \$max_len_tomerge,
    'min_diff_fold=i' => \$min_diff_fold,
    'sv_min_dp=i' => \$sv_min_dp,
    'enable_norm_alle=i' => \$enable_norm_alle,
    'min_lendiff_tomerge=i' => \$min_lendiff_tomerge,
    'max_diff_fold=f' => \$max_diff_fold,
    'debug!' => \$debug,
    'mode=s' => \$mode,
    'by_type_force_pav=i' => \$by_type_force_pav,
    'by_length_merge_big=i' => \$by_length_merge_big,
    'by_length_merge_small=i' => \$by_length_merge_small,
);

$max_len_tomerge = 999999999 if $max_len_tomerge == 0;
$max_len_tomerge = 0 if $max_len_tomerge <= -1;
&usage() if $opt_help or not $in or not $out;

my $I = open_in_fh($in);
my $O = open_out_fh($out, 8);

while(<$I>) {
    if(/^##/) {
        print $O $_;
        next;
    } elsif (/^#/) {
        say $O qq(##CommandLine="$0 $ARGVs");
        print $O $_;
        last;
    }
}

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

#close $I;
exit 0;

sub process_line {
    my ($line) = @_;
    chomp $line;
    return undef unless $line;
    if ($line=~/^#/) {
        return $line;
    }
    my @F = split /\t/, $line;
    my $ref = uc($F[3]);
    my @alts = split /,/, uc($F[4]);
    die $line unless @alts;
    my @ref_alts = ($ref, @alts);
    my $replace;
    if($mode eq 'by_length') {
        $replace = &gen_replace_by_length_only($ref, \@alts, \@ref_alts);
    } elsif ($mode eq 'by_type') {
        $replace = &gen_replace_by_type($ref, \@alts, \@ref_alts);
    } else {
        die "unknown mode: $mode";
    }
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
        my @newalts = @ref_alts[@$newalts];
        $F[4] = join ",", @newalts;
        &norm_alle_alles([$ref, @newalts], \@F) if $enable_norm_alle==1;
    }
    if($debug==1 and $F[1] eq '90651') { # debug
        say STDERR "replace: " . Dumper $replace;
        say STDERR "replace2: " . Dumper $replace2;
        say STDERR "newalts: " . Dumper $newalts;
        die;
    }
    return join("\t", @F);
}

sub norm_alle_alles {
    my ($alles, $F) = @_;
    my ($skip_start, $skip_end) = &norm_alle_alles_cal_length($alles);
    my $reflen = length $alles->[0];
    foreach my $ialt (0..$#$alles) {
        my $seq = $alles->[$ialt];
        my $len = length $seq;
        $seq = substr($seq, $skip_start, $len-$skip_start-$skip_end);
        $alles->[$ialt] = $seq;
    }
    $$F[1] += $skip_start;
    $$F[3] = $alles->[0];
    $$F[4] = join ",", @$alles[1..$#$alles];
}

sub norm_alle_alles_cal_length {
    my ($alles) = @_;
    my @tie_seqs;
    my ($skip_start, $skip_end) = (0, 0);
    my $reflen = length $alles->[0];
    my @alt_lens = map {length $_} @$alles;
    foreach my $ialt (0..$#$alles) {
        my $seq = $alles->[$ialt];
        tie $tie_seqs[$ialt]->@*, 'Tie::CharArray', $seq;
    }
    # skip from end
    LOOP_SKIP_END: while(1) {
        my $last = $tie_seqs[0]->[-$skip_end-1];
        foreach my $ialt (1..$#$alles) {
            my $len = $alt_lens[$ialt];
            if ($len <= $skip_end+1 or $tie_seqs[$ialt][-$skip_end-1] ne $last) {
                last LOOP_SKIP_END;
            }
        }
        if($skip_end+1 >= $reflen) {
            last LOOP_SKIP_END;
        } else {
            $skip_end++;
        }
    }
    # skip from start
    LOOP_SKIP_START: while(1) {
        my $first = $tie_seqs[0][$skip_start] // die $skip_start;
        foreach my $ialt (1..$#$alles) {
            my $len = $alt_lens[$ialt];
            if ($len-$skip_end <= 1+$skip_start or $tie_seqs[$ialt][$skip_start] ne $first) {
                #die;
                last LOOP_SKIP_START;
            }
        }
        if($skip_start+1 >= $reflen-$skip_end) {
            last LOOP_SKIP_START;
        } else {
            $skip_start++;
        }
    }

    return ($skip_start, $skip_end);
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
        # 将很短的ref/alt合并
        my $max_len_tomerge_now = min($max_len_tomerge, int($ref_alt_max_len/$min_diff_fold),
                                        max($ref_alt_max_len - $min_lendiff_tomerge,
                                            int($ref_alt_max_len*$max_diff_fold) ) );
        #say STDERR "$max_len_tomerge_now = min($max_len_tomerge, int($ref_alt_max_len/$min_diff_fold));";
        my $min_len_tomerge_now = max($ref_alt_max_len - $min_lendiff_tomerge, 
                                        int($ref_alt_max_len*$max_diff_fold),
                                        min($max_len_tomerge, int($ref_alt_max_len/$min_diff_fold)) );
        my $shortest_ialt;
        my $shortest_len = 999999999;
        my $longest_ialt;
        my $longest_len = 0;
        my @short_ialts;
        my @big_ialts;
        foreach my $ialt (0..$ref_alts_maxi) {
            my $len = $lengths[$ialt];
            if($len <= $max_len_tomerge_now) {
                push @short_ialts, $ialt;
                if ($len < $shortest_len) {
                    $shortest_ialt = $ialt;
                    $shortest_len = $len;
                }
            } elsif ($len >= $min_len_tomerge_now) {
                push @big_ialts, $ialt;
                if ($len > $longest_len) {
                    $longest_ialt = $ialt;
                    $longest_len = $len;
                }
            }
        }
        if ($by_length_merge_small==1 and scalar(@short_ialts) >= 2) {
            if (0 ~~ @short_ialts) {
                $shortest_ialt = 0;
            }
            foreach my $ialt (@short_ialts) {
                $replace{$ialt} = $shortest_ialt if $ialt != $shortest_ialt;
            }
        }
        if ($by_length_merge_big==1 and scalar(@big_ialts) >= 2) {
            if (0 ~~ @big_ialts) {
                $longest_ialt = 0;
            }
            foreach my $ialt (@big_ialts) {
                $replace{$ialt} = $longest_ialt if $ialt != $longest_ialt;
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
        # die "$reflen >= $sv_min_dp and $alt_max_len>=$sv_min_dp";
    } elsif ($reflen <= $sv_min_dp and $alt_max_len<=$sv_min_dp) {
        # cannot determine type
        # die "$reflen <= $sv_min_dp and $alt_max_len<=$sv_min_dp";
    } elsif ($reflen <= $max_len_tomerge and $alt_max_len >= $sv_min_dp and 
            $reflen*$min_diff_fold <= $alt_max_len) { # ins
        my $max_len_tomerge_now = min($max_len_tomerge, int($alt_max_len/$min_diff_fold));
        my ($longest_ialt, @ins_ialts);
        my $longest_ialt_len = 0;
        foreach my $ialt (1..$ref_alts_maxi) {
            my $len = length $$ref_alts[$ialt];
            if ($len <= $max_len_tomerge_now) {
                $replace{$ialt} = 0 if $ialt != 0;
            } elsif($by_type_force_pav==1) {
                # force convert to PAV, force convert to only one INS
                $longest_ialt = $ialt if $len > $longest_ialt_len;
                push @ins_ialts, $ialt;
            }
        }
        if($by_type_force_pav==1) {
            if (0 ~~ @ins_ialts) {
                die;
            }
            foreach my $ialt (@ins_ialts) {
                $replace{$ialt} = $longest_ialt if $ialt != $longest_ialt;
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
            if ($len <= $max_len_tomerge_now) { # del allele
                $replace{$ialt} = $shortest_ialt if $ialt != $shortest_ialt;
            } elsif ($by_type_force_pav==1) {
                # force convert to PAV, non-del to REF
                $replace{$ialt} = 0 if $ialt != 0;
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


