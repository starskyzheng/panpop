#!/usr/bin/env perl
#===============================================================================
#
#         FILE: flt_maf_per_alle.pl
#
#        USAGE: ./flt_maf_per_alle.pl  
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
#      CREATED: 12/31/2021 05:41:54 PM
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
use zzIO;
use MCE::Flow;
use MCE::Candy;
use Getopt::Long;


my ($in, $out, $threads, $min_maf, $max_miss_freq, $opt_help);

sub help {
    print STDERR <<EOF;
usage: perl xx.pl --in xx.vcf --out xx.flt.vcf 
        [--threads 1]
        [--min_maf 0.05]
        [--max_miss_freq 0.3]
EOF
    say STDERR @_ if @_;
exit(-1);
}

my $ARGVs = join ' ', @ARGV;

GetOptions (
        'help|h!' => \$opt_help,
        'in|i=s' => \$in,
        'threads|t=i' => \$threads,
        'min_maf=s' => \$min_maf,
        'max_miss_freq=s' => \$max_miss_freq,
        'out|o=s' => \$out,
);

&help() if $opt_help;
&help() unless $in and $out;

$min_maf //= 0.05;
$max_miss_freq //= 0.3;
$threads //= 1;

my $I = open_in_fh($in);
my $O = open_out_fh($out, 8);

my @header;
while(<$I>) {
    chomp;
    next unless $_;
    if (/^##/) {
        say $O $_;
        next;
    } elsif (/^#/) {
       @header = split(/\t/, $_);
       say $O '##INFO=<ID=AF,Number=1,Type=String,Description="Allele count">';
       say $O qq(##CommandLine="$0 $ARGVs");
       say $O $_;
       last; 
    } else {
        die;
    }
}

my $max_idi = scalar(@header)-1;
my $id_count = $max_idi-8;
$min_maf = $min_maf * $id_count * 2;
$min_maf = int($min_maf)+1 if $min_maf != int($min_maf);

say STDERR "Total `$id_count` ids. Min alle is `$min_maf`";


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
    max_workers => $threads,  
    #max_workers => 1,
        #use_slurpio => 0,
    gather => MCE::Candy::out_iter_fh($O)
    },
    \&flt_mce,
    $I;

exit;

#while(<$I>) {
sub process_line {
    chomp;
    return undef unless $_;
    my @F = split(/\t/, $_);
    my $update_ia = process_alts_freq(\@F);
    return undef unless $update_ia;
    for(my $i=9; $i<=$max_idi; $i++) {
        $F[$i]=~m#^(\d+)[/|](\d+)# or next;
        my ($ia1new, $ia2new) = ($$update_ia{$1}, $$update_ia{$2});
        if ($ia1new<0 or $ia2new<0) {
            $F[$i]=~s#^(\d+)[/|](\d+)#./.#;
        } else {
            $F[$i]=~s#^(\d+)[/|](\d+)#$ia1new/$ia2new#;
        }
    }
    return join "\t", @F;
}

sub process_alts_freq {
    my ($F) = @_;
    my @ref_alts = ($$F[3], split(/,/, $$F[4]));
    my %freqs;
    for(my $i=9; $i<=$max_idi; $i++) {
        $$F[$i]=~m#^(\d+)[/|](\d+)# or next; # miss
        #my ($a1, $a2) = ($ref_alts[$1], $ref_alts[$2]);
        my ($ia1, $ia2) = ($1, $2); # allele in interger
        $freqs{$ia1}++;
        $freqs{$ia2}++;
    }
    my %ia_del;
    my $ref_alts_len_ori = $#ref_alts;
    for(my $i=$ref_alts_len_ori; $i>=1; $i--) { # reverse
        my $freq = $freqs{$i} // 0;
        if ($freq < $min_maf) {
            #push @ia_del, $i;
            $ia_del{$i}++;
            splice @ref_alts, $i, 1;
        }
    }
    if (scalar(@ref_alts)<2) {
        return undef; #### ref only
    }
    $$F[3] = $ref_alts[0];
    $$F[4] = join ',', @ref_alts[1..$#ref_alts];
    my $del_num=0;
    my %update_ia;
    my $del_freq=0;
    my @info_AF;
    my $all_not_miss_count=0;
    for(my $i=0; $i<=$ref_alts_len_ori; $i++) {
        my $freq = $freqs{$i} // 0;
        if (exists $ia_del{$i}) {
            $del_num++;
            $update_ia{$i} = -1; ##
            $del_freq += $freq;
            next;
        }
        $all_not_miss_count += $freq;
        my $newi = $i - $del_num;
        push @info_AF, "$newi:$freq";
        $update_ia{$i} = $newi;
    }
    return undef if $all_not_miss_count < $id_count - $max_miss_freq;
    push @info_AF, "DEL:$del_freq:".scalar(keys %ia_del) if $del_freq;
    my $F7_append = "AF=" . join(',', @info_AF);
    if ($$F[7] eq '.') {
        $$F[7] = $F7_append;
    } else {
        $$F[7] = $$F[7] . ";$F7_append";
    }
    #say join ' ', %ia_del;
    return(\%update_ia);
}
