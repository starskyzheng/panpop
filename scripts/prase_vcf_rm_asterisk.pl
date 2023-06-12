#!/usr/bin/env perl
#===============================================================================
#
#         FILE: .pl
#
#        USAGE: .pl  
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
#      CREATED: 6/04/2023 11:54:32 AM
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
use Getopt::Long;

use prase_vcf qw/read_ref read_vcf_header/;


my ($infile, $ref_fa, $outfile);

sub USAGE {
    my ($info) = @_;
    my $usage = <<USAGE;
    Options:
    -i|--infile <str>   infile
    -r|--ref_fa <str>   ref_fa
    -o|--outfile <str>  outfile
    -h|--help         help
USAGE
    if($info) {
        say STDERR "===============";
        say STDERR $info;
        say STDERR "===============";
    }
    print $usage;
    exit;
}

GetOptions(
    "i|infile=s" => \$infile,
    "r|ref_fa=s" => \$ref_fa,
    "o|outfile=s" => \$outfile,
    "h|help!" => \&USAGE,
);

&USAGE("No --infile") unless $infile;
&USAGE("No --ref_fa") unless $ref_fa;
&USAGE("No --outfile") unless $outfile;

say STDERR "Now reading ref fasta : ";
my ($REF_SEQS, $REF_SEQS_STAT) = read_ref($ref_fa);
say STDERR "Done reading ref fasta $REF_SEQS_STAT";


my $I = open_in_fh($infile);
my $O = open_out_fh($outfile);

my @old_stat; # chr   end  [old_line_split]   [seqs_for_ids] 
while(<$I>) {
    chomp;
    if(/^#/) {
        say $O $_;
        next;
    }

    my @F = split /\t/;
    my $chr = $F[0];
    my $pos = $F[1];
    my $ref = $F[3];
    my $alts = $F[4];
    my @alts = split /,/, $alts;
    my $end = $pos + length($ref) - 1;
    if(@old_stat and $chr eq $old_stat[0] and $pos <= $old_stat[1]) {
        # two muts were overlap, merge them
        die "$chr $pos $old_stat[1]" unless $pos == $old_stat[1];
        if($alts=~m/\*/) {
            my $next_base = substr($REF_SEQS->{$chr}, $end, 1);
            #say $next_base;
            &alts_append_base(\@alts, $next_base);
            $F[3] .= $next_base;
            $F[4] = join ',', @alts;
            $old_stat[1] = 1+$end;
        }
        my $newline = &merge_two_line($old_stat[2], \@F);
        $old_stat[2] = $newline;
    } else {
        if(@old_stat) {
            # not overlap, print old_stat to new vcf
            my $old_line = join "\t", $old_stat[2]->@*;
            say $O $old_line;
            @old_stat = ();
        }
        if($alts=~m/\*/) {
            # not overlap, but have asterisk, save it
            my $next_base = substr($REF_SEQS->{$chr}, $end, 1);
            &alts_append_base(\@alts, $next_base);
            $ref .= $next_base;
            $F[3] = $ref;
            $F[4] = join ',', @alts;
            @old_stat = ($chr, 1+$end, \@F, []);
        } else {
            # not overlap, no asterisk, just print it
            say $O $_;
        }
    }
}

if(@old_stat) {
    my $old_line = join "\t", $old_stat[2]->@*;
    say $O $old_line;
    @old_stat = ();
}

exit;


sub merge_two_line {
    my ($line1, $line2) = @_;
    my ($new_alts_id_seqs, $newref) = &merge_alleles_each_id($line1, $line2);
    my $newline = &rebuild_line($line1, $newref, $new_alts_id_seqs);
    return($newline);
}


sub rebuild_line {
    my ($line, $newref, $new_alts_id_seqs) = @_;
    my $idi_max = scalar(@$line)-1;
    my @new_line = ();
    $new_line[0] = $line->[0];
    $new_line[1] = $line->[1];
    $new_line[2] = $line->[2];
    $new_line[3] = $newref;
    my ($new_alts, $new_gts) = &build_new_alts($new_alts_id_seqs, $newref);
    $new_line[4] = join ",", @$new_alts;
    $new_line[5] = '.'; # QUAL;
    $new_line[6] = '.'; # FILTER;
    $new_line[7] = '.'; # INFO;
    $new_line[8] = 'GT'; # FORMAT;
    foreach my $idi (9..$idi_max) {
        if(! defined( $$new_gts{$idi} )) {
            $new_line[$idi] = './.';
            next;
        }
        my $newgt_id = join("/", $$new_gts{$idi}->@*);
        $new_line[$idi] = $newgt_id;
    }
    return (\@new_line);
}


sub build_new_alts {
    my ($new_alts_id_seqs, $newref) = @_;
    my %new_refalts = ($newref=>0);
    my %new_gts;
    my $new_refalts_count_now = 1;
    foreach my $idi (sort {$a<=>$b} keys %$new_alts_id_seqs) {
        my @gts;
        if(! defined $$new_alts_id_seqs{$idi}) {
            $new_gts{$idi} = undef;
            next;
        }
        my ($a1, $a2) = @{$new_alts_id_seqs->{$idi}};
        foreach my $seq($a1, $a2) {
            if(exists $new_refalts{$seq}) {
                my $gt = $new_refalts{$seq};
                push @gts, $gt;
            } else {
                my $gt = $new_refalts_count_now;
                $new_refalts{$seq} = $gt;
                push @gts, $gt;
                $new_refalts_count_now++;
            }
        }
        $new_gts{$idi} = \@gts;
    }
    my @new_refalts = sort {$new_refalts{$a}<=>$new_refalts{$b}} keys %new_refalts;
    my @new_alts = @new_refalts[1..$#new_refalts];
    return (\@new_alts, \%new_gts);
}

sub merge_alleles_each_id {
    my ($line1, $line2) = @_; # old new
    my $idi_max = scalar(@$line1)-1;
    my %new_alts;
    foreach my $idi (9..$idi_max) {
        my ($seq1_a1, $seq1_a2) = &line2id_seqs($line1, $idi);
        my ($seq2_a1, $seq2_a2) = &line2id_seqs($line2, $idi);
        if(defined($seq1_a1) and defined($seq2_a1)) {
            my $rm1 = substr($seq1_a1, -1, 1, "");
            my $rm2 = substr($seq1_a2, -1, 1, "");
            unless($rm1 eq $rm2) {
                say STDERR "@$line1";
                die "Error: $idi, $rm1, $rm2";
            }
            my $seqm_a1 = $seq1_a1 . $seq2_a1;
            my $seqm_a2 = $seq1_a2 . $seq2_a2;
            $new_alts{$idi} = [$seqm_a1, $seqm_a2];
        } else {
            $new_alts{$idi} = undef;
        }
    }
    my $ref1 = $$line1[3];
    my $ref2 = $$line2[3];
    substr($ref1, -1, 1, "");
    my $newref = $ref1 . $ref2;
    return(\%new_alts, $newref);
}

sub line2id_seqs {
    my ($line, $idi, $remove_asterisk) = @_;
    my @alts = split(/,/, $$line[4]);
    if(defined $remove_asterisk and $remove_asterisk==1) {
        foreach my $ialt (0..$#alts) {
            $alts[$ialt] =~ s/\*//;
        }
    }
    my $ref = $$line[3];
    my @ref_alts = ($ref, @alts);
    my $gts = $$line[$idi];
    $gts=~m#^([\d.]+)[/|]([\d.]+)# or die;
    my @gts = ($1, $2);
    if($1 eq '.' or $2 eq '.') {
        return(undef);
    }
    my @seqs = map {$ref_alts[$_]} @gts;
    return(@seqs);
}



sub alts_append_base {
    my ($alts, $base) = @_;
    my $ialts_max = scalar(@$alts)-1;
    foreach my $ialts (0..$ialts_max) {
        if($$alts[$ialts] eq '*') {
            $$alts[$ialts] = $base;
        } else {
            $$alts[$ialts] .= $base;
        }
        die if $$alts[$ialts] =~ m/\*/;
    }
}

