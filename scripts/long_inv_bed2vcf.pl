#!/usr/bin/env perl
use strict;
use warnings;
use v5.24;
use Getopt::Long;
no warnings 'experimental::smartmatch';
use Data::Dumper;

use FindBin qw($Bin);
use lib "$Bin/../lib";
use zzBed;
use zzIO;


my ($ref_fasta_file, $out, $opt_help, $not_rename);
my ($sample_name );

my $max_len = 50000;
my $min_dp = 5;

sub help {
    my ($info) = @_;
    say STDERR $info if $info;
    say <<EOF;
Usage: $0 -r <ref.fa> -o <out.vcf> <in1.bed> <in2.bed> ...
    -r, --ref <ref.fa>      reference fasta file
    -o, --out <out.vcf>     output vcf file
    -h, --help              print this help message
    -s, --sample            Sample name
    --max_len <int>         max length of SV, default $max_len
    --min_dp <int>          min reads depth of SV, default $min_dp
EOF
    exit(-1);
}

my $ARGVs = join ' ', $0, @ARGV;

GetOptions (
        'help|h!' => \$opt_help,
        'r|ref=s' => \$ref_fasta_file,
        'o|out=s' => \$out,
        's|sample=s' => \$sample_name,
        'max_len=s' => \$max_len,
        'd|min_dp=i' => \$min_dp,
);


my @in_beds = @ARGV;

&help() if $opt_help;
&help() unless @in_beds;
&help() unless defined $out;
&help() unless defined $ref_fasta_file;



my $bedobj = new zzBed({max_len=>$max_len});
my $ref_fasta = &read_ref_fasta($ref_fasta_file);

my $O = open_out_fh($out);
my @header;

my $vcfheader_chrlen = &gen_vcfheader_chrlen($ref_fasta);

print $O <<EOF;
##fileformat=VCFv4.2
##source=Assemblytics
$vcfheader_chrlen
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SUPP,Number=1,Type=String,Description="Number of samples supporting the variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##CommandLine="$ARGVs"
EOF
say $O join "\t", qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT /, $sample_name;

foreach my $bed (@in_beds) {
    $bedobj->read_add_bed($bed);    
}

$bedobj->sort_bed();
#die Dumper $bedobj->{array};
$bedobj->cal_cov_count();
$bedobj->cov_count_filter();
my $final_bed = $bedobj->cov_count_filter_merge();
foreach my $chr (sort keys %$final_bed) {
    foreach my $pos (sort {$a <=> $b} keys %{$$final_bed{$chr}}) {
        my ($end, $supp, $gt) = $$final_bed{$chr}{$pos}->@*;
        my $len = $end - $pos + 1;
        my $ref = substr($$ref_fasta{$chr}, $pos - 1, $len);
        my $alt = &revcomp($ref);
        my @newline = ($chr, $pos, '.', $ref, $alt, 
                        '.', '.', "SVTYPE=INV;SVLEN=$len;SUPP=$supp", 
                        'GT', $gt);
        say $O join "\t", @newline;
    }
}

exit;




sub gen_vcfheader_chrlen {
    my ($fa) = @_;
    my @outs;
    foreach my $chr(sort keys %$fa) {
        my $seq = $$fa{$chr};
        my $len = length($seq);
        push @outs, "##contig=<ID=$chr,length=$len>";
    }
    return join "\n", @outs;
}

sub max {
    my ($a, $b) = @_;
    return $a > $b ? $a : $b;
}


sub revcomp {
    my ($seq) = @_;
    my $rev = reverse $seq;
    $rev =~ tr/ACGTacgt/TGCAtgca/;
    $rev =~ s/[^ACGTacgt]/N/g;
    return $rev;
}



exit;

sub read_ref_fasta {
    my $ref_fasta_file = shift;
    my $I = open_in_fh($ref_fasta_file);
    my %ref;
    my $id;
    while(<$I>) {
        chomp;
        if(/^>(\S+)/) {
            $id = $1;
            next;
        }
        $ref{$id} .= $_;
    }
    return \%ref;
}

