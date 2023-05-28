#!/usr/bin/env perl
use strict;
use warnings;
use v5.10;
use Getopt::Long;
use MCE::Flow;
use MCE::Candy;
no warnings 'experimental::smartmatch';

# Possible attempt to put comments in qw() list
no warnings 'qw';

use FindBin qw($Bin);
use lib "$Bin/../lib";
use zzIO;


my ($in, $ref_fasta_file, $out, $opt_help, $software, $not_rename);
my ($query_fasta_file, $out_inv);

my $max_len = 50000;
my $min_dp = 5;
my @known_softwares = qw/pbsv svim cutesv sniffles2/;
my $threads = 1;

sub help {
    my ($info) = @_;
    say STDERR $info if $info;
    say <<EOF;
Usage: $0 -i <in.bed> -r <ref.fa> -q <query.fa> -o <out.vcf>
    -i, --in <in.bed>       input bed file
    -r, --ref <ref.fa>      reference fasta file
    -q, --query <query.fa>  query fasta file
    -o, --out <out.vcf>     output vcf file
    -t, --threads <int>     threads, default $threads
    -h, --help              print this help message
    --max_len <int>         max length of SV, default $max_len
    --min_dp <int>          min reads depth of SV, default $min_dp
EOF
    exit(-1);
}

my $ARGVs = join ' ', $0, @ARGV;

GetOptions (
        'help|h!' => \$opt_help,
        'i|in=s' => \$in,
        'r|ref=s' => \$ref_fasta_file,
        'q|query=s' => \$query_fasta_file,
        'o|out=s' => \$out,
        #'O|out_inv=s' => \$out_inv,
        'max_len=s' => \$max_len,
        'd|min_dp=i' => \$min_dp,
        't|threads=i' => \$threads,
        's|software=s' => \$software,
);

&help() if $opt_help;
&help("No -i input") unless defined $in;
&help("No -o input") unless defined $out;
#&help() unless defined $software;
&help("No -r input") unless defined $ref_fasta_file;

$software = lc($software) if defined $software;
&help("software error!") if(defined $software and ! $software ~~ @known_softwares);


say STDERR "Now reading reference fasta file: $ref_fasta_file";
my $ref_fasta = &read_ref_fasta($ref_fasta_file);
say STDERR "Now reading query fasta file: $query_fasta_file";
my $query_fasta = &read_ref_fasta($query_fasta_file);

my $I = open_in_fh($in);
my $O = open_out_fh($out);
my @header;
say STDERR "Now process bed file: $in";
my $vcfheader_chrlen = &gen_vcfheader_chrlen($ref_fasta);

print $O <<EOF;
##fileformat=VCFv4.2
##source=Assemblytics
$vcfheader_chrlen
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
say $O join "\t", qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Assemblytics/;
my $bedheader = <$I>;

if($threads==1) {
    while(my $line = <$I>) {
        my $vcf_line = &bed2vcf_line($line);
        next unless defined $vcf_line;
        say $O $vcf_line;
    }
} else {
    MCE::Flow->init(
        max_workers => $threads,
        chunk_size => 1,
        #use_slurpio => 1,
        gather => MCE::Candy::out_iter_fh($O),
    );
    mce_flow_f sub {
        my ($mce, $chunk_ref, $chunk_id) = @_;
        my $vcf_line = &bed2vcf_line($$chunk_ref[0]);
        $mce->gather($chunk_id) unless defined $vcf_line;
        $mce->gather($chunk_id, "$vcf_line\n");
    }, $I;
    MCE::Flow->finish();
}

say STDERR "Done!";

exit;


sub bed2vcf_line {
    my ($line) = @_;
    chomp $line;
    my @F = split /\t/, $line;
    my $dp;
    # svid CHROM1, POS1, CHROM2, POS2, SVTYPE, SVLEN, SAMPLE_ID, SV_ID
    #  0      1      2      3      4       5      6         7      8
    my ($svid, $chr1, $pos1, $chr2, $pos2, $svtype, $sample_id, $sv_id_ori) = @F;
    my $sv_type = $F[6];
    next if $chr1 ne $chr2;
    return undef if $svlen > $max_len;
    #                #CHROM  POS  ID    REF     ALT   QUAL  FILTER  INFO  FORMAT sample
    my @newline = ($chr1, $pos1, $svid, 'REF', 'ALT' , '.', 'PASS', '.', 'GT', '1/1');
    if( $svtype eq 'DEL' or $svtype eq 'INS') {
        my $svlen_ref = 1 + abs($ref_stop - $ref_start);
        $F[9]=~/(^[^:]+):(\d+)-(\d+):([\-+])$/ or die;
        my $qry_chr = $1;
        my $qry_start = $2;
        my $qry_end = $3;
        my $qry_strand = $4;
        my $qry_len = 1 + abs($qry_end - $qry_start);
        say STDERR "$qry_start, $qry_len, $svlen_ref";
        exists $$ref_fasta{$ref_chr} // die "! Ref not exists $ref_chr ! $F[9] ! $line";
        die if $$ref_fasta{$ref_chr} eq '';
        exists $$query_fasta{$qry_chr} // die "! Ref not exists $qry_chr ! $F[9] ! $line";
        my $ref_chrnow = \$$ref_fasta{$ref_chr};
        my $ref_seqnow = substr($$ref_chrnow, $ref_start-1, $svlen_ref);
        my $alt_chrnow = \$$query_fasta{$qry_chr};
        my $alt_seqnow = substr($$alt_chrnow, $qry_start-1, $qry_len);
        $newline[3] = $ref_seqnow;
        $newline[4] = $alt_seqnow;
        if($F[6] eq 'Deletion') {
            $newline[7] = "SVTYPE=DEL;SVLEN=$svlen";
        } elsif($F[6] eq 'Insertion') {
            $newline[7] = "SVTYPE=INS;SVLEN=$svlen";
        }
        if($qry_strand ne $F[5]) {
            $newline[4] = &revcomp($newline[4]);
        }
    } else {
        return undef;
    }
    #$F[7] = '.'; # remove all infos
    return join "\t", @newline;
}

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

