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
#      CREATED: 10/21/2021 02:16:06 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
no warnings qw( experimental::smartmatch );
use v5.24;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Carp qw/confess/;
use zzIO;
use Data::Dumper;
#use MCE::Loop;
use Getopt::Long;
use prase_vcf qw/read_ref/;

sub help() {
    say STDERR "@_" if @_;
        die <<EOF;
usage: perl x.pl [options]
-i|--in_vcf <str>  Input vcf file
-o|--out_vcf <str> Output vcf file
-d|--dp_file <str> Depth file
-r|--ref <str>     Reference.fa
--min_dp <int>     Minimum depth [default: 3]
--min_cov <int>    Minimum coverage [default: 0.8]
-s|--sort_vcf      Output sort vcf file? [default: False]
--chr <str>        Only prase this chromsome
EOF
}



my $APPEND_HEADER = <<'EOF';
##FORMAT=<ID=MDP,Number=1,Type=Float,Description="Mean Read Depth"
##FORMAT=<ID=RACOV,Number=1,Type=Float,Description="Coverage of this position in reference allele"
##FORMAT=<ID=DPSOURCE,Number=1,Type=String,Description="source of Depth. can be both or original or append"
EOF

my ($in_vcf, $out_vcf, $dp_file, $ref_file);
my $mincov = 0.8;
my $mindp = 3;
my $threads = 1;
my $sort_vcf = 0;
my $need_chr;

my ($opt_help);

GetOptions (
        'help|h!' => \$opt_help,
        'i|in_vcf=s' => \$in_vcf,
        'o|out_vcf=s' => \$out_vcf,
        #'t=i' => \$threads,
        'd|dp_file=s' => \$dp_file,
        'min_dp=i' => \$mindp,
        'min_cov=s' => \$mincov,
        'r|ref=s' => \$ref_file,
        's|sort_vcf!' => \$sort_vcf,
        'chr=s' => \$need_chr,
);

$need_chr = undef if $need_chr eq 'all';

unless ($in_vcf && $out_vcf && $dp_file && $ref_file) {
    help();
}

say STDERR "Now: read ref";
my ($refs) = &read_ref($ref_file);
my $sortcmd = &run($dp_file, $in_vcf, $out_vcf);
say STDERR "All done";

if ($sort_vcf==1) {
    say STDERR "Now exit script and then start to sort VCF";
    my $cmd = "bcftools sort -m 40G -o $out_vcf.sort.vcf.gz -O z $out_vcf > $out_vcf.sort.vcf.gz.log 2>&1 ";
    say STDERR $cmd;
    #exec($cmd);
}

exit;

sub run {
    my ($dp_file, $vcf_in, $vcf_out) = @_;
    say STDERR "Now: read depth";
    my ($dps) = &read_dp($dp_file);
    my $I = open_in_fh($vcf_in);
    my $O = open_out_fh($vcf_out, 6);
    say STDERR "Now: process VCF";
    &process_header($I, $O);
    while(my $line = <$I>) { # normal vcf line
        chomp $line;
        next unless $line;
        die if $line=~/^#/;
        my (@F) = split(/\t/, $line);
        my $chr = $F[0];
        if (defined $need_chr) {
            next unless $chr eq $need_chr;
        }
        my $pos = $F[1];
        my $ref = $F[3];
        my $alles = $F[9];
        if (exists $$dps{$chr}{$pos}) { # if this position has depth
            my $reflen = length($ref);
            my ($dp, $cov) = (-9, -9);
            my $str = delete $$dps{$chr}{$pos}{$reflen};
            my $DPSOURCE = 'original';
            if($str) {
                ($dp, $cov) = split(/:/, $str);
                $DPSOURCE = 'both';
            } else {
                #say STDERR "Warn: $chr : $pos : $reflen not in depth file";
            }
            $F[8] .= ":MDP:RACOV:DPSOURCE";
            $F[9] .= ":$dp:$cov:$DPSOURCE";
            if ($alles=~m#^\.[/\|]\.#) {
                substr($F[9], 0, 3) = '0/0';
            }
            say $O join("\t", @F);
            next;
            if (scalar( keys $$dps{$chr}{$pos}->%* ) == 0) {
                delete $$dps{$chr}{$pos};
            }
        } else { # if this position has NO depth
            $F[8] .= ":MDP:RACOV:DPSOURCE";
            $F[9] .= ":-9:-9:original";
            say $O join("\t", @F);
        }
    }
    say STDERR "Now: output left depth of dp_file";
    foreach my $chr (sort keys %$dps) {
        foreach my $pos (sort {$a<=>$b} keys %{$dps->{$chr}}) {
            foreach my $reflen (sort {$a<=>$b} keys %{$dps->{$chr}{$pos}}) {
                my $str = $dps->{$chr}{$pos}{$reflen} // next;#die "?? {$chr}{$pos}{$reflen}";
                my ($dp, $cov) = split(/:/, $str);
                my $ref_now = $$refs{$chr};
                my $refseq_now = substr($ref_now, $pos-1, $reflen);
                my $dpint = sprintf("%.0f", $dp);
                say $O join("\t", $chr, $pos, '.', $refseq_now, '.', '.',
                '.', '.', 'GT:DP:MDP:RACOV:DPSOURCE', "0/0:$dpint:$dp:$cov:append");
            }
        }
    }
}

sub process_header {
    my ($I, $O) = @_;
    my %chr_in_vcf;
    while(my $line = <$I>) {
        chomp $line;
        next unless $line;
        if($line=~/^##/) {
            say $O $line;
            #if()
            next;
        } elsif ($line=~/^#/) {
            print $O $APPEND_HEADER;
            say STDERR [split(/\t/,$line)]->[-1];
            say $O $line;
            last;
        } else {
            die;
        }
    }
}

sub read_dp {
    my ($in) = @_;
    my $I = open_in_fh($in);
    my %ret;
    while(<$I>) {
        chomp;
        next unless $_;
        my ($chr, $pos, $reflen, $sum_dp, $cov_count) = split(/\t/, $_);
        my $dp = $sum_dp / $reflen;
        my $cov = $cov_count / $reflen;
        $cov = 1 if $cov > 1;
        $cov = sprintf("%.3f", $cov);
        $dp = sprintf("%.1f", $dp);
        if ($dp < $mindp or $cov < $mincov) {
            next;
        }
        $ret{$chr}{$pos}{$reflen} = "$dp:$cov";
    }
    return(\%ret);
}


