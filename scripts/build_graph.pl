#!/usr/bin/env perl
use v5.10;
use warnings;
use strict;
use File::Basename;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use zzIO;
use Getopt::Long;
use read_config qw/read_config_yaml/;

my $config = read_config_yaml("$Bin/../config.yaml");

my $vg = $$config{vg} or die "vg not defined in config file!";
my $minigraph = $$config{minigraph} or die "minigraph not defined in config file!";


my ($backbone, $outdir, $opt_help, $rgfa_in, $gfa_in);
$outdir = '.';
my $threads = 24;

sub usage {
    my $msg = shift;
    print STDERR "$msg\n" if $msg;
    print STDERR <<"EOF";
USAGE: perl $0 [options] [genomes]
options:
    --rgfa <file>                  rgfa file by minigraph as input (do not generate by fasta)
    --gfa <file>                   gfa file by vg as input (do not generate by fasta)
    -b | --backbone <file>         Backbone genome file (in fasta)
    -o | --out_dir <file>          Output dir
    -t | --threads <int>           Number of threads (default: $threads)
    -h | --help                    Print this help

The name of chromosome will be specified.
Backbone genome should be numeric only.
Non-Backbone genome should NOT be numeric only.
For example: perl $0 -b $Bin/../example/build_graph.example/MT-human.fa -o $Bin/../example/build_graph.example/out -t 2 $Bin/../example/build_graph.example/MT-chimp.fa $Bin/../example/build_graph.example/MT-orangA.fa
EOF
    exit(1);
}


GetOptions (
        'h|help!' => \$opt_help,
        'b|backbone=s' => \$backbone,
        'o|out_dir=s' => \$outdir,
        't|threads=i' => \$threads,
        'gfa=s' => \$gfa_in,
        'rgfa=s' => \$rgfa_in,
);

my @moregenomes = @ARGV;



die "can't use --gfa and --rgfa at the same time!" if defined $rgfa_in and defined $gfa_in;


if (! -e $outdir) {
    mkdir $outdir or die "$outdir can not be created";
}

my $gfa1 = "$outdir/1.original.rgfa";
my $gfa2 = "$outdir/2.vg.gfa";
my $gfa3 = "$outdir/3.final.gfa";

$gfa1 = $rgfa_in if defined $rgfa_in;
$gfa2 = $gfa_in if(defined $gfa_in);

unless(defined $rgfa_in or defined $gfa_in) {
    &usage() unless $outdir and $backbone and @moregenomes;
    &check_backbone($backbone);
    &check_dup_chrs($backbone, @moregenomes);
    my $cmd_build = "$minigraph -o $gfa1 --inv no -xggs -L 10 -K 4G -t $threads $backbone @moregenomes";
    say STDERR "Now start build graph, command is:";
    say STDERR $cmd_build;
    system($cmd_build);
}

if(!defined $gfa_in) {
    my $cmd_gfa2vg = "$vg convert --gfa-out --gfa-in $gfa1 > $gfa2";
    say STDERR "Now covert gfa, command is:";
    say STDERR $cmd_gfa2vg;
    system($cmd_gfa2vg);
}




say STDERR "Now modify gfa";
&mod_vg_gfa($gfa2, $gfa3);

say STDERR "===========================";
say STDERR "===========================";
say STDERR "===========================";
say STDERR "Done. final gfa file is $gfa3";
say STDERR "===========================";
say STDERR "===========================";

exit;


sub mod_vg_gfa {
    my ($in, $out) = @_;
    my $I = open_in_fh($in);
    my $O = open_out_fh($out);
    while(<$I>) {
        if (/^P/) {
            /^P\t(\S+)\t/ or die;
            my $chr = $1;
            print $O $_ if $chr=~/^\d+$/;
        } else {
            print $O $_;
            next;
        }
    }
}

sub check_backbone {
    my ($infa) = @_;
    my $I = open_in_fh($infa);
    while(<$I>) {
        chomp;
        next unless $_;
        next unless /^>/;
        /^>(\S+)/ or die "Chromosome of Backbone error!";
        my $chr = $1;
        $chr=~/^\d+$/ or die "Chromosome of Backbone should be numeric only! : $_";
        $chr=~/^0/ and die "Chromosome of Backbone should not start with 0! : $_";
    }
}

sub check_dup_chrs {
    my ($backbone, @moregenomes) = @_;
    my %chrs;
    foreach my $fa ($backbone, @moregenomes) {
        my $I = open_in_fh($fa);
        while(<$I>) {
            chomp;
            next unless $_;
            next unless /^>/;
            /^>(\S+)/ or die "Chromosome name error for $fa line $. : $_";
            my $chr = $1;
            my $isint = $chr=~/^\d+$/;
            if ($isint==1 and $fa ne $backbone) {
                die "Chromosome name of none-backbone genome should NOT be numeric only";
            }
            if (exists $chrs{$chr}) {
                die "Dumplicated chromosome name! $fa : $chr";
            }
            $chrs{$chr}++;
        }
    }
}

