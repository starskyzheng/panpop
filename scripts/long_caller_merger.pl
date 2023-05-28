#!/usr/bin/env perl
use strict;
use warnings;
use v5.10;
use Getopt::Long;
use Data::Dumper;

use FindBin qw($Bin);
use lib "$Bin/../lib";
use zzIO;

my ($in, $out, $min_supporting, $samplename, $opt_help);

$min_supporting = 2;
my $force_merge_insertion = 0;
my $force_merge_deletion = 0;
my $force_merge_fold = 5;

sub help {
    my ($info) = @_;
    say STDERR $info if $info;
    say <<EOF;
usage: perl $0 -i <in.vcf> -o <out.vcf> -s <sample_name>
    -i, --in <in.vcf>           input vcf file
    -o, --out <out.vcf>         output vcf file
    -m, --min_supporting <int>  min supporting samples, default $min_supporting
    --force_merge_insertion     force merge insertion
    --force_merge_deletion      force merge deletion
    --force_merge_fold <int>    force merge fold for insertion or deletion, default $force_merge_fold
    -s, --sample <str>          sample name
    -h, --help                  print this help message
EOF
    exit(-1);
}

my $ARGVs = join ' ', $0, @ARGV;

GetOptions (
        'help|h!' => \$opt_help,
        'i|in=s' => \$in,
        'o|out=s' => \$out,
        'force_merge_insertion!' => \$force_merge_insertion,
        'force_merge_deletion!' => \$force_merge_deletion,
        'force_merge_fold=i' => \$force_merge_fold,
        'm|min_supporting=s' => \$min_supporting,
        's|sample=s' => \$samplename,
);

&help("force_merge_fold < 1: $force_merge_fold") if $force_merge_fold < 1;
&help() if $opt_help;
&help() unless defined $in;
&help() unless defined $out;

my $I = open_in_fh($in);
my $O = open_out_fh($out);
my @header;

while(my $line = <$I>) { # vcf
    chomp $line;
    if($line=~/^#/) {
        if($line=~/^##/) {
            say $O $line;
            next;
        }
        say $O q|##INFO=<ID=SUPP,Number=1,Type=String,Description="Number of samples supporting the variant">|;
        say $O qq+##CommandLine="$ARGVs"+;
        @header = split /\t/, $line;
        say $O join("\t", @header[0..8], $samplename);
        next;
    }
    my @F = split /\t/, $line;
    $F[4] = [split(/,/, $F[4])];
    my ($chr, $pos, $svid, $ref, $alts) = @F[0,1,2,3,4];
    #next unless $pos == 4474552;
    my $is_ins = $force_merge_insertion==1 ? &cal_is_ins(\@F) : 0;
    my $is_del = $force_merge_deletion==1 ? &cal_is_del(\@F) : 0;
    my ($max_alle, $new_gt, $supp) = &cal_alle_freq(\@F, $is_ins, $is_del);
    #say $O join "\t", $F[0], $F[1], $max_alle if defined $max_alle;
    say $O &rebuild_line(\@F, $max_alle, $new_gt, $supp) if defined $max_alle;
}

exit;

sub cal_is_del {
    my ($F) = @_;
    my $ref = $$F[3];
    my $alts = $$F[4];
    my $min_len = length($ref)/$force_merge_fold;
    $min_len = 1 if $min_len < 1;
    foreach my $alt (@$alts) {
        if(length($alt) > $min_len) {
            return 0;
        }
    }
    return 1;
}

sub cal_is_ins {
    my ($F) = @_;
    my $ref = $$F[3];
    my $alts = $$F[4];
    my $min_len = $force_merge_fold * length($ref);
    foreach my $alt (@$alts) {
        if(length($alt) < $min_len) {
            return 0;
        }
    }
    return 1;
}

sub rebuild_line {
    my ($F, $alle, $new_gt, $supp) = @_;
    my @newline = $F->@[0..8];
    my $alts = $$F[4];
    my $newalt = $$alts[$alle-1];
    $newline[4] = $newalt;
    if($newline[7] eq '.') {
        $newline[7] = "SUPP=$supp";
    } else {
        $newline[7] .= ";SUPP=$supp";
    }
    $newline[8] = 'GT';
    $newline[9] = $new_gt;
    return join "\t", @newline;
}

sub cal_max_alle {
    my ($alle_freq, $hetero_freq) = @_;
    my @alles = sort{$$alle_freq{$b} <=> $$alle_freq{$a}} keys %{$alle_freq};
    my $max_alle = $alles[0] // return undef;
    my $max_alle_freq = $$alle_freq{$max_alle};
    if($max_alle_freq < $min_supporting) {
        return undef;
    }
    if(scalar(@alles) >= 2) {
        my $max_alle2 = $alles[1];
        my $max2_alle_freq = $$alle_freq{$max_alle2};
        if($max2_alle_freq == $max_alle_freq) {
            $$hetero_freq{$max_alle} //= 0;
            $$hetero_freq{$max_alle2} //= 0;            
            if($$hetero_freq{$max_alle} < $$hetero_freq{$max_alle2}) {
                # 频率相同，输出纯合多的
                return $max_alle;
            } elsif($$hetero_freq{$max_alle} > $$hetero_freq{$max_alle2}) {
                return $max_alle2;
            } else {
                # ??? select first
                return $max_alle;
            }
        }
    }
    return $max_alle;
}

sub cal_alle_freq {
    my ($F, $is_ins, $is_del) = @_;
    my $is_force = $is_ins || $is_del;
    my %alle_freq;
    my %heteros;
    my $not_missing_ref = 0;
    foreach my $i (9..$#{$F}) {
        my $gt = $F->[$i];
        if($gt=~/^\./) {
            next;
        }
        
        $gt=~m#(\d+)/(\d+)# or die;
        my ($gt1, $gt2) = ($1, $2);
        $not_missing_ref++ if $gt1!=0 or $gt2!=0;
        if($gt1!=0 and $gt2!=0 and $gt1!=$gt2) { # 1/2
            #die "Error: $gt";
            $alle_freq{$gt1} += 1;
            $alle_freq{$gt2} += 1;
            $heteros{$gt1}++;
            $heteros{$gt2}++;
        } elsif($gt1!=$gt2) { # 0/1 or 1/0 or 1/2 or 0/2
            foreach my $gt ($gt1, $gt2) {
                next if $gt==0;
                $alle_freq{$gt}++;
                $heteros{$gt}++;
            }
        } elsif($gt1==$gt2 and $gt1!=0) { # 1/1 or 2/2
            $alle_freq{$gt1}++;
        }
    }
    my $max_alle = &cal_max_alle(\%alle_freq, \%heteros);
    if(! defined $max_alle) {
        if($is_force==1 and $not_missing_ref >= $min_supporting) {
            my $longest_allei;
            if($is_ins==1) {
                $longest_allei = 1 + &get_longest_alle_i($F->[4]);
            } elsif ($is_del==1) {
                $longest_allei = 1 + &get_shortest_alle_i($F->[3]);
            } else {
                die;
            }
            my $gt_now = '1/1';
            if($alle_freq{$longest_allei}==1 and exists $heteros{$longest_allei} and $heteros{$longest_allei}==1) {
                $gt_now = '0/1';
            }
            return($longest_allei, $gt_now, $not_missing_ref);
        } else {
            return(undef);
        }
    }
    my $supp = $alle_freq{$max_alle};
    my $hetero_alle = $heteros{$max_alle} // 0;
    my $hetero_rate = $hetero_alle / $supp;
    my $gt;
    if($hetero_rate > 0.5) {
        $gt = '0/1';
    } else {
        $gt = '1/1';
    }
    if($is_force==1) {
        return($max_alle, $gt, $min_supporting);
    } else {
        return($max_alle, $gt, $supp);
    }
}

sub get_longest_alle_i {
    my ($alles) = @_;
    my $max_len = 0;
    my $max_i = 0;
    foreach my $i (0..$#{$alles}) {
        my $len = length($$alles[$i]);
        if($len > $max_len) {
            $max_len = $len;
            $max_i = $i;
        }
    }
    return $max_i;
}

sub get_shortest_alle_i {
    my ($alles) = @_;
    my $min_len = 100000000;
    my $min_i = 0;
    foreach my $i (0..$#{$alles}) {
        my $len = length($$alles[$i]);
        if($len < $min_len) {
            $min_len = $len;
            $min_i = $i;
        }
    }
    return $min_i;
}




