#
#===============================================================================
#
#         FILE: prase_vcf.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zeyu Zheng (LZU), zhengzy2014@lzu.edu.cn
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 12/01/2021 10:11:37 PM
#     REVISION: ---
#===============================================================================

package prase_vcf;
use strict;
use warnings;
use v5.24;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use zzIO;
use List::Util qw/max min/;
use Carp qw/confess carp/; # carp=warn;confess=die

require Exporter;

use vars qw(
  @ISA
  %EXPORT_TAGS
  @EXPORT_OK
  @EXPORT
);

@ISA = qw(Exporter);

%EXPORT_TAGS = (
    'all' => [
        qw(
            read_ref
            read_vcf_header
            find_union
            get_min_max_pos
            read_fai_file
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw(read_ref read_vcf_header find_union get_min_max_pos read_fai_file);



sub read_fai_file {
    my ($file) = @_;
    $file .= ".fai" if $file=~/\.(fa|fasta)$/;
    my $I = open_in_fh($file);
    my %ret;
    while(<$I>) {
        chomp;
        next unless $_;
        my @F = split(/\t/);
        $ret{$F[0]} = $F[1];
    }
    return \%ret;
}



sub read_ref {
    my ($in) = @_;
    my %ret;
    my $I = open_in_fh($in);
    my $chr_now;
    while(<$I>){
        chomp;
        next unless $_;
        if (/^>(\S+)/) {
            $chr_now = $1;
            next;
        }
        $ret{$chr_now} .= $_;
    }
    my $total_length;
    my $total_chrs = scalar(keys %ret);
    $total_length += length($ret{$_}) foreach keys %ret;
    my $stat = "read $total_chrs contig / total $total_length bp";
    return (\%ret, $stat);
}



sub read_vcf_header {
    my ($I, $O_VCF, $append) = @_;
    my @vcf_header;
    while(<$I>) {
        chomp;
        next unless $_;
        if (/^##/) {
            say $O_VCF $_ if $O_VCF;
            next;
        } elsif (/^#/) {
            @vcf_header = split(/\t/, $_);
            print $O_VCF $append if ($O_VCF and $append);
            say $O_VCF $_ if $O_VCF;
            last;
        } else {
            die;
        }
    }
    return(\@vcf_header);
}




sub find_union {
    my ($snps, $pos2snpid) = @_;
    my @unions;
    my %all_snpids; 
    $all_snpids{$_}=1 foreach (sort {$$snps{$a}[0]<=>$$snps{$b}[0]} keys %$snps);
    #die Dumper $pos2snpid;
    while (1) {
        say STDERR join(' ', "????? ", sort keys %all_snpids) if $main::debug;
        my ($snpid) = %all_snpids;  # one random pos which needs merge
        last unless defined $snpid;
        my $pos = $$snps{$snpid}[0] // last;
        last unless defined $pos;
        die "!? $pos ".Dumper $snps unless (exists $$pos2snpid{$pos});
        my @snpids = (sort {$$snps{$a}[0]<=>$$snps{$b}[0]} $$pos2snpid{$pos}->@*);

        exists $all_snpids{$_} and delete $all_snpids{$_} foreach @snpids;
        
        my %snpids;
        $snpids{$_}=1 foreach @snpids;
        my ($min_start, $max_end) = &get_min_max_pos(\@snpids, $snps);
        REDO1:
        say STDERR "REDO1: $min_start ~ $max_end " if $main::debug;
        foreach my $i ( List::Util::shuffle($min_start..$max_end) ) {
            my ($min_start_new, $max_end_new);
            my $is_update=0;
            if ( exists $$pos2snpid{$i} ) {
                $is_update=1;
                foreach my $snpid ($$pos2snpid{$i}->@*) {
                    exists $all_snpids{$snpid} and delete $all_snpids{$snpid};
                    $snpids{$snpid}++;
                }
                ($min_start_new, $max_end_new) = 
                    &get_min_max_pos($$pos2snpid{$i}, $snps, $min_start, $max_end);
            }
            if ($is_update==1 and ($min_start_new != $min_start or $max_end_new != $max_end)) {
                $min_start = $min_start_new;
                $max_end = $max_end_new;
                goto REDO1;
            }
        }
        foreach my $snpid (@snpids) {
            delete $all_snpids{$snpid} if exists $all_snpids{$snpid};
        }
        @snpids = sort {$a<=>$b} keys %snpids;
        push @unions, [\@snpids, $min_start, $max_end];
    }
    return(\@unions);
}



sub get_min_max_pos {
    my ($snpids, $snps, $old_min, $old_max) = @_;
    my @start_poss;
    my @end_poss;
    foreach my $snpid (@$snpids) {
        my ($start, $len) = $$snps{$snpid}->@*;
        my $end = $start + $len-1;
        push @start_poss, $start;
        push @end_poss, $end;
    }
    if (defined $old_min and defined $old_max) {
        return(min(@start_poss, $old_min), max(@end_poss, $old_max));
    } else {
        return(min(@start_poss), max(@end_poss));
    }
}


1;
