#!/usr/bin/env perl
use warnings;
use strict;
use v5.24;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Data::Dumper;
use Carp qw/confess carp/; # carp=warn;confess=die
use MCE;
use MCE::Channel;
use MCE::Flow;
use MCE::Candy;
use Getopt::Long;
use zzIO;
no warnings 'qw';

sub help {
    print <<EOF;
Usage: $0  [options]
Options:
    -i | --invcf <file>	         Input vcf file
    -o | --outvcf <file>         Output file, default: -
    -p | --threads <int>         Thread number, default: auto
    -F | --force_recode	<bool>   Force recode for each line, default: False
EOF
}

my ($in_vcf, $out_vcf, $threads, $force_recode, $ignore_dp);# = @ARGV;

my $header_append_info = "##merge_vcf_samepos=perl $0 @ARGV";

GetOptions (
    'i|invcf=s' => \$in_vcf,
    'o|outvcf=s' => \$out_vcf,
    'p|threads=i' => \$threads,
    'F|force_recode' => \$force_recode,
    'h|help' => sub { help(); exit(0); },
    'ignore_dp!' => \$ignore_dp,
);

&help() unless $in_vcf and -e $in_vcf and $out_vcf;


$threads //= 1;
$force_recode //= 0;

my $I = open_in_fh($in_vcf);
my $O = open_out_fh($out_vcf, 8);
my $OM = open_out_fh("$out_vcf.allmising", 8);
say $OM join("\t", qw/#chr pos ref_len/);
#my $infos = 'GT:DP';

my @vcf_header;
while(<$I>) {
    chomp;
    next unless $_;
    if (/^#/) {
        if (/^##/) {
            say $O $_; ###########
            next;
        }
        @vcf_header = split(/\t/);
        say $O $header_append_info;
        say $O $_;
        last;
    }
    die;
}

my $idi_max = scalar(@vcf_header)-1;
my $id_count = $idi_max - 8;
my $queue = MCE::Channel->new(impl => 'Mutex', mp => 0);


MCE::Flow->init(
    max_workers => [ 1, $threads ],
    task_name   => [ 'mce_genertor', 'consumer' ],
    gather => MCE::Candy::out_iter_fh( $O ),
    task_end    => sub {
        my ($mce, $task_id, $task_name) = @_;
        if ( $task_name eq 'mce_genertor' ) {
            $queue->end;
        }
    }
);

mce_flow \&mce_genertor, \&consumer;
MCE::Flow->finish;

exit;


sub consumer {
   # receive items
   while ( my $parms = $queue->dequeue() ) {
        my ( $lineid, $lines ) = @$parms;
        last unless ($lines and @$lines);
        my $new_lines = &print_lens($lines);
        MCE->gather( $lineid , $new_lines );
   }
}

sub mce_genertor {
    my $chr_old;
    my $pos_old;
    my @dat_old;
    my $lineid=1;
    while(my $line = <$I>) {
        #say STDERR $lineid;
        chomp $line;
        next unless $line;
        $line=~m#^(\S+)\t(\S+)\t# or die;
        #my @F = split(/\t/, $line);
        #my ($chr, $pos) = ($F[0], $F[1]);
        my $chr = $1;
        my $pos = $2;
        $chr_old //= $chr;
        $pos_old //= $pos;
        if ($chr_old ne $chr or $pos_old != $pos) {
            $queue->enqueue( [$lineid, \@dat_old] );
            $chr_old = $chr;
            $pos_old = $pos;
            @dat_old = ();
            $lineid++;
            redo;
        }
        say STDERR "Now line $. at $chr : $pos" if ($. % 1e5 == 0);
        #push @dat_old, \@F;
        push @dat_old, $line;
    }
    $queue->enqueue( [$lineid, \@dat_old] ) if @dat_old;
    $queue->end();
}



sub print_lens {
    my ($lines) = @_;
    my %id2seqs;
    my %refs;
    my %infos;
    if ($force_recode==1 and scalar(@$lines)==1) { # no need to merge lines
        return $$lines[0] . "\n";
    }
    my %id2types;
    LINE:foreach my $line (@$lines) {
        my @line = split(/\t/, $line);
        my $ref = $line[3];
        die unless defined $ref and $ref;
        #die unless $line[8]=~/^$infos/;
        my ($DP_num, $GQ_num, $MAD_num, $DPSOURCE_num, $MDP_num);
        my @info=split(/:/,$line[8]);
        for (my $ii=0;$ii<@info;$ii++){
            if ($info[$ii] eq 'DP'){
                $DP_num = $ii ;
                next;
            } elsif($info[$ii] eq 'GQ'){
                $GQ_num = $ii;
                next;
            } elsif($info[$ii] eq 'MAD'){
                $MAD_num = $ii;
                next;
            } elsif ($info[$ii] eq 'DPSOURCE'){
                $DPSOURCE_num = $ii;
                next;
            }
        }
        # die unless defined $DP_num;# and defined $GQ_num and defined $MAD_num;
        my $ref_len = length($ref);
        $infos{$ref_len} = [ $line[0], $line[1] ];
        $refs{$ref_len} //= $ref;
        die if $refs{$ref_len} ne $ref;
        my @ref_alts = ($ref, split(/,/, $line[4]));
        for(my $idi=9; $idi<=$idi_max; $idi++) {
            #say STDERR "$idi $vcf_header[$idi] " ;
            my @F = split(/:/, $line[$idi]);
            $F[0]=~m#^([.\d]+)[/|]([.\d]+)$# or die $F[0];
            my ($a1, $a2) = ($1, $2);
            my $dp = defined $DP_num ? $F[$DP_num] : '.';
            my $dpsource = defined $DPSOURCE_num ? $F[$DPSOURCE_num] : undef;
            my $gq = defined $GQ_num ? $F[$GQ_num] : '.';
            my $mad = defined $MAD_num ? $F[$MAD_num] : '.';
            if ( defined $dpsource and $dpsource eq 'append' ) {
                $mad //= $dp;
            }
            if($ignore_dp==1) {
                $dp = 9;
            }
            if ($a1 eq '.' or $a2 eq '.' or 
                    $dp eq '.' or 
                    $dp==0 #or 
                    #$gq eq '.' or 
                    #$mad eq '.'
                    ) {
                #$id2seqs{$ref_len}{$idi} = [ -1 ];
                next;
            }
            my $mut_now = -1;
            my $mut_old = -1;
            if (exists $id2seqs{$ref_len}{$idi}) {
                $mut_old = $id2seqs{$ref_len}{$idi}[0] // die;
            }
            if($a1==0 and $a2==0) {
                # ref
                $mut_now = 0;
            } else {
                # alt
                $mut_now = 1;
            }
            if ( $mut_now==1 and $mut_old==1) {
                # ERROR! same pos with two mut
                my $old_seq1 = $id2seqs{$ref_len}{$idi}[1];
                my $old_seq2 = $id2seqs{$ref_len}{$idi}[2];
                carp " exists! {$ref_len}{$idi} @{$id2seqs{$ref_len}{$idi}} $vcf_header[$idi]
                    !!! line !!! :  @$line ";
            } elsif ($mut_now == -1 and $mut_old != -1) {
                # missing now but not missing before
                next;
            } elsif ($mut_now == -1 and $mut_old == -1) {
                # both missing
                #$id2seqs{$ref_len}{$idi} = [ -1 ];
                next;
            } elsif ( ($mut_now != -1 and $mut_old == -1) or # not missing replace missing
                    ($mut_now == 1 and $mut_old == 0) ) { # mut replace ref
                my $seq1 = $ref_alts[$a1];
                my $seq2 = $ref_alts[$a2];
                $id2seqs{$ref_len}{$idi} = [ $mut_now, $seq1, $seq2, $dp, $gq, $mad ];
            } elsif($mut_now == 0 and $mut_old == 0 ) {
                # both not missing and ref, use maxmal dp
                if ($dp > $id2seqs{$ref_len}{$idi}[3]) {
                    $id2seqs{$ref_len}{$idi} = [ $mut_now, $ref, $ref, $dp, $gq, $mad ];
                } else {
                    next;
                }
            } elsif($mut_now == 0 and $mut_old == 1 ) {
                # ref not replace alt
                next;
            } else {die;}
        }
    }
    #die Dumper \%id2seqs;

    my @final_prints;
    PROC_BY_LEN:foreach my $ref_len (sort {$a<=>$b} keys %id2seqs) {
        my $seqs = $id2seqs{$ref_len};
        my $ref = $refs{$ref_len};
        my %new_ids2alles;
        my %new_ref_alts = ($ref=>0);
        my $miss_count = 0;
        for(my $idi=9; $idi<=$idi_max; $idi++) {
            unless ( exists $$seqs{$idi} and defined $$seqs{$idi} ) { 
                $new_ids2alles{$idi} = './.:.:.:.';
                $miss_count++;
                next;
            }
            #my $mut_type = $$seqs{$idi}[0];
            my $dp = $$seqs{$idi}[3] // die;
            my $gq = $$seqs{$idi}[4] // die;
            my $mad = $$seqs{$idi}[5] // die;
            my @idi_alles;
            foreach my $alle_seq ( $$seqs{$idi}[1], $$seqs{$idi}[2] ) {
                if (exists $new_ref_alts{$alle_seq}) {
                    my $alle = $new_ref_alts{$alle_seq};
                    push @idi_alles, $alle;
                } else {
                    my $now_new_alles_count = scalar(keys %new_ref_alts);
                    $new_ref_alts{$alle_seq} = $now_new_alles_count;
                    push @idi_alles, $now_new_alles_count;
                }
            }
            $new_ids2alles{$idi} = join('/', @idi_alles) . ":$dp:$gq:$mad";
        }
        if ($id_count <= $miss_count) {
            # all missing
            say $OM join("\t", $infos{$ref_len}[0], $infos{$ref_len}[1], length($ref));
            next PROC_BY_LEN;
        }
        my @ref_alts_seqs = sort {$new_ref_alts{$a}<=>$new_ref_alts{$b}} keys %new_ref_alts;
        shift @ref_alts_seqs; # remove ref
        if (! @ref_alts_seqs) {
            # no alt!
            @ref_alts_seqs = ('.');
        }
        my @prints = ($infos{$ref_len}[0], $infos{$ref_len}[1], '.', $ref, 
            join(',', @ref_alts_seqs), 
            '.', 'PASS', '.', "GT:DP:GQ:MAD");
        push @prints, $new_ids2alles{$_} foreach (9..$idi_max);
        #say $O join "\t", @prints;
        push @final_prints, join("\t", @prints);
    }
    #print STDERR join("\n", @final_prints)."\n";
    return join("\n", @final_prints)."\n" if @final_prints;
    return '';
}

