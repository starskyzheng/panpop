#!/usr/bin/env perl
#===============================================================================
#
#         FILE: sv2snp_indel.pl
#
#        USAGE: ./sv2snp_indel.pl  
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
use IPC::Open2;
use Tie::CharArray;
use MCE::Flow;
use MCE::Candy;
use Getopt::Long;
use MCE::Channel;
use Coro::Generator;
use File::Temp;
use List::Util qw/max min/;

use prase_vcf qw/read_ref read_vcf_header/;
use read_config qw/read_config_yaml/;
use realign_alts qw/process_alts/;


my ($infile, $outfile, $opt_help, $ref_fasta_file);
our $tmp_dir_def = '/run/user/' . `id -u`; chomp $tmp_dir_def;
our $debug = 0;
our $verb = 0;
our $align_level = 1;
my $is_aug = 0;
my $ext_bp_max = 10;
my $ext_bp_min = 1;
my $print_all = 0;
our $threads = 1; # default threads
our $tmp_dir = $tmp_dir_def;

#my @SKIPCHR = qw/1 2 3 5 6/;

sub usage {
    my $msg = shift;
    print STDERR "$msg\n" if $msg;
    print STDERR <<"EOF";
USAGE: perl $0 [options]
options:
    -i | --in_vcf <file>           Input vcf file
    -o | --out_vcf <file>          Output vcf file
    -r | --ref_fasta_file <file>   Reference fasta file
    -t | --threads <int>           Number of threads (default: $threads)
    -E | --ext_bp_max <int>        Max extension bp (default: $ext_bp_max)
    -e | --ext_bp_min <int>        Min extension bp (default: $ext_bp_min)
    --tmpdir <dir>                 Temprory directory (default: $tmp_dir)
    -h | --help                    Print this help
    --verb <bool>
    --all <bool>                   Print lines even no mutation.
    --level <1|2>
EOF
    exit(1);
}


GetOptions (
        'h|help!' => \$opt_help,
        'i|in_vcf=s' => \$infile,
        'o|out_vcf=s' => \$outfile,
        'debug!' => \$debug,
        'verb!' => \$verb,
        't|threads=i' => \$threads,
        'r|ref_fasta_file=s' => \$ref_fasta_file,
        'aug!' => \$is_aug,
        'tmpdir=s' => \$tmp_dir,
        'E|ext_bp_max=i' => \$ext_bp_max,
        'e|ext_bp_min=i' => \$ext_bp_min,
        'all!' => \$print_all,
        'level=i' => \$align_level,
);


&usage() if $opt_help or !$infile or !$outfile or !$ref_fasta_file;

$tmp_dir = $tmp_dir_def if $tmp_dir eq 'Default';
if (! -e $tmp_dir) {
    mkdir $tmp_dir or die "cannot create tmpdir: $tmp_dir!";
}

our $config = read_config_yaml("$Bin/../config.yaml");
realign_alts::init();

my $I_VCF = open_in_fh($infile);
my $O_VCF = open_out_fh($outfile);
my $queue = MCE::Channel->new( impl => 'Mutex', mp => 0 );

say STDERR "Now reading ref fasta : ";
my ($REF_SEQS, $REF_SEQS_STAT, $NONREF_SEQS);
if ($is_aug==0) {
    ($REF_SEQS, $REF_SEQS_STAT, $NONREF_SEQS) = &read_ref_noaug($ref_fasta_file);
} elsif ($is_aug==1) {
    ($REF_SEQS, $REF_SEQS_STAT) = &read_ref($ref_fasta_file);
} else {
    die "?? Error";
}
say STDERR "Done reading ref fasta $REF_SEQS_STAT";

my $VCF_APPEND = <<'EOF';
##INFO=<ID=ZORIPOS,Number=1,Type=Integer,Description="original start position before Realign of Variant">
##INFO=<ID=ZORIEND,Number=1,Type=Integer,Description="original end position before Realign of Variant">
##INFO=<ID=ZLEN,Number=1,Type=Integer,Description="original length before Realign of Variant">
EOF

my $vcf_header = &read_vcf_header($I_VCF, $O_VCF, $VCF_APPEND);
my $idi_max = scalar(@$vcf_header)-1;

my $zzgenerator = generator {
    my $I = $I_VCF;
    my @lines;
    my $min_start;
    my $max_end;
    my $end_real;
    my $chr_old;
    while(<$I>) {
        chomp;
        next if /^#/;
        next unless $_;
        my $snpid = $.;
        #my $snpid = "$F[0]:$F[1]"; # may not unique
        my @F = split(/\t/, $_);
        next if $F[4] eq '';
        my $chr = $F[0];
        my $pos = $F[1];
        my $ref_seq = $F[3];
        my @alts = split(/,/, $F[4]);
        $F[4] = \@alts;
        # next if $pos < 14999102; ######### debug
        #next if exists $SKIPCHR{$F[0]};
        unless (defined $chr_old) { # init
            $chr_old = $chr;
            $min_start = $pos;
            ($max_end, $end_real) = &cal_max_ext_range($ref_seq, \@alts, $pos, -1);
            redo;
        }
        
        if ($chr_old ne $chr or $max_end+1 < $pos) {
            #print(STDERR "2!! $chr_old\t$min_start\t$max_end\t$pos\n");
            die unless defined $min_start;
            die unless defined $end_real;
            yield([\@lines, $chr_old, $min_start, $end_real]);
            @lines=();
            if ($chr_old ne $chr) {
                $chr_old = undef;
                $min_start = undef;
                $max_end = undef;
                $end_real = undef;
            } else {
                $max_end = $pos;
                $end_real = $pos;
                $min_start = $pos;
            }
            redo;
        }
        $min_start //= $pos;
        ($max_end, $end_real) = &cal_max_ext_range($ref_seq, \@alts, $pos, $max_end);
        push @lines, \@F;
    }
    close $I;
    yield([\@lines, $chr_old, $min_start, $end_real]) if @lines;
    yield(undef);
};

sub zz_MCE_OUT {
    my $fh = $_[0];
    $fh = \$_[0] if (!ref $fh && ref \$_[0]);
    if ($fh->can('print')) {
        return sub {
            my ($chunk_id, $data) = @_;
            $fh->print($data);
        };
    } else {
        return sub {
            my ($chunk_id, $data) = @_;
            print {$fh} $data;
        };
    }
}

MCE::Flow->init(
    chunk_size => 1,
    max_workers => [ 1, $threads ],
    task_name => [ 'producer', 'consumer' ],
    gather => &zz_MCE_OUT($O_VCF),
    #gather => MCE::Candy::out_iter_fh($O_VCF), # output sorted vcf but consume too much memory
    task_end => sub {
        my ( $mce, $task_id, $task_name ) = @_;
        if ( $task_name eq 'producer' ) {
            $queue->end;
        }
    }
);

mce_flow \&zz_mce_producer, \&mce_run;

MCE::Flow->finish;
close $O_VCF;

say STDERR "All Done";


exit;


sub cal_max_ext_range {
    my ($ref, $atls, $pos, $max_end) = @_;
    my $reflen = length($ref);
    my $alt_max_len = &cal_max_length_array($atls);
    my $end_real = $pos + $reflen - 1; # not consider '*'
    my $end = $end_real;
    my $max_len = $alt_max_len > $reflen ? $alt_max_len : $reflen;
    my $ext = 0;
    if ($max_len <= 1) { # snp
        $ext = 1;
    } else { # sv/indel
        $ext = int($max_len * 0.2)+1;
    }
    if ($ext > $ext_bp_max) {
        $ext = $ext_bp_max;
    } elsif ($ext < $ext_bp_min) {
        $ext = $ext_bp_min;
    }
    $end += $ext;
    if ($end > $max_end) {
        return ($end, $end_real);
    } else {
        return ($max_end, $end_real);
    }
}

sub cal_max_length_array {
    my $A = shift;
    my $max_len = 0;
    foreach my $a (@$A) {
        my $len = length($a);
        $max_len = $len if $len > $max_len;
    }
    return $max_len;
}


sub mce_run {
    while ( my ( $chunk_id, $parms ) = $queue->dequeue(2) ) {
        last unless defined $parms;
        my ($lines, $chr, $start, $end) = @$parms;
        #say STDERR Dumper $parms;
        last unless defined $lines;
        if ($chunk_id % 1000==0) {
            my ($oneline) = @$lines;
            my $chr = $$oneline[0] or die Dumper $oneline;
            my $pos = $$oneline[1] or die Dumper $oneline;
            say STDERR "Now: $chr:$pos";#, snpids: ", join ',', sort {$a<=>$b} keys %$snps;
        }
        my $results = &process_lines_new($lines, $chr, $start, $end);

        if ($results) {
            MCE->gather( $chunk_id, join('', $results));
        } else {
            MCE->gather( $chunk_id, '');
        }
    }
}

sub zz_mce_producer {
    my $iline=1;
    while ( my ($parms) = $zzgenerator->() ) {
        unless (defined $parms) { # finsh read vcf
            last;
        }
        $queue->enqueue( $iline, $parms ) ;
        $iline++;
    }
}



sub process_lines_new {
    my ($lines, $chr, $win_start, $win_end) = @_;
    my @lines_split;
    foreach my $F (@$lines) {
        #my @F = split(/\t/, $line);
        my $alts = $$F[4];
        for(my $i=0; $i<scalar(@$alts); $i++) {
            $$alts[$i]='' if $$alts[$i] eq '*';
        }
        #$$F[4] = \@alts;
        for (my $idi = 9; $idi <= $idi_max; $idi++) {
            if (! defined $$F[$idi]) {
                my $l = scalar(@$F);
                die "Error: len: $l; idi:$idi; line: @$F ; ";
            }
            if ($$F[$idi]=~m#^(\d+)/(\d+)#) {
                $$F[$idi] = [$1, $2];
            } else {
                $$F[$idi] = undef;
            }
        }
    }
    my ($new_ref_alts, $new_ids2alles) = &rebuild_cons_seqs($lines, $chr, $win_start, $win_end);
    my $all_lines = &process_line_new($new_ref_alts, $new_ids2alles, $chr, $win_start);
    return undef unless defined $all_lines;
    my $final='';
    my $maxi = scalar(@$all_lines)-1;
    for(my $i=0; $i<=$maxi; $i++) {
        $final .= join("\t", $$all_lines[$i]->@*) . "\n";
    }
    return($final);
}


sub process_line_new {
    my ($ref_alts, $ids2alles, $chr, $win_start) = @_;
    my $is_snp=1;
    my $alt_max_length = -1;
    my $alt_min_length = 999999;
    foreach my $alt (@$ref_alts) {
        my $len = length ($alt);
        if ( $len != 1 ) {
            $is_snp=0;
        }
        $alt_max_length = $len if $alt_max_length < $len;
        $alt_min_length = $len if $alt_min_length > $len;
    }
    if ($is_snp==1 or $alt_min_length<=1) { ########## SNP
        my $retline = &gen_lines_before_process($ids2alles, $ref_alts, 
                                                    $chr, $win_start);
        if ($retline) {
            return([$retline])
        } else {return undef}
    }
    my $ref_len = length($$ref_alts[0]);
    my $max_alts = scalar(@$ref_alts)-1;
    #say STDERR $F[1];
    if ($max_alts==0) {
        say STDERR "No alts at $chr $win_start" if $debug;
        return undef;
    }
    my ($muts) = &process_alts($ref_alts, $max_alts, $alt_max_length);
    if (! defined $muts) {
        say STDERR "Fail to cal aln !! : ";
        return undef;
    }
    my $new_lines = &gen_lines($muts, $ids2alles, $max_alts, $ref_len, $chr, $win_start);
    return($new_lines);
}

sub gen_lines_before_process {
    my ($ids2alles, $ref_alts, $chr, $win_start) = @_;
    for(my $i=0; $i<scalar(@$ref_alts); $i++) {
        $$ref_alts[$i]='*' if $$ref_alts[$i] eq '';
    }
    my $ref_seq = shift @$ref_alts // die;
    my $new_alts_join = join(',', @$ref_alts);
    return undef if ($print_all==0 and ! @$ref_alts); # do not print if no alts
    $new_alts_join = '.' if $new_alts_join eq ''; # no alts
    my @new_line = ($chr, $win_start, '.', $ref_seq, $new_alts_join,
            '.', 'PASS', '.', 'GT');
    foreach my $idi (sort {$a<=>$b} keys %$ids2alles) {
        my $alle_array = $$ids2alles{$idi};
        my $alle;
        if (! defined($alle_array)) {
            $alle = './.';
        } else {
            $alle = join('/', @$alle_array);
        }
        push @new_line, $alle;
    }
    return(\@new_line);
}

sub get_ref_substr_seq {
    my ($chr, $win_start, $win_end) = @_;
    die unless defined $win_end;
    die unless defined $win_start;
    if ( exists $$REF_SEQS{$chr} ){
        my $ref = substr($$REF_SEQS{$chr}, $win_start-1, $win_end - $win_start + 1);
        return $ref;
    } elsif ( $is_aug==0 and exists $$NONREF_SEQS{$chr} ) {
        my $dist_min=999999999;
        my $dist_min_ref_pos;
        foreach my $ref_pos_t (keys $$NONREF_SEQS{$chr}->%*) {
            my $dist = $win_start - $ref_pos_t;
            next if $dist < 0;
            if ($dist < $dist_min) {
                $dist_min = $dist;
                $dist_min_ref_pos = $ref_pos_t;
            }
        }
        die "nonref: $chr:$win_start-$win_end not exists in fasts! $ref_fasta_file"
            unless defined $dist_min_ref_pos;
        my $ref_all = $$NONREF_SEQS{$chr}{$dist_min_ref_pos};
        my $len_all = length $ref_all;
        my $ref_pos_end = $dist_min_ref_pos + $len_all - 1;
        if ($ref_pos_end < $win_end) {
            die "????";
        }
        die unless defined $win_end;
        die unless defined $win_start;
        #say STDERR " $chr $win_start-$dist_min_ref_pos, $win_end - $win_start + 1";
        my $ref = substr($ref_all, $win_start-$dist_min_ref_pos, $win_end - $win_start + 1);
        return $ref;
    } else {
        die "chr: $chr not exists in fasta: $ref_fasta_file";
    }
}

sub rebuild_cons_seqs {
    my ($lines, $chr, $win_start, $win_end) = @_;
    my $ref = &get_ref_substr_seq($chr, $win_start, $win_end);
    say STDERR "$chr $win_start-$win_end $ref " . length($ref) if $debug;
    my $iline_max = scalar(@$lines)-1;
    my @phase_issue_ids;
    my %new_ref_alts = ($ref=>0);
    my %new_ids2alles;
    say STDERR "Now start rebuild_cons_seqs" if $debug;
    rebuild_cons_seqs_IDI:for (my $idi = 9; $idi <= $idi_max; $idi++) {
        my %sites_status; # -1:miss   0:ref    1:alt
        my $a1 = $ref;
        my $a2 = $ref;
        my $last_start = 99999999999999;
        my $alt_svs=0;
        for(my $iline=$iline_max; $iline>=0; $iline--) { # from last to first
            my $line = $$lines[$iline];
            my $sv_start = $$line[1];
            die if $last_start < $sv_start;
            $last_start = $sv_start;
            my $ref_len = length($$line[3]);
            my $sv_end = $sv_start + $ref_len-1;

            my $alles = $$lines[$iline][$idi];
            if (!defined $alles) { # miss
                for(my $p=$sv_start; $p<=$sv_end; $p++) {
                    $sites_status{$p} = -1 if ( !exists $sites_status{$p} );
                }
                next;
            } elsif ($$alles[0]==0 and $$alles[1]==0) { # ref
                for(my $p=$sv_start; $p<=$sv_end; $p++) {
                    if ( (!exists $sites_status{$p}) or $sites_status{$p} == -1 ) {
                        $sites_status{$p} = 0;
                    }
                }
                next;
            }
            ## alt
            $alt_svs++;
            for(my $p=$sv_start; $p<=$sv_end; $p++) {
                if (exists $sites_status{$p} and $sites_status{$p} == 1) {
                    die "DIE! chr:$chr, wins:$win_start, wine:$win_end, pos:$p, idi:$idi, line:$iline, alles:$$alles[0], $$alles[1]";
                } else {
                    $sites_status{$p}=1;
                }
            }
            if ($$alles[0] != 0) {
                my $alt = $$line[4][ $$alles[0]-1 ];
                say STDERR "$a1, $chr : $sv_start, $win_start~$win_end, $ref_len, $alt" if $debug;
                substr($a1, $sv_start-$win_start, $ref_len, $alt);
            }
            if ($$alles[1] != 0) {
                my $alt = $$line[4][ $$alles[1]-1 ];
                substr($a2, $sv_start-$win_start, $ref_len, $alt);
            }
        }
        my ($count_all, $count_miss) = (0,0);
        foreach my $p (keys %sites_status) {
            $count_all++;
            if ($sites_status{$p} == -1) { # has miss
                $count_miss++;
            }
        }
        if ($count_all == $count_miss) {
            # all missing
            $new_ids2alles{$idi} = undef;
            next rebuild_cons_seqs_IDI;
        }

        foreach my $alle_seq ($a1, $a2) {
            if (exists $new_ref_alts{$alle_seq}) {
                my $alle = $new_ref_alts{$alle_seq};
                push $new_ids2alles{$idi}->@*, $alle;
            } else {
                my $now_new_alles_count = scalar(keys %new_ref_alts);
                $new_ref_alts{$alle_seq} = $now_new_alles_count;
                push $new_ids2alles{$idi}->@*, $now_new_alles_count;
            }
        }
        push @phase_issue_ids, $idi if ($alt_svs>1);
    }
    say STDERR "Debug: phase: @phase_issue_ids" if ($debug>=1 and @phase_issue_ids);
    my @new_ref_alts;
    push @new_ref_alts, $_ foreach sort {$new_ref_alts{$a}<=>$new_ref_alts{$b}} keys %new_ref_alts;
    return(\@new_ref_alts, \%new_ids2alles);
}


sub string_array_max_len {
    my ($a) = @_;
    my $maxlen = -1;
    my $maxi = scalar(@$a)-1;
    for(my $i = 0; $i <= $maxi; $i++) {
        my $len = length($$a[$i]);
        $maxlen = $len if $len>$maxlen;
    }
    return($maxlen);
}


sub gen_lines {
    my ($muts, $ids2alles, $max_alts, $seqlen, $chr, $old_pos) = @_;
    my @new_lines;
    my $old_pos_end = $old_pos + $seqlen -1;
    my @oldline_ori = (('.') x 8);
    $oldline_ori[7] = "ZORIPOS=$old_pos;ZORIEND=$old_pos_end;ZLEN=$seqlen";
    my $new_info = 'GT';
    foreach my $nmut (sort {$a <=> $b} keys %$muts) {
        my $muts_now = $$muts{$nmut};
        my @alt_seqs = $muts_now->@[0..$max_alts];
        my ($remap_alts, $new_ref, $new_alts) = &thin_alts(\@alt_seqs);
        my $ref_seq = shift @alt_seqs;
        next if ($print_all==0 and ! @$new_alts); # do not print if no alts
        my $new_alts_join = join(',', @$new_alts);
        $new_alts_join = '.' if $new_alts_join eq ''; # no alts
        my @line = ($chr, $old_pos+$nmut, '.', $new_ref, $new_alts_join, 
            @oldline_ori[5..7], $new_info);

        for (my $i = 9; $i<=$idi_max; $i++) {
            if (! defined $$ids2alles{$i}) {
                push @line, "./.";
                next;
            }
            my ($g1, $g2) = $$ids2alles{$i}->@*;
            my @new_gts;
            foreach my $g($g1, $g2) {
                if (exists $$remap_alts{$g}) {
                    $g = $$remap_alts{$g};
                }
                push @new_gts, $g;
            }
            push @line, "$new_gts[0]/$new_gts[1]";# . "$oldline_ori[$i]";
        }
        say STDERR "@line" if $debug;
        push @new_lines, [@line];
    }
    return(\@new_lines);
}

sub thin_alts {
    my ($alts) = @_;
    my $max_alts = scalar(@$alts)-1;
    my %seq2i;
    my %seq2remapi;
    my %remapi=();
    my @new_alts;
    my $now_maxalle=0;
    for (my $i = 0; $i <= $max_alts ; $i++) {
        my $seq = $$alts[$i];
        if (exists $seq2remapi{$seq}) {
            $remapi{$i} = $seq2remapi{$seq};
            next;
        }
        $remapi{$i} = $now_maxalle;
        $seq2remapi{$seq} = $now_maxalle;
        $now_maxalle++;
        push $seq2i{$seq}->@*, $i;
        push @new_alts, $seq;
    }
    for(my $i=0; $i<@new_alts; $i++) {
        $new_alts[$i]='*' if $new_alts[$i] eq '';
        $new_alts[$i] = uc $new_alts[$i];
    }
    my $new_ref = shift @new_alts;
    return(\%remapi, $new_ref, \@new_alts);
}




sub read_ref_noaug {
    my ($in) = @_;
    my %ref;
    my %nonref;
    my $seq_now;
    my $I = open_in_fh($in);
    my $chr_now;
    while(<$I>){
        chomp;
        next unless $_;
        if (/^>(\S+)(\s|$)/) {
            my $chr_raw = $1;
            if ($chr_raw=~/^\d+$/) {
                $seq_now = \$ref{$chr_raw};
            } else {
                $chr_raw=~/^([\w_]+)\[(\d+)\]$/ or die; # zzperl vg2gfa2fa
                #$chr_raw=~/^(\S+)_(\d+)_\d+(\s|$)/ or die; # gfatools gfa2fa
                my $chr = $1;
                my $pos = $2;
                $seq_now = \$nonref{$chr}{$pos};
            }
            $$seq_now = '';
            next;
        }
        $$seq_now .= $_;
    }
    my $total_length=0;
    my $total_chrs = scalar(keys %ref);
    $total_length += length($ref{$_}) foreach keys (%ref);
    my $stat = "read $total_chrs contig / total $total_length bp";
    return (\%ref, $stat, \%nonref);
}


