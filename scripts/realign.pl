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
#use IPC::Open2;
use MCE::Flow;
use MCE::Candy;
use Getopt::Long;
use MCE::Channel;
#use File::Temp;
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
my $skip_mut_at_same_pos = 2;
our $threads = 1; # default threads
our $tmp_dir = $tmp_dir_def;
my $first_merge = 0;
my $chr_tolerance = 0;
our $force_realign = 0;
our $not_use_merge_alle_afterall = 0;
my $mask_bed_file;
my $mask_bed;
my $skip_snp = 0;

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
    --mask_bed_file <file>         Mask bed file
    --tmpdir <dir>                 Temprory directory (default: $tmp_dir)
    --chr_tolerance                Tolerance of chr name (default: $chr_tolerance)
    -h | --help                    Print this help
    --verb <int>
    --all <bool>                   Print lines even no mutation.
    --level <BITCODE>              Level &2 will not merge lines. Level &4 will split missing alles.
    --skip_snp <bool>              Skip snp
EOF
    exit(1);
}

my $ARGVs = join(" ", $0, @ARGV);

GetOptions (
        'h|help!' => \$opt_help,
        'i|in_vcf=s' => \$infile,
        'o|out_vcf=s' => \$outfile,
        'debug!' => \$debug,
        'verb=i' => \$verb,
        't|threads=i' => \$threads,
        'r|ref_fasta_file=s' => \$ref_fasta_file,
        'aug!' => \$is_aug,
        'tmpdir=s' => \$tmp_dir,
        'E|ext_bp_max=i' => \$ext_bp_max,
        'e|ext_bp_min=i' => \$ext_bp_min,
        'all!' => \$print_all,
        'level=i' => \$align_level,
        'skip_mut_at_same_pos=i' => \$skip_mut_at_same_pos,
        'first_merge!' => \$first_merge,
        'chr_tolerance!' => \$chr_tolerance,
        'force_realign!' => \$force_realign,
        'not_use_merge_alle_afterall!' => \$not_use_merge_alle_afterall,
        'mask_bed_file=s' => \$mask_bed_file,
        'skip_snp!' => \$skip_snp,
);


&usage() if $opt_help or !$infile or !$outfile or !$ref_fasta_file;

$tmp_dir = $tmp_dir_def if $tmp_dir eq 'Default';
if (! -e $tmp_dir) {
    mkdir $tmp_dir or die "cannot create tmpdir: $tmp_dir!";
}

our $config = read_config_yaml("$Bin/../config.yaml");
realign_alts::init();
$mask_bed = read_mask_bed_file($mask_bed_file) if $mask_bed_file;

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

my $VCF_APPEND = <<"EOF";
##INFO=<ID=ZORIPOS,Number=1,Type=Integer,Description="original start position before Realign of Variant">
##INFO=<ID=ZORIEND,Number=1,Type=Integer,Description="original end position before Realign of Variant">
##INFO=<ID=ZLEN,Number=1,Type=Integer,Description="original length before Realign of Variant">
##CommandLine="$ARGVs"
EOF

my $vcf_header = &read_vcf_header($I_VCF, $O_VCF, $VCF_APPEND);
my $idi_max = scalar(@$vcf_header)-1;


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
        $ext = int($max_len * 0.3)+1;
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

sub check_is_in_mask {
    my ($chr, $pos, $len) = @_;
    #die Dumper $mask_bed;
    #say STDERR "$chr $pos $len";
    return 1 if ! defined $mask_bed_file; # not enable mask
    return 0 if ! exists $mask_bed->{$chr};
    my $end = $pos + $len - 1;
    my $arrays = $mask_bed->{$chr};
    my $array_n = scalar(@$arrays);
    my @windows_in;
    foreach my $arrayi (0..$array_n-1) {
        my $array = $arrays->[$arrayi];
        my ($mask_start, $mask_end) = @$array;
        if ($mask_end < $pos) {
            next;
        } elsif ($pos < $mask_end && $end > $mask_start) {
            #return $arrayi;
            push @windows_in, $arrayi;
        } elsif ($mask_start > $end) {
            last;
        }
    }
    if (@windows_in) {
        my $min_windowsid = min(@windows_in);
        return $min_windowsid;
    } else {
        return 0;
    }
}

sub zz_mce_producer {
    my $I = $I_VCF;
    my @lines;
    my $min_start;
    my $max_end;
    my $end_real;
    my $chr_old;
    my $iline=1;
    while(my $line = <$I>) {
        chomp $line;
        next if $line=~/^#/;
        next unless $line;
        my $snpid = $.;
        #my $snpid = "$F[0]:$F[1]"; # may not unique
        my @F = split(/\t/, $line);
        next if $F[4] eq '';
        my $chr = $F[0];
        my $pos = $F[1];
        my $ref_seq = $F[3];
        my $reflen = length($ref_seq);
        my @alts = split(/,/, $F[4]);
        $F[4] = \@alts;
        if ($skip_snp==1) {
            my $alt_len_max = max(map {length($_)} @alts);
            my $is_snp = $reflen == 1 && $alt_len_max == 1;
            next if $is_snp==1; ##### skip snp
        }
        if ($align_level & 2) { # align_level==2, not merge lines
            $end_real = $pos + $reflen - 1;
            $queue->enqueue( $iline, [[\@F], $chr, $pos, $end_real]);
            $iline++;
            next;
        }
        # next if $pos < 14999102; ######### debug
        unless (defined $chr_old) { # init or new chr
            $chr_old = $chr;
            $min_start = $pos;
            ($max_end, $end_real) = &cal_max_ext_range($ref_seq, \@alts, $pos, -1);
            redo;
        }
        if ($chr_old ne $chr or $max_end+1 < $pos) { # new chr or not overlap, treated as new chunk
            #print(STDERR "2!! $chr_old\t$min_start\t$max_end\t$pos\n");
            die unless defined $min_start;
            die unless defined $end_real;
            my $lines_count = scalar @lines;
            #say STDERR "enqueue: $chr_old:$min_start-$end_real ($lines_count lines)";
            if(defined $mask_bed_file) {
                &split_queue_by_mask_bed(\@lines, \$iline);
            } else {
                $queue->enqueue( $iline, [\@lines, $chr_old, $min_start, $end_real] ) if @lines;
            }
            @lines=();
            if ($chr_old ne $chr) {
                $chr_old = undef;
                $min_start = undef;
                $max_end = undef;
                $end_real = undef;
            } else {
                $min_start = $pos;
                $max_end = $pos;
                $end_real = $pos;
            }
            redo;
        }
        $min_start //= $pos;
        my ($max_end_now, $end_real_now) = &cal_max_ext_range($ref_seq, \@alts, $pos, $max_end);
        $end_real = $end_real_now if((!defined $end_real) or $end_real < $end_real_now);
        $max_end = $max_end_now if((!defined $max_end) or $max_end < $max_end_now);
        #say STDERR "$pos: $end_real $max_end";
        push @lines, \@F;
        $iline++;
    }
    close $I;
    if(defined $mask_bed_file) {
        &split_queue_by_mask_bed(\@lines, \$iline);
    } else {
        $queue->enqueue( $iline, [\@lines, $chr_old, $min_start, $end_real] ) if @lines;
    }
    $queue->end();
}

sub split_queue_by_mask_bed {
    my ($lines, $iline_) = @_;
    my %windowi2lines;
    foreach my $line (@$lines) {
        my @F = @$line;
        my $chr = $F[0];
        my $pos = $F[1];
        my $reflen = length($F[3]);
        my $window_id = &check_is_in_mask($chr, $pos, $reflen);
        if($window_id==0) {
            my $end_real = $pos + $reflen - 1;
            $queue->enqueue( $$iline_, [[$line], $chr, $pos, $end_real] );
            $$iline_++;
        } else {
            push $windowi2lines{$window_id}->@*, $line;
        }
    }
    #die Dumper \%windowi2lines;
    foreach my $window_id (keys %windowi2lines) {
        my $lines = $windowi2lines{$window_id};
        my $line0 = $lines->[0];
        my $chr = $$line0[0];
        my $pos = $$line0[1];
        my $end_real = $pos;
        foreach my $line (@$lines) {
            my $chr = $$line[0];
            my $pos = $$line[1];
            my $ref_seq = $$line[3];
            my $reflen = length($ref_seq);
            my $end_real_now = $pos + $reflen - 1;
            $end_real = $end_real_now if $end_real < $end_real_now;
        }
        $queue->enqueue( $$iline_, [$lines, $chr, $pos, $end_real] );
        $$iline_++;
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
    if ($force_realign==0 and($is_snp==1 or scalar(@$ref_alts)==2)) { ########## SNP #  or $alt_min_length<=1
        my $retline = &gen_lines_before_process($ids2alles, $ref_alts, $chr, $win_start);
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
    my ($muts);
    my $ref_alt_max_length = $alt_max_length > $ref_len ? $alt_max_length : $ref_len;
    if( ($ref_alt_max_length>1e3 and $max_alts>100)  or
                      ($ref_alt_max_length>1e5)  or
                      ($ref_alt_max_length>1000 and $alt_max_length>400) )  {
        ($muts) = &process_alts($ref_alts, $max_alts, $alt_max_length);
    } else {
        ($muts) = &process_alts($ref_alts, $max_alts, $alt_max_length);
    }
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
        my $ref_seq_now = \$$REF_SEQS{$chr};
        my $ref = substr($$ref_seq_now, $win_start-1, $win_end - $win_start + 1);
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
        #say STDERR "Now : $idi"; # debug
        #my %sites_status1; # -1:miss   0:ref    1:alt
        #my %sites_status2; # -1:miss   0:ref    1:alt
        my %sites_status; # -1/undef:miss   0:ref    1:alt_hetero  2:alt_hetero(changed_phase)  8:alt_homo  9:alt_already_overlaped
        my $a1 = $ref;
        my $a2 = $ref;
        my $last_start = 99999999999999;
        my $alt_svs=0;
        my @old_stat = (); # phase_status, mut_is_ref_mut_hetero, is_first
        ILINE:for(my $iline=$iline_max; $iline>=0; $iline--) { # from last to first
            my $line = $$lines[$iline];
            my $alts = $$line[4];
            my $sv_start = $$line[1];
            die if $last_start < $sv_start;
            $last_start = $sv_start;
            my $ref_len = length($$line[3]);
            my $sv_end = $sv_start + $ref_len-1;
            my $alles = $$lines[$iline][$idi];
            if (!defined $alles) { # miss
                next ILINE; # skip this position for this sample
            } elsif ($$alles[0]==0 and $$alles[1]==0) { # ref
                for(my $p=$sv_start; $p<=$sv_end; $p++) {
                    if ( (!exists $sites_status{$p}) or $sites_status{$p} == -1 ) {
                        $sites_status{$p} = 0;
                    }
                }
                next ILINE; # no mutation in this sample, skip this position for this sample
            }

            my $phase_status = -1; # 1 for 0/x; 2 for x/0; -1 for other
            my $mut_is_ref_mut_hetero = 0; 
            if($$alles[0] != $$alles[1] and ($$alles[0]==0 or $$alles[1]==0)) {
                $mut_is_ref_mut_hetero = 1;
                if($$alles[0]==0) {
                    $phase_status = 1;
                } else {
                    $phase_status = 2;
                }
            }
            if(!@old_stat) { # INIT
                @old_stat = ($phase_status, $mut_is_ref_mut_hetero, 1);
            }
            ## alt
            $alt_svs++;
            # check mutation pos confict
            my $confict_status = undef; # undef:unknown  -1:confict  
                                        # 1:no_confict_until_now(no change phase)  
                                        # 2:no_confict_until_now(changed phase)  
                                        # 0:no_confict_all
            CHECKPOS:for(my $p=$sv_start; $p<=$sv_end; $p++) {
                if (exists $sites_status{$p} and $sites_status{$p} >= 1) {
                    my $id_now = $$vcf_header[$idi];
                    say STDERR "Warn: Mutation overlaped! chr:$chr, wins:$win_start, wine:$win_end pos:$p id:$id_now" if $verb>=2;
                    if($skip_mut_at_same_pos == 1) {
                        # 只要有一个位点有问题，就跳过这个个体的这一行
                        next ILINE;
                    } elsif($skip_mut_at_same_pos == 2) {
                        # TODO: 根据杂合与否，判断自动修复/反转Allele并redo
                        # 跳过冲突的位点 (自动修复)
                        if($mut_is_ref_mut_hetero==1) {
                            if((!defined $confict_status) and exists $sites_status{$p} )  { # init
                                my $sites_status_now = $sites_status{$p};
                                unless($sites_status_now == 1 or $sites_status_now == 2) {
                                    goto FIX_ALTS;
                                }
                                $confict_status //= $sites_status_now;
                                die if $confict_status == -1; # should not happen
                            }
                            CHECKPOS2:foreach my $p2($p..$sv_end) {
                                if(exists $sites_status{$p2}) {
                                    my $sites_status_now = $sites_status{$p2};
                                    if($sites_status_now != -1 and $sites_status_now != $confict_status) {
                                        $confict_status = -1;
                                        last CHECKPOS2;
                                    }
                                } else {
                                    $confict_status = 0;
                                }
                            }
                            $confict_status=0 if $confict_status>=1;
                            if($confict_status==0) { # change phase
                                say STDERR "change phase at $p" if $verb>=2;
                                ($$alles[0], $$alles[1]) = ($$alles[1], $$alles[0]);
                                $phase_status = 3-$phase_status; # 1->2, 2->1
                                last CHECKPOS;
                            } else {
                                goto FIX_ALTS;
                            }
                        } else {
                            FIX_ALTS:
                            say STDERR "Fixing alts at $p" if $verb>=1;
                            my $ref_len_max = $p-$sv_start;
                            # deep copy
                            my $alts_ = [@$alts];
                            my $ialt = $$alles[0]-1;
                            $$alts_[ $ialt ] = substr($$alts_[ $ialt ], 0, $ref_len_max);
                            $ref_len = $ref_len_max;
                            $alts = $alts_;
                            #die $$alts[ $ialt ];
                            last CHECKPOS;
                        }
                    } else { # skip_mut_at_same_pos == 0
                        # 有冲突位点就报错退出!
                        say STDERR "$ref ;\n$a1 ; \n$a2 $id_now";
                        die "DIE! chr:$chr, wins:$win_start, wine:$win_end, pos:$p, idi:$idi, sv_start:$sv_start, sv_end:$sv_end line:$iline, ref_len:$ref_len, alles:$$alles[0], $$alles[1]";
                    }
                } else {
                    if($mut_is_ref_mut_hetero==1) {
                        $sites_status{$p} = $phase_status; # mut_is_ref_mut_hetero
                    } else {
                        $sites_status{$p} = 8; # homo with mut
                    }
                }
            }
            if($first_merge==1 and $old_stat[2]!=1) { # not first SV in window, 隔行changed_phase
                if($old_stat[1] == 1 and $mut_is_ref_mut_hetero == 1) { # mut_is_ref_mut_hetero for both this site and old site
                    if($old_stat[0] == $phase_status) { # both 0/x or both x/0
                        # change phase
                        say STDERR "change phase for $sv_start" if $verb>=2;
                        ($$alles[0], $$alles[1]) = ($$alles[1], $$alles[0]);
                        #say STDERR "$$alles[0], $$alles[1]";
                        $phase_status = 3-$phase_status; # 1->2, 2->1
                        for(my $p=$sv_start; $p<=$sv_end; $p++) {
                            $sites_status{$p} = $phase_status;
                        }
                    }
                }
            }
            @old_stat = ($phase_status, $mut_is_ref_mut_hetero, 0);
            if ($$alles[0] != 0) {
                my $alt = $$alts[ $$alles[0]-1 ];
                say STDERR "a1:($idi) $a1, $chr : $sv_start, $win_start~$win_end, $sv_start-$win_start, $ref_len, $alt" if $debug;
                substr($a1, $sv_start-$win_start, $ref_len, $alt);
            }
            if ($$alles[1] != 0) {
                my $alt = $$alts[ $$alles[1]-1 ];
                say STDERR "a2:($idi) $a2, $chr : $sv_start, $win_start~$win_end, $sv_start-$win_start, $ref_len, $alt" if $debug;
                substr($a2, $sv_start-$win_start, $ref_len, $alt);
            }
        }
        my ($count_all, $count_miss) = (0,0);
        for(my $p=$win_start; $p<=$win_end; $p++) {
            if((! exists $sites_status{$p}) or $sites_status{$p} == -1) {
                $count_miss++;
            }
        }
        $count_all = $win_end-$win_start+1;
        if ($count_all <= $count_miss) {
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
    #say STDERR Dumper \%new_ref_alts;
    #die;
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
            if ($chr_tolerance==1 or $chr_raw=~/^\d+$/) {
                $seq_now = \$ref{$chr_raw};
            } elsif ($chr_raw=~/^([\w_]+)\[(\d+)\]$/) {  # zzperl vg2gfa2fa
                #$chr_raw=~/^(\S+)_(\d+)_\d+(\s|$)/ or die; # gfatools gfa2fa
                my $chr = $1;
                my $pos = $2;
                $seq_now = \$nonref{$chr}{$pos};
            } else {
                die "Error: Chr name error! without --chr_tolerance, Chr must be numeric only: $in : $chr_raw";
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

sub read_mask_bed_file {
    my ($in) = @_;
    say STDERR "Reading mask file: $in";
    my %mask;
    my $I = open_in_fh($in);
    while(<$I>){
        chomp;
        next unless $_;
        my ($chr, $start, $end) = split /\t/;
        $mask{$chr}{$start+1} = $end;
    }
    my %mask2;
    foreach my $chr (keys %mask) {
        my @array;
        foreach my $start (sort {$a <=> $b} keys %{$mask{$chr}}) {
            my $end = $mask{$chr}{$start};
            push @array, [$start, $end];
        }
        $mask2{$chr} = \@array;
    }
    return(\%mask2);
}
