#!/usr/bin/env perl
#===============================================================================
#
#         FILE: realign_gen_mask.pl
#
#        USAGE: ./realign_gen_mask.pl  
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
#      CREATED: 3/24/2023 11:54:32 AM
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

use realign_alts qw/aln_halign process_alts_aln_only/;
use read_config qw/read_config_yaml/;
use Identity qw/seq2identity/;



my ($infile, $outfile, $opt_help);

my $min_identity = 0.5;
my $force_merge_fold = 5;
our $tmp_dir_def = '/run/user/' . `id -u`; chomp $tmp_dir_def;
our $tmp_dir = $tmp_dir_def;
our $debug = 0;
my $ext_size = 1000;
my $ext_size_print = 20;

sub usage {
    my $msg = shift;
    print STDERR "$msg\n" if $msg;
    print STDERR <<"EOF";
USAGE: perl $0 [options]
options:
    -i | --in_vcf <file>           Input vcf file
    -o | --out_vcf <file>          Output vcf file
    --min_identity <float>         Min identity for indel, default $min_identity
    --force_merge_fold <int>       Force merge fold for insertion or deletion, default $force_merge_fold
    --tmpdir <dir>                 Temprory directory (default: $tmp_dir)
    --ext_size <int>               Ext size for indel, default $ext_size
    -h | --help                    Print this help message
    --ext_size_print <int>         Ext size for print, default $ext_size_print
EOF
    exit(1);
}

my $ARGVs = join(" ", $0, @ARGV);

GetOptions (
        'h|help!' => \$opt_help,
        'i|in_vcf=s' => \$infile,
        'o|out_vcf=s' => \$outfile,
        'min_identity=f' => \$min_identity,
        'force_merge_fold=i' => \$force_merge_fold,
        'tmpdir=s' => \$tmp_dir,
        'ext_size=i' => \$ext_size,
        'ext_size_print=i' => \$ext_size_print,
);

usage('HELP: -h') if $opt_help;
usage("in_vcf not defined") unless defined $infile;
usage("out_vcf not defined") unless defined $outfile;

our $config = read_config_yaml("$Bin/../config.yaml");
realign_alts::init();

my $I = open_in_fh($infile);
my $O = open_out_fh($outfile);
my $zzw_ins = new zzwindow({ext_size=>$ext_size, outfh=>$O, min_identity=>$min_identity, type=>'ins'});
my $zzw_del = new zzwindow({ext_size=>$ext_size, outfh=>$O, min_identity=>$min_identity, type=>'del'});


while(my $line = <$I>) { # vcf
    next if $line=~/^#/;
    chomp $line;
    my @F = split /\t/, $line;
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = @F;
    my @alts = split /,/, $alt;
    $F[4] = \@alts;
    if( &cal_is_del(\@F) ){
        $zzw_del->add(\@F);
    } elsif( &cal_is_ins(\@F) ){
        $zzw_ins->add(\@F);
    } else {
        say STDERR "Unclasfied line: $line";
    }
}

if($zzw_ins->{old_window_Fs}) {
    $zzw_ins->print_window();
}

if($zzw_del->{old_window_Fs}) {
    $zzw_del->print_window();
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


{
    package zzwindow;
    use List::Util qw/max min/;
    sub new {
        my ($class,$args) = @_;
        my $self = {
            old_chr => undef,
            old_pos => undef,
            old_pos_maxcov => undef,
            ext_size => $$args{ext_size} // 20,
            type => $$args{type} // die,
            min_identity => $$args{min_identity} // die,
            old_window_Fs => [],
            outfh => $$args{outfh} // die,
            force_merge_fold => $$args{force_merge_fold} // 5,
        };
        $self->{identity} = new Identity({mcl_group_threshold_identity=>$self->{min_identity}, 
                                            seqs_maxi=>1,
                                            mcl_group_threshold_diff=>99999});
        die "type ~~ ins/del" unless $self->{type} eq 'del' or $self->{type} eq 'ins';
        bless $self, $class;
        return $self;
    }


    sub add {
        my ($self, $F) = @_;
        my ($chr, $pos) = ($$F[0], $$F[1]);
        my $reflen = length($$F[3]);
        if(!defined $self->{old_chr}) {
            # init
            $self->create_new_window($F);
        } elsif($chr ne $self->{old_chr}) {
            # new chr
            $self->print_window();
            $self->create_new_window($F);
        }
        die unless defined $self->{old_pos};
        my $dis = $pos - $self->{old_pos};
        die if $dis<0;
        my $dis_maxcov = $pos - $self->{old_pos_maxcov};
        if($dis_maxcov > $self->{ext_size} or $self->check_pos_compat($F)==0) { # new window
            $self->print_window();
            $self->create_new_window($F);
            return;
        }
        if($self->{type} eq 'del') { # ext window
            if($self->check_is_same_del($F, $self->{old_window_Fs}->[-1])) {
                $self->add_to_window($F);
            } else {
                $self->print_window();
                $self->create_new_window($F);
            }
        } elsif($self->{type} eq 'ins') { # ext window
            if($self->check_is_same_ins($F, $self->{old_window_Fs}->[-1])) {
                $self->add_to_window($F);
            } else {
                $self->print_window();
                $self->create_new_window($F);
            }
        } else {
            die;
        }
    }

    sub print_window {
        my ($self) = @_;
        my $old_window_Fs = $$self{old_window_Fs};
        my $chr = $$old_window_Fs[0][0];
        my $min_pos = $$old_window_Fs[0][1];
        my $max_pos = -1;
        my $max_pos_maxcov = -1;
        foreach my $F (@$old_window_Fs) {
            my $pos = $$F[1];
            my $reflen = length($$F[3]);
            my $pos_maxcov = $pos + $reflen - 1;
            $max_pos = $pos if $pos > $max_pos;
            $max_pos_maxcov = $pos_maxcov if $pos_maxcov > $max_pos_maxcov;
        }
        my $outfh = $$self{outfh};
        my $ext_size = $$self{ext_size};
        my $window_contain_svs = scalar(@$old_window_Fs);
        say $outfh join "\t", $chr, $min_pos-$ext_size_print-1, $max_pos_maxcov+$ext_size_print, $self->{type}, $window_contain_svs;
        $$self{old_window_Fs} = [];
    }

    sub add_to_window {
        my ($self, $F) = @_;
        my $old_window_Fs = $$self{old_window_Fs};
        my ($chr, $pos) = ($$F[0], $$F[1]);
        my $reflen = length($$F[3]);
        my $new_pos_maxcov = $pos + $reflen - 1;
        $self->{old_pos_maxcov} = $new_pos_maxcov if $new_pos_maxcov > $self->{old_pos_maxcov};
        $self->{old_pos} = $pos if $pos > $self->{old_pos};
        push $$self{old_window_Fs}->@*, $F;
        return 0;
    }

    sub create_new_window {
        my ($self, $F) = @_;
        my ($chr, $pos) = ($$F[0], $$F[1]);
        my $reflen = length($$F[3]);
        $self->{old_chr} = $$F[0];
        $self->{old_pos} = $$F[1];
        my $pos_maxcov = $pos + $reflen - 1;
        $self->{old_pos_maxcov} = $pos_maxcov;
        $self->{old_window_Fs} = [$F];
    }

    sub check_pos_compat { # for indel, check is overlaped
        my ($self, $F) = @_;
        my ($chr, $pos) = ($$F[0], $$F[1]);
        my $reflen = length($$F[3]);
        my $pos_maxcov = $pos + $reflen - 1;
        my $old_pos_maxcov = $self->{old_pos_maxcov};
        my $old_window_Fs = $$self{old_window_Fs};
        my $old_pos_min = $$old_window_Fs[0][1];
        my $ext_size = $$self{ext_size};
        if($chr ne $self->{old_chr}) {
            return 0;
        }
        if($pos < $old_pos_maxcov + $ext_size and
                $pos_maxcov > $old_pos_min - $ext_size) {
            return 1;
        } else {
            return 0;
        }
    }

    sub check_del_range_compat {
        my ($self, $F1, $F2) = @_;
        my $start1 = $$F1[1];
        my $len1 = length($$F1[3]);
        my $end1 = $start1 + $len1 - 1;
        my $start2 = $$F2[1];
        my $len2 = length($$F2[3]);
        my $end2 = $start2 + $len2 - 1;
        my $len_max = max($len1, $len2);
        # cal overlap len
        my $overlap_start = max($start1, $start2);
        my $overlap_end = min($end1, $end2);
        my $overlap_len = $overlap_end - $overlap_start + 1;
        if($overlap_len > $len_max * $self->{min_identity}) {
            return 1;
        } else {
            return 0;
        }
    }

    sub check_is_same_del {
        my ($self, $F1, $F2) = @_;
        my $seq1 = $$F1[3]; # ref
        my $seq2 = $$F2[3]; # ref
        return 0 if 0 == $self->check_len_compat_2seq($seq1, $seq2);
        return 1 if 1 == $self->check_del_range_compat($F1, $F2);
        return 0 if 0 == $self->check_is_same_2seq($seq1, $seq2);
        return 1;
    }

    sub check_is_same_ins {
        my ($self, $F1, $F2) = @_;
        my $alts1 = $$F1[4]; # alts
        my $alts2 = $$F2[4]; # alts
        my $alts1_maxlen = max(map {length($_)} @$alts1);
        my $alts2_maxlen = max(map {length($_)} @$alts2);
        foreach my $alt1 (@$alts1) {
            next if 2*length($alt1) < $alts1_maxlen;
            foreach my $alt2 (@$alts2) {
                next if 2*length($alt2) < $alts2_maxlen;
                if (1 == $self->check_len_compat_2seq($alt1, $alt2) and 
                            1 == $self->check_is_same_2seq($alt1, $alt2)) {
                    return 1;
                }
            }
        }
        return 0;
    }

    sub check_len_compat_2seq { # for del
        my ($self, $seq1, $seq2) = @_;
        my $len1 = length($seq1);
        my $len2 = length($seq2);
        my $len_min = $len1 < $len2 ? $len1 : $len2;
        my $len_diff = abs($len1 - $len2);
        if ($len_diff / $len_min > $self->{min_identity}) {
            return 0; # not compat
        } else {
            return 1; # compat
        }
    }

    sub check_is_same_2seq {
        my ($self, $seq1, $seq2) = @_;
        my $seqs = [$seq1, $seq2];
        #use Data::Dumper;
        #say STDERR Dumper $seqs;
        my ($aln_seqs) = realign_alts::process_alts_aln_only([$seq1, $seq2], -1, -1);
        #say STDERR Dumper $aln_seqs;
        my ($identity, $diff, $gap) = $self->{identity}->seq2identity($aln_seqs);
        #say STDERR "$seq1, $seq2";
        #say STDERR "identity: $identity, diff: $diff, gap: $gap";
        my $is_same = $identity >= $self->{min_identity} ? 1 : 0;
        return $is_same;
    }
}


