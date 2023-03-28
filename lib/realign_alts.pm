#
#===============================================================================
#
#         FILE: realign_alts.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zeyu Zheng (LZU), zhengzy2014@lzu.edu.cn
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 10/01/2021 10:11:37 PM
#     REVISION: ---
#===============================================================================

package realign_alts;

use strict;
use warnings;
use v5.24;

use FindBin qw($Bin);
use lib "$Bin/../lib";

#use Bio::SeqIO;
use zzIO;
use List::Util qw/max min/;
use Data::Dumper;
use Carp qw/confess carp/; # carp=warn;confess=die
use IPC::Open2;
#use Clone qw(clone);
use Storable qw(dclone);

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
            process_alts
            aln_halign
            process_alts_aln_only
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw(process_alts aln_halign process_alts_aln_only);


my $famsa;
my $HAlignC;
my $muscle3;
my $mafft;
my $tmp_dir;
my $debug;
my $align_level;
my $ALN_PARAMS_max_tryi;
my %ALN_PARAMS;
my $init = 0;

sub init {
    return if $init == 1;
    $famsa = $main::config->{famsa} // confess "famsa path not defined in config file";
    $HAlignC = $main::config->{stmsa} // confess "stmsa path not defined in config file";
    $muscle3 = $main::config->{muscle3} // undef;
    $mafft = $main::config->{mafft} // undef;
    $tmp_dir = $main::tmp_dir // confess "tmp_dir not defined";
    $debug = $main::debug // confess "debug not defined";
    $align_level = $main::align_level // 1;

    foreach my $refs ([$HAlignC, 'stmsa'], 
                      [$muscle3, 'muscle3'],
                      [$mafft, 'mafft'],
                      [$famsa, 'famsa'] )  {
        my ($exe, $name) = @$refs;
        &check_bin_path($exe, $name) if defined $exe;
    }

    $ALN_PARAMS_max_tryi = $main::config->{realign_max_try_times_per_method};
    #$ALN_PARAMS{famsa} = "$famsa -t 1 -medoidtree -gt upgma STDIN STDOUT -keep-duplicates 2>/dev/null";
    $ALN_PARAMS{famsaP} = "$famsa -t 8 -medoidtree -gt upgma STDIN STDOUT -keep-duplicates 2>/dev/null";
    $ALN_PARAMS{HAlignC} = "C";
    $ALN_PARAMS{mafft} = "$mafft --auto --thread 4 - 2>/dev/null" if $mafft;
    $ALN_PARAMS{muscle} = "$muscle3 -maxiters 16 -diags 2>/dev/null" if $muscle3;
    $init = 1;
}

sub check_bin_path {
    my ($bin, $name) = @_;
    return 1 if -x $bin;
    system("which $bin >/dev/null 2>&1") == 0 or confess "$name path not found : $bin";
    return 1;
}


sub select_aln_software {
    my ($lefti, $length, $force_aln_method) = @_;
    my $next;
    if (defined $force_aln_method and exists $$lefti{$force_aln_method} and $$lefti{$force_aln_method} > 0) {
        my $tried = $ALN_PARAMS_max_tryi-$$lefti{$force_aln_method}+1;
        $$lefti{$force_aln_method}--;
        return($force_aln_method, $tried);
    } elsif ( $length<100 and exists $$lefti{muscle} and $$lefti{muscle}>0) {
        my $tried = $ALN_PARAMS_max_tryi-$$lefti{muscle}+1;
        $$lefti{muscle}--;
        return('muscle', $tried);
    } elsif ( $$lefti{famsaP}>0 ) {
        my $tried = $ALN_PARAMS_max_tryi-$$lefti{famsaP}+1;
        $$lefti{famsaP}--;
        return('famsaP', $tried);
    } elsif ( $$lefti{HAlignC}>0 ) { ##########
        my $tried = $ALN_PARAMS_max_tryi-$$lefti{HAlignC}+1;
        $$lefti{HAlignC}--;
        return('HAlignC', $tried);
    } elsif ( $length>50 and exists $$lefti{mafft} and $$lefti{mafft}>0 ) {
        my $tried = $ALN_PARAMS_max_tryi-$$lefti{mafft}+1;
        $$lefti{mafft}--;
        return('mafft', $tried);
    } elsif ( exists $$lefti{muscle} and $$lefti{muscle}>0 ) {
        my $tried = $ALN_PARAMS_max_tryi-$$lefti{muscle}+1;
        $$lefti{muscle}--;
        return('muscle', $tried);
    } elsif ( exists $$lefti{mafft} and $$lefti{mafft}>0 ) {
        my $tried = $ALN_PARAMS_max_tryi-$$lefti{mafft}+1;
        $$lefti{mafft}--;
        return('mafft', $tried);
    } elsif ( $$lefti{HAlignC}>0 ) {
        my $tried = $ALN_PARAMS_max_tryi-$$lefti{HAlignC}+1;
        $$lefti{HAlignC}--;
        return('HAlignC', $tried);
    } else {
        return(undef)
    }
}



sub aln_halign {
    my ($selsoft, $alts, $max_alts) = @_;
    my $algorithm = $ALN_PARAMS{$selsoft} // confess();
    my $fastafh = File::Temp->new(DIR=>"$tmp_dir");
    my $fasta = $fastafh->filename;
    my $fastaaln = "$fasta.aln";
    foreach my $ialt(0..$max_alts) {
        my $alt = $$alts[$ialt];
        say $fastafh ">$ialt";
        say $fastafh $alt;
    }
    close $fastafh;
    my $cmd;
    if ($algorithm eq 'C') {
        $cmd = "$HAlignC --in $fasta --out $fastaaln >/dev/null 2>&1";
    } else {
        confess();
    }
    my $aln_exit_status = system($cmd);
    undef $fastafh;
    unless (-e $fastaaln) {
        #`cp $fasta $fasta.err`;
        #say STDERR "??  $fasta.err";
        return(-1, undef);
    }
    my $alnfh = open_in_fh($fastaaln);
    my $aln_alts = &read_fa_fh($alnfh);
    unlink $fastaaln if $debug==0;
    return($aln_exit_status, $aln_alts);
}

sub aln_pipeline {
    my ($selsoft, $alts, $max_alts) = @_;
    my $aln_param = $ALN_PARAMS{$selsoft} // confess();
    my $aln_pid = open2(my $chld_out, my $chld_in, "$aln_param");
    #die "$aln_pid: $chld_in $chld_out";
    #die $aln_param;
    my @empty_ialt;
    foreach my $ialt(0..$max_alts) {
        my $alt_seq = $$alts[$ialt];
        if($alt_seq eq "" or $alt_seq eq '-' or $alt_seq eq '*') {
            push @empty_ialt, $ialt;
            next;
        }
        #die $aln_pid;
        say $chld_in ">$ialt";
        say $chld_in $alt_seq;
    }
    close $chld_in;
    my $aln_alts = &read_fa_fh($chld_out, \@empty_ialt);
    waitpid( $aln_pid, 0 );
    my $aln_exit_status = $? >> 8;
    return($aln_exit_status, $aln_alts);
}

sub cal_alts_max_len {
    my ($alts) = @_;
    my $max_len = 0;
    foreach my $alt(@$alts) {
        my $len = length($alt);
        $max_len = $len if $len>$max_len;
    }
    return($max_len);
}

sub process_alts_aln_only {
    my ($alts, $max_alts, $alt_max_length, $force_aln_method) = @_;
    my $t=$alt_max_length; #######
    undef $alt_max_length;
    if($max_alts == -1) {
        $max_alts = scalar(@$alts)-1;
    } elsif($max_alts==1) { # biallic, to skip calculation
        return({0=>$alts}, length($$alts[0]));
    } elsif (!defined $max_alts) {
        $max_alts = scalar(@$alts)-1;
        return({0=>$alts}, length($$alts[0])) if $max_alts==1;
    }
    if(! defined $alt_max_length) {
        $alt_max_length = &cal_alts_max_len($alts);
    }
    my %lefti;
    $lefti{$_}=$ALN_PARAMS_max_tryi foreach keys %ALN_PARAMS; # max try 3 times
    for(my $i=0; $i<scalar(@$alts); $i++) {
        # prase for Halign and mafft (not muscle)
        $$alts[$i]='-' if $$alts[$i] eq '';
    }
    REDO_process_alts:
    my ($selsoft, $tryi) = &select_aln_software(\%lefti, $alt_max_length, $force_aln_method);
    #say STDERR $selsoft;
    return undef unless defined $selsoft;
    my ($aln_exit_status, $aln_alts);
    if ($selsoft =~/^HAlign/) {
        ($aln_exit_status, $aln_alts) = &aln_halign($selsoft, $alts, $max_alts);
    } else {
        ($aln_exit_status, $aln_alts) = &aln_pipeline($selsoft, $alts, $max_alts);
    }
    if( $aln_exit_status!=0 ) {
        carp "Warn! aln software($selsoft) exit status($aln_exit_status)! redo! No. $tryi";
        goto REDO_process_alts;
    }
    unless (defined $aln_alts and exists $$aln_alts[0]) {
        carp "Warn! no alt[0] from aln software($selsoft)! redo! No. $tryi";
        goto REDO_process_alts;
    }
    return($aln_alts);
}


sub process_alts {
    my ($alts, $max_alts, $alt_max_length, $force_aln_method) = @_;
    if($max_alts==1) { # biallic, to skip calculation
        my %ret = (0=>$alts);
        return(\%ret, length($$alts[0]));
    }
    say STDERR "Ori no align alts: " . Dumper $alts if $debug;
    my $aln_alts = &process_alts_aln_only($alts, $max_alts, $alt_max_length, $force_aln_method);
    say STDERR "alignd alts: " . Dumper($aln_alts) if $debug;
    my ($muts, $aln_seqlen) = &alt_alts_to_muts($aln_alts, $max_alts);
    return($muts, $aln_seqlen);
}



sub alt_alts_to_muts {
    # $mut{ipos}[ialt] = seq
    my ($alts, $max_alts) = @_; # aln_alts
    #say STDERR Dumper $alts;
    my @tie_seqs;
    my %muts;
    my $aln_seqlen = length( $$alts[0] );
    foreach my $ialt (0..$max_alts) {
        my $seq = $$alts[$ialt] // die $ialt;
        my $l = length($seq);
        confess("length:  $l != $aln_seqlen") if $l != $aln_seqlen;
        tie $tie_seqs[$ialt]->@*, 'Tie::CharArray', $seq;
    }
    my $miss_start = -9;
    my $is_miss=0;
    my $ref_missn=0;
    my $old_sarray=[];
    my $is_ref_miss_at_start=0;
    my ($last_nomiss_pos, $last_nomiss_sarray) = (0, []);
    for (my $i = 0; $i < $aln_seqlen; $i++) {
        #if(exists $$old_sarray[0] and $$old_sarray[0] eq 'CGGT') {die;}
        #say STDERR $$old_sarray[0] if exists $$old_sarray[0];
        my $sarray = &get_sarray(\@tie_seqs, $i, $max_alts, 1);
        my ($is_same_now, $is_miss_now, $is_ref_miss_now) = &sarray_is_same_miss($sarray, $max_alts);
        say STDERR "!!" . " $i " . ($i-$ref_missn) . ' : ' . $tie_seqs[0][$i] . " is_same_now$is_same_now, is_miss_now$is_miss_now, is_ref_miss_now$is_ref_miss_now, is_miss$is_miss" if $debug;
        if ($is_miss_now==1) {
            $ref_missn++ if $is_ref_miss_now;
            if ($is_miss==1) {
                # continue missing
                if ($align_level & 4) {
                    my $diff_is_same = &cal_sarray_is_compatible($old_sarray, $sarray);
                    if ( $diff_is_same==0 and length($$old_sarray[0])>0 and 
                            $is_ref_miss_now==0 ) {
                        # not compatible, end missing. then start new missing
                        $muts{$miss_start} = $old_sarray;
                        #say STDERR Dumper \%muts;
                        $miss_start = $i - $ref_missn;
                        $old_sarray = [];
                        $is_miss = 0;
                        $last_nomiss_pos = $i - $ref_missn;
                        $last_nomiss_sarray = $sarray;
                    }
                }
                &append_sarray($old_sarray, $sarray);
            } elsif ($is_miss==0) {
                # start missing
                $is_miss=1;
                if ($i==0) {
                    # ref start with missing
                    $miss_start = 0;
                    $old_sarray = &get_sarray(\@tie_seqs, 0, $max_alts, 0);
                    $is_ref_miss_at_start=1 if $is_ref_miss_now==1;
                } else {
                    # start with missing in seq
                    $miss_start = $i - $ref_missn;
                    if (exists $muts{$miss_start} or
                        (defined $last_nomiss_pos and $last_nomiss_pos==$i-$ref_missn) ) {
                        # already have a non-missing site at this pos
                        # cover and replase this site
                        delete $muts{$miss_start};
                        say STDERR "delete mut $miss_start" if $debug;
                        $old_sarray = &get_sarray(\@tie_seqs, $i-1, $max_alts, 0);
                        &append_sarray($old_sarray, $sarray);
                    } else {
                        $old_sarray = &get_sarray(\@tie_seqs, $i, $max_alts, 0);
                    }
                }
            } else {confess();}
        } elsif ($is_miss==1) { # with $is_miss_now==0
            # end missing
            ($last_nomiss_pos, $last_nomiss_sarray) = (undef, undef);
            $is_miss = 0;
            my $ref_len = length($$old_sarray[0]);
            if ($ref_len>0) {
                $muts{$miss_start} = $old_sarray;
                $old_sarray = []; $miss_start = -9;
                redo;
            } else {
                &append_sarray($old_sarray, $sarray);
                $is_ref_miss_at_start = 0;
                $muts{0} = $old_sarray;
                $old_sarray = []; $miss_start = -9;
            }
        } else {
            # not missing
            $last_nomiss_pos = $i - $ref_missn;
            $last_nomiss_sarray = $sarray;
            if ($is_same_now==0) { # with $is_miss_now==0 and $is_miss==0 and $is_ref_miss_at_start==0
                # not same
                my $now_i = $i - $ref_missn;
                $muts{$now_i} = $sarray;
                $is_miss = 0 if $is_miss==1;
            } elsif ($is_same_now==1) {
                # do nothing;
            } else {
                confess();
            }
        }
    }
    if (@$old_sarray and $miss_start>=0) {
        if ( length($$old_sarray[0])==0 ) {
            # 最后一个SV的ref是空的，将序列补到倒数第二个上
            my $last_pos_start = max(keys %muts);
            if (! defined $last_pos_start) {
                ($last_nomiss_pos, $last_nomiss_sarray) = (undef, undef);
                if (defined $last_nomiss_pos) {
                    &append_sarray($last_nomiss_sarray , $old_sarray);
                    $muts{$last_nomiss_pos} = $last_nomiss_sarray;
                } else {
                    confess "last_pos_start is undef @$old_sarray";
                }
            }
            &append_sarray( $muts{$last_pos_start} , $old_sarray);
        } else {
            $muts{$miss_start} = $old_sarray;
        }
    }
    # say STDERR "!!muts: " . Dumper \%muts;
    &merge_alle_afterall(\%muts, $$alts[0]); # 合并可以合并的兼容的位点
    #die Dumper $alts;
    return (\%muts, $aln_seqlen);
}


sub merge_alle_afterall {
    my ($muts, $refseq) = @_;
    $refseq=~s/-//g;
    my %i2groups;
    my $igroups = 0;
    my @is = sort {$a<=>$b}keys %$muts;
    #say STDERR Dumper $muts;
    return if scalar(@is)<=1;
    my %spectrums;
    my @groups2i;
    my $old_spectrum = &cal_sarray_mut_spectrum($$muts{$is[0]});
    push $groups2i[0]->@*, $is[0];
    foreach my $i (1..$#is) {
        my $pos = $is[$i];
        my $spectrums = &cal_sarray_mut_spectrum($$muts{$pos});
        if($spectrums eq $old_spectrum) {
            push $groups2i[$#groups2i]->@*, $pos;
        } else {
            push $groups2i[1+$#groups2i]->@*, $pos;
            $old_spectrum = $spectrums;
        }
    }
    #say STDERR "groups2i: " . Dumper \@groups2i;
    #say STDERR "Before, muts: " . Dumper $muts;;
    ### start to merge
    foreach my $g (reverse 0..$#groups2i) {
        my @is_now = $groups2i[$g]->@*;
        next if scalar(@is_now)==1;
        &merge_alle_afterall_update($muts, \@is_now, $refseq);
    }
    #say STDERR "After, muts: " . Dumper $muts;;
    #die Dumper $muts;
}

sub merge_alle_afterall_update {
    my ($muts, $i_to_merge, $refseq) = @_;
    my $min_i = min(@$i_to_merge);
    my $max_i = max(@$i_to_merge);
    my $nalts = scalar($$muts{$min_i}->@*);
    my $max_ilen = length($$muts{$max_i}[0]);
    my $len = $max_i + $max_ilen - $min_i;
    my $refseq_now = substr($refseq, $min_i, $len);
    my @mut;
    push @mut, $refseq_now for(0..$nalts-1);
    #say STDERR Dumper $muts;
    #say STDERR "REF: $refseq_now";
    foreach my $i (sort {$b<=>$a} @$i_to_merge) {
        my $mut_reflen = length($$muts{$i}[0]);
        for my $allei (0..$nalts-1) {
            #say STDERR join(" - ",$mut[$allei], $i, $min_i, $len, $$muts{$i}[$allei]);
            substr($mut[$allei], $i-$min_i, $mut_reflen, $$muts{$i}[$allei]);
        }
        delete $$muts{$i};
    }
    $$muts{$min_i} = \@mut;
}


=error not use
sub alt_alts_to_muts_notuse {
    # $mut{ipos}[ialt] = seq
    my ($alts, $max_alts) = @_; # aln_alts
    #say STDERR Dumper $alts;
    my @tie_seqs;
    my %muts;
    my $aln_seqlen = length( $$alts[0] );
    foreach my $ialt (0..$max_alts) {
        my $seq = $$alts[$ialt];
        my $l = length($seq);
        confess("length:  $l != $aln_seqlen") if $l != $aln_seqlen;
        tie $tie_seqs[$ialt]->@*, 'Tie::CharArray', $seq;
    }
    my $old_sarray = [];
    my $mut_in_ref_pos = 0;
    my $now_pos_in_ref = -1;
    my $toadd_sarray;
    my $toadd_posinref;
    my $toadd_sarray_ext;
    for (my $i = 0; $i < $aln_seqlen; $i++) {
        my $sarray = &get_sarray(\@tie_seqs, $i, $max_alts, 0);
        #say STDERR Dumper $sarray;
        #say STDERR Dumper $old_sarray;
        #say STDERR '--';
        my ($is_same_now, $is_miss_now, $is_ref_miss_now) = &sarray_is_same_miss($sarray, $max_alts);
        $now_pos_in_ref++ if $is_ref_miss_now==0; # now_pos_in_ref 不对 23x?
        #say STDERR "!!" . " $i " . ' : ' . $tie_seqs[0][$i] . " is_same_now$is_same_now, is_miss_now$is_miss_now, is_ref_miss_now$is_ref_miss_now, now_pos_in_ref$now_pos_in_ref" ;
        if($is_same_now==1) { # this site has no mut
            if(scalar(@$old_sarray)>0) { # has mut before
                if($$old_sarray[0] eq '') { # ref is start with missing, must extend
                    &append_sarray($old_sarray, $sarray);
                } else { # to end mut
                    if(defined $toadd_posinref) { # continue extend 战未来
                        &append_sarray($toadd_sarray_ext, $sarray);
                    } else { # init 战未来
                        #die Dumper $old_sarray;
                        $toadd_sarray = dclone($old_sarray);
                        $toadd_sarray_ext = dclone($old_sarray);
                        $toadd_posinref = $mut_in_ref_pos;
                        &append_sarray($toadd_sarray_ext, $sarray);
                    }
                    $old_sarray = [];
                    $mut_in_ref_pos = -9;
                }
                say STDERR Dumper $sarray;
                say STDERR "2>; " . Dumper $toadd_sarray_ext;
            }
        } else { # now in mut
            if(scalar(@$old_sarray)>0) { # has mut before
                my $compat = &cal_sarray_is_compatible($old_sarray, $sarray);
                my $compat_toadd = &cal_sarray_is_compatible($toadd_sarray_ext, $sarray);
                say STDERR "Now: $i " . Dumper $toadd_sarray_ext;
                if ($compat==1) { # compat, to extend
                    &append_sarray($old_sarray, $sarray);
                    &append_sarray($toadd_sarray_ext, $sarray);
                } else { # not compat
                    if($is_ref_miss_now==1) { # ref miss, must extend
                        #say STDERR "! $is_ref_miss_now " . Dumper $old_sarray;
                        &append_sarray($old_sarray, $sarray);
                        &append_sarray($toadd_sarray_ext, $sarray);
                    } else { # refseq not miss now, to end mut
                        if($$old_sarray[0] eq '') { # refseq miss since start
                            # merge with this mut
                            &append_sarray($old_sarray, $sarray);
                            &append_sarray($toadd_sarray_ext, $sarray);
                        } else { # to end old mut, start new mut
                            if(defined $toadd_posinref and ! @$old_sarray) {
                                $muts{$toadd_posinref} = $toadd_sarray;
                                say STDERR "290: toadd_sarray: " . Dumper $toadd_sarray;
                                $toadd_sarray = undef;
                                $toadd_posinref = undef;
                            } elsif(defined $toadd_posinref and @$old_sarray) {
                                # combine two mut
                                say STDERR "??". Dumper $toadd_sarray_ext;
                                $muts{$toadd_posinref} = $toadd_sarray_ext;
                                say STDERR "295: toadd_sarray: " . Dumper $toadd_sarray;
                                $toadd_sarray = undef;
                                $toadd_posinref = undef;
                            }
                            $old_sarray = $sarray;
                            $mut_in_ref_pos = $now_pos_in_ref;
                        }
                    }
                }
            } elsif($i==0) { # mut in first pos
                $old_sarray = $sarray;
                $mut_in_ref_pos = 0;
            } else { # no mut before
                $old_sarray = $sarray;
                &append_sarray($toadd_sarray_ext, $sarray);
                $mut_in_ref_pos = $now_pos_in_ref;
            }
        }
    }

    if(defined $toadd_posinref) {
        #die $toadd_posinref;
        say STDERR "!!1" .  Dumper $toadd_sarray;
        $muts{$toadd_posinref} = $toadd_sarray;
        $toadd_sarray = undef;
        $toadd_posinref = undef;
    }
    if (@$old_sarray) {
        say STDERR "!!2" .  Dumper $old_sarray;
        if ( length($$old_sarray[0])==0 ) {
            # 最后一个SV的ref是空的，将序列补到倒数第二个上
            my $last_pos_start = max(keys %muts);
            if (! defined $last_pos_start) {
                confess "last_pos_start is undef @$old_sarray";
            }
            &append_sarray( $muts{$last_pos_start} , $old_sarray);
        } else {
            $muts{$now_pos_in_ref} = $old_sarray;
        }
    }
    say STDERR "!!muts: " . Dumper \%muts if $debug;
    say STDERR Dumper \%muts;
    #die;
    return (\%muts, $aln_seqlen);
}
=cut

sub cal_sarray_mut_spectrum {
    my ($sarray) = @_;
    my %alle2i;
    my $now_allei = 0;
    my $spectrum = '';
    foreach my $alle (@$sarray) {
        if(!exists $alle2i{$alle}) { # new alle
            $alle2i{$alle} = $now_allei;
            $spectrum .= $now_allei;
            $now_allei++;
        } else { # exists alle
            my $i = $alle2i{$alle};
            $spectrum .= $i;
        }
    }
    return($spectrum);
}

sub cal_sarray_is_compatible {
    my ($sarray1, $sarray2) = @_;
    my $spectrum1 = &cal_sarray_mut_spectrum($sarray1);
    my $spectrum2 = &cal_sarray_mut_spectrum($sarray2);
    if ($spectrum1 eq $spectrum2) {
        return 1;
    } else {
        return 0;
    }
}


sub append_sarray {
    my ($sarray, $append, $max_alts) = @_;
    unless (defined $max_alts) {
        $max_alts = scalar(@$append) - 1;
    }
    if (! @$sarray) {
        for (my $i=0; $i<=$max_alts ; $i++) {
            my $app_seq = $$append[$i];
            $app_seq='' if $app_seq eq '-'; # skip miss
            $$sarray[$i] = $app_seq;
        }
        return;
    } else {
        for (my $i=0; $i<=$max_alts ; $i++) {
            my $app_seq = $$append[$i];
            next if $app_seq eq '-'; # skip miss
            $$sarray[$i] .= $app_seq;
        }
        return;
    }
}


sub get_sarray {
    # sarray{ialt} = seq
    my ($tie_seqs, $i, $max_alts, $no_sub_miss) = @_;
    $no_sub_miss //= 0;
    my @sarray;
    foreach my $ialt (0..$max_alts) {
        my $seq = $$tie_seqs[$ialt][$i];
        if ($no_sub_miss==0) {
            $seq='' if $seq eq '-';
        }
        @sarray[$ialt] = $seq;
    }
    return \@sarray;
}


sub sarray_is_same_miss {
    my ($sarray, $max_alts) = @_;
    unless (defined $max_alts) {
        $max_alts = scalar(@$sarray) - 1;
    }
    my $first = $$sarray[0];

    my $is_ref_miss=0;
    if ($first eq '-' or $first eq '') {
        $is_ref_miss=1;
    }
    my %seqs;
    for(my $i = 0; $i<= $max_alts; $i++) {
        my $seq = $$sarray[$i];
        $seqs{$seq}++;
    }
    my $is_miss = 0;
    if (exists $seqs{'-'} or exists $seqs{''}) {
        $is_miss=1;
        #delete $seqs{'-'};
    }
    my $is_same = scalar(keys %seqs)>1 ? 0 : 1 ;
    #die Dumper \%seqs;
    return ($is_same, $is_miss, $is_ref_miss);
}



sub sarray_diff_array {
    my ($sarray1, $sarray2) = @_; # sarray1 is old_array may contain multiple chrs
    my @diff;
    for (my $i = 0; $i < @$sarray1; $i++) {
        my $s1 = $$sarray1[$i];
        my $s2 = $$sarray2[$i];
        my $s1_last_word = substr($s1, -1);
        my $s2_last_word = substr($s2, -1);
        push @diff, $s1_last_word ne $s2_last_word;
    }
    return \@diff;
}

sub numarray_is_same {
    my ($numarray1, $numarray2) = @_;
    if (! @$numarray1) {return 0;}
    for (my $i=0; $i<@$numarray1; $i++) {
        return 0 if $$numarray1[$i] ne $$numarray2[$i];
    }
    return 1;
}




sub read_fa_fh {
    # seqids MUST be integers
    my ($fh, $empty_ialts) = @_;
    my %seqs;
    my $seqid_now='NA';
    while (<$fh>) {
        chomp;
        #say STDERR $_;
        if (/^>(\S+)/) {
            $seqid_now = $1;
            next;
        }
        $seqs{$seqid_now} .= $_;
    }
    close $fh;
    my $i_not_empty = 0;
    my $len;
    while(!defined $len) {
        next unless exists $seqs{$i_not_empty};
        $len = length($seqs{$i_not_empty});
        last;
    }
    die "??? all empty?" if !defined $len;
    #my $len = length($seqs{0});
    foreach my $ialt (@$empty_ialts) {
        $seqs{$ialt} = '-' x $len;
    }
    my @seqids = sort {$a<=>$b} keys %seqs;
    # check len must same
    foreach my $seqid (@seqids) {
        my $len_now = length($seqs{$seqid});
        if( $len_now != $len ) {
            #die "seqid $seqid len $len not same as 0 len $len";
            return undef;
        }
    }
    return( [ @seqs{@seqids} ] );
}


return 1;
