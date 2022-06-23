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
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw(process_alts aln_halign);



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
    $HAlignC = $main::config->{stmsa} // confess "stmsa path not defined in config file";
    $muscle3 = $main::config->{muscle3} // undef;
    $mafft = $main::config->{mafft} // undef;
    $tmp_dir = $main::tmp_dir // confess "tmp_dir not defined";
    $debug = $main::debug // confess "debug not defined";
    $align_level = $main::align_level // 1;

    foreach my $refs ([$HAlignC, 'stmsa'], 
                      [$muscle3, 'muscle3'],
                      [$mafft, 'mafft'])  {
        my ($exe, $name) = @$refs;
        &check_bin_path($exe, $name) if defined $exe;
    }

    $ALN_PARAMS_max_tryi = $main::config->{realign_max_try_times_per_method};

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
    my ($lefti, $length) = @_;
    my $next;
    if ( $length > 0 and $$lefti{HAlignC}>0 ) { ##########
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
    foreach my $ialt(0..$max_alts) {
        my $alt = $$alts[$ialt];
        say $chld_in ">$ialt";
        say $chld_in $alt;
    }
    close $chld_in;
    my $aln_alts = &read_fa_fh($chld_out);
    waitpid( $aln_pid, 0 );
    my $aln_exit_status = $? >> 8;
    return($aln_exit_status, $aln_alts);
}




sub process_alts {
    my ($alts, $max_alts, $alt_max_length) = @_;
    my %lefti;
    $lefti{$_}=$ALN_PARAMS_max_tryi foreach keys %ALN_PARAMS; # max try 3 times
    for(my $i=0; $i<scalar(@$alts); $i++) {
        # prase for Halign and mafft (not muscle)
        $$alts[$i]='-' if $$alts[$i] eq '';
    }
    REDO_process_alts:
    my ($selsoft, $tryi) = &select_aln_software(\%lefti, $alt_max_length);
    return undef unless defined $selsoft;
    my ($aln_exit_status, $aln_alts);
    if ($selsoft =~/^HAlign/) {
        ($aln_exit_status, $aln_alts) = &aln_halign($selsoft, $alts, $max_alts);
    } else {
        ($aln_exit_status, $aln_alts) = &aln_pipeline($selsoft, $alts, $max_alts);
    }
    unless (exists $$aln_alts[0]) {
        carp "Warn! no alt[0] from aln software($selsoft)! redo! No. $tryi";
        goto REDO_process_alts;
    }
    say STDERR Dumper($aln_alts) if $debug;
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
        my $seq = $$alts[$ialt];
        my $l = length($seq);
        confess("length:  $l != $aln_seqlen") if $l != $aln_seqlen;
        tie $tie_seqs[$ialt]->@*, 'Tie::CharArray', $seq;
    }
    my $miss_start = -9;
    my $is_miss=0;
    my $ref_missn=0;
    my $old_sarray=[];
    my $is_ref_miss_at_start=0;
    my $old_changearray = [];
    my ($last_nomiss_pos, $last_nomiss_sarray) = (0, []);
    for (my $i = 0; $i < $aln_seqlen; $i++) {
        my $sarray = &get_sarray(\@tie_seqs, $i, $max_alts, 1);
        my ($is_same_now, $is_miss_now, $is_ref_miss_now) = &sarray_is_same_miss($sarray, $max_alts);
        say STDERR "!!" . " $i " . ($i-$ref_missn) . ' : ' . $tie_seqs[0][$i] . " is_same_now$is_same_now, is_miss_now$is_miss_now, is_ref_miss_now$is_ref_miss_now, is_miss$is_miss" if $debug;
        my $changearray = [];
        $changearray = &sarray_diff_array($old_sarray, $sarray) if $i>0;
        if ($is_miss_now==1) {
            $ref_missn++ if $is_ref_miss_now;
            if ($is_miss==1) {
                # continue missing
                if ($align_level==2) {
                    my $diff_is_same = &numarray_is_same($old_changearray, $changearray);
                    if ( $diff_is_same==0 and length($$old_sarray[0])>0 and 
                            $is_ref_miss_now==0 ) {
                        $muts{$miss_start} = $old_sarray;
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
                $muts{$i - $ref_missn} = $sarray;
                $is_miss = 0 if $is_miss==1;
            } elsif ($is_same_now==1) {
                # do nothing;
            } else {
                confess();
            }
        }
        $old_changearray = $changearray;
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
    say STDERR "!!muts: " . Dumper \%muts if $debug;
    return (\%muts, $aln_seqlen);
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
    if ($first eq '-') {
        $is_ref_miss=1;
    }
    my %seqs;
    for(my $i = 0; $i<= $max_alts; $i++) {
        my $seq = $$sarray[$i];
        $seqs{$seq}++;
    }
    my $is_miss = 0;
    if (exists $seqs{'-'}) {
        $is_miss=1;
        delete $seqs{'-'};
    }
    my $is_same = scalar(keys %seqs)>1 ? 0 : 1 ;
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
    my ($fh) = @_;
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
    my @seqids = sort {$a<=>$b} keys %seqs;
    return( [ @seqs{@seqids} ] );
}


return 1;
