#!/usr/bin/env perl
use strict;
use warnings;
use v5.10;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/../lib";
use zzIO;
use zzBed;
use List::Util qw(sum);
no warnings 'experimental::smartmatch';
use Data::Dumper;

my ($in, $ref_fasta_file, $opt_help, $software, $not_rename, $newname);

my ($out, $out_inv, $out_ins, $out_del);
my $clear_info = 1;
my $remove_no_mut = 1;
my $remove_missing = 1;
my $filter_by_PASS = 1;
my $not_prase = 0; # only filter
my $force = 0;

my $max_len = 50000;
my $min_dp = 5;
my @known_softwares = qw/pbsv svim cutesv sniffles2 assemblytics guess_per_line svimmer svim_asm/;

my $known_softwares_str = join ', ', @known_softwares;

sub help {
    my ($info) = @_;
    say STDERR $info if $info;
    say <<EOF;
Usage: $0 -i <in.vcf> -r <ref.fa> -o <out.vcf>
    -i, --in <in.vcf>       input vcf file
    -r, --ref <ref.fa>      reference fasta file
    -o, --out <out.vcf>     output vcf file
    --out_inv <xx.bed> output inv bed file (if not defined, will not split inversions)
    --out_ins <xx.vcf> output ins vcf file (if not defined, will not split insertions)
    --out_del <xx.vcf> output del vcf file (if not defined, will not split deletions)
    -s, --software <str>    software name: $known_softwares_str
    -h, --help              print this help message
    --max_len <int>         max length of SV, default $max_len
    --min_dp <int>          min reads depth of SV. 0 to disable filter default $min_dp
    --not_rename            do not rename the sample name in vcf header
    --newname <str>         rename the sample name in vcf header
    --clear_info <int>      clear info field, default $clear_info
    --remove_no_mut <int>   remove no mutation sites(0/0), default $remove_no_mut
    --remove_missing <int>  remove missing sites(./.), default $remove_missing
    --not_prase             only filter, do not parse
    --filter_by_PASS        Filter by FILTER : PASS. default $filter_by_PASS
    --force
EOF
    exit(-1);
}

my $ARGVs = join ' ', $0, @ARGV;

GetOptions (
        'help|h!' => \$opt_help,
        'i|in=s' => \$in,
        'r|ref=s' => \$ref_fasta_file,
        'o|out=s' => \$out,
        'out_inv=s' => \$out_inv,
        'out_ins=s' => \$out_ins,
        'out_del=s' => \$out_del,
        'max_len=s' => \$max_len,
        'd|min_dp=i' => \$min_dp,
        's|software=s' => \$software,
        'not_rename!' => \$not_rename,
        'clear_info=s' => \$clear_info,
        'remove_no_mut=i' => \$remove_no_mut,
        'remove_missing=i' => \$remove_missing,
        'newname=s' => \$newname,
        'not_prase!' => \$not_prase,
        'filter_by_PASS=i' => \$filter_by_PASS,
        'force!' => \$force,
);

&help() if $opt_help;
&help("--in not defined") unless defined $in;
&help("--out not defined") unless defined $out;
#&help("--out_inv not defined") unless defined $out_inv;
unless (defined $out_inv) {
    say STDERR "Info: --out_inv not defined, will not split inversions";
}
#&help("--software not defined") unless defined $software;
&help("--ref not defined") unless defined $ref_fasta_file;

$software = lc($software) if defined $software;
&help() if(defined $software and ! $software ~~ @known_softwares);


my $ref_fasta;
say STDERR "Now reading reference fasta file: $ref_fasta_file";
$ref_fasta = &read_ref_fasta($ref_fasta_file) if defined $ref_fasta_file;
say STDERR "Done reading reference fasta file";
my $I = open_in_fh($in);
my $O = open_out_fh($out);
my ($inv_bed,$ins_fh,$del_fh);
$inv_bed = new zzBed({outfile=>$out_inv, max_len=>9999999999}) if defined $out_inv;
$ins_fh = open_out_fh($out_ins) if defined $out_ins;
$del_fh = open_out_fh($out_del) if defined $out_del;
#$inv_bed = new zzBed({outfile=>$out_inv, max_len=>$max_len});
my @header;

LINE:while(my $line = <$I>) { # vcf-header
    chomp $line;
    my $software_guessed;
    if($line=~/^#/) {
        if($line=~/^##/) {
            if($line=~/^##source=([a-zA-Z0-9]+)/) {
                $software_guessed = $1;
                $software_guessed = lc $software_guessed;
                if(!defined $software) {
                    if($line=~/SVIM-asm-v1/) {
                        $software = 'svim_asm';
                        say STDERR "software guessed: $software";
                    } elsif ($software_guessed ~~ @known_softwares) {
                        $software = $software_guessed;
                        say STDERR "software guessed: $software";
                    }else {
                        die "Fail to guess software and software name not defined, please use -s to specify it";
                    }
                }
            }
            say $O $line;
            say $ins_fh $line if defined $ins_fh;
            say $del_fh $line if defined $del_fh;
            next;
        }
        say $O qq+##CommandLine="$ARGVs"+;
        say $ins_fh qq+##CommandLine="$ARGVs"+ if defined $ins_fh;
        say $del_fh qq+##CommandLine="$ARGVs"+ if defined $del_fh;
        @header = split /\t/, $line;
        if(defined $not_rename) {
            say $O $line;
            say $ins_fh $line if defined $ins_fh;
            say $del_fh $line if defined $del_fh;
        } else {
            $header[8] //= 'GT';
            my $new_name_now = defined $newname ? $newname : $software;
            say $O join "\t", @header[0..8], $new_name_now;
            say $ins_fh join "\t", @header[0..8], $new_name_now if defined $ins_fh;
            say $del_fh join "\t", @header[0..8], $new_name_now if defined $del_fh;
        }
        last;
    }
}

LINE:while(<$I>) { # vcf
    chomp;
    my @F = split /\t/;
    my ($chr, $pos, $svid, $ref, $alts) = @F[0,1,2,3,4];
    my $infos = $F[7];
    my $gts = $F[9];
    next if $remove_no_mut==1 and $gts=~m#^0/0#; # skip non-mutation sites, especially for sniffles2
    next if $remove_missing==1 and $gts=~m#^\.\/\.#; # skip missing sites
    my $dp;
    my $type;
    my $software_ = $software;
    if($software eq 'guess_per_line') {
        $software_ = &guess_software_by_svid($svid);
        next LINE unless defined $software_;
        #say STDERR "software guessed: $software";
    }
    if($filter_by_PASS==1 and $F[6] ne 'PASS') { # FILTER
        next LINE;
    }
    if($software_ eq 'sniffles2') {
        $dp = &get_DR_DV(\@F);
        if( $F[4] eq '<DEL>' ) {
            $type = 'DEL';
            next LINE if $min_dp>0 and $dp < $min_dp;
            &Update_ref_alt(\@F, 'DEL');
        } elsif( $F[4] eq '<INS>' ) {
            $type = 'INS';
            next LINE if $min_dp>0 and $dp < $min_dp;
            next LINE;
        } elsif( $svid =~ /Sniffles2.INS\.\w+$/ ) {
            $type = 'INS';
            next LINE if $min_dp>0 and $dp < $min_dp;
            &Update_ref_alt(\@F, 'INS');
        } elsif( $F[4] eq '<INV>' ) {
            $type = 'INV';
            if(defined $out_inv) {
                $inv_bed->add_inv_bed(\@F, $svid); next LINE;
            } else {
                &Update_ref_alt(\@F, 'INV');
            }
        } elsif( $F[4] eq '<DUP>' ) {
            $type = 'INS';
            next LINE if $min_dp>0 and $dp < $min_dp;
            &Update_ref_alt(\@F, 'DUP');
            &Update_ref_alt(\@F, 'INS');
        } else {
            next LINE;
        }
    } elsif($software_ eq 'pbsv') {
        $dp = &get_DP(\@F);
        if($svid =~ /pbsv.DEL\.\w+$/) {
            $type = 'DEL';
            next LINE if $min_dp>0 and $dp < $min_dp;
            # do nothing
        } elsif($svid =~ /pbsv.INS\.\w+$/) {
            $type = 'INS';
            next LINE if $min_dp>0 and $dp < $min_dp;
        } elsif( $F[4] eq '<DUP>') {
            $type = 'INS';
            next LINE if $min_dp>0 and $dp < $min_dp;
            &Update_ref_alt(\@F, 'DUP');
        } elsif( $F[4] eq '<INV>') {
            $type = 'INV';
            if(defined $out_inv) {
                $inv_bed->add_inv_bed(\@F, $svid); next LINE;
            } else {
                &Update_ref_alt(\@F, 'INV');
            }
        } else {
            # say STDERR "unknown type: @F";
            next LINE;
        }
    } elsif($software_ eq 'svim') {
        $dp = &get_DP(\@F);
        next LINE if $F[9]=~m#^0\/0#;
        next LINE if $F[9]=~m#^\./\.#;        
        if($svid =~ /svim.DEL\.\w+$/) {
            $type = 'DEL';
            next LINE if $min_dp>0 and $dp < $min_dp;
            &Update_ref_alt(\@F, 'DEL') if $F[4] eq '<DEL>';
            # do nothing
        } elsif($svid =~ /svim.INS\.\w+$/) {
            $type = 'INS';
            next LINE if $min_dp>0 and $dp < $min_dp;
            if ($F[4] eq '<INS>') {
                my $ins_seq = &svim_ins_cal_seq(\@F);
                if (! defined $ins_seq) {
                    if(!$force) {
                        say STDERR "========================";
                        say STDERR "========================";
                        say STDERR "Error: No ins_seq: svim_ins_cal_seq error";
                        say STDERR "================== Line: $.";
                        die "@F";
                        say STDERR "========================";
                    }
                } else {
                    $F[4] = $ins_seq;
                }
                
                &Update_ref_alt(\@F, 'INS');
            }
            &Update_ref_alt(\@F, 'INS') if($F[3] eq 'N');
        } elsif( $svid =~ /svim.INV\.\w+$/ ) {
            $type = 'INV';
            if (defined $out_inv) {
                $inv_bed->add_inv_bed(\@F, $svid); next LINE;
            } else {
                &Update_ref_alt(\@F, 'INV');
            }
        } elsif( $F[4]=~m/^<DUP/) {
            $type = 'INS';
            next LINE if $min_dp>0 and $dp < $min_dp;
            &Update_ref_alt(\@F, 'DUP');
            if($F[3] eq 'N') {
                &Update_ref_alt(\@F, 'INS');
            }
        } else {
            next LINE;
        }
    } elsif($software_ eq 'assemblytics') {
        # do nothing
    } elsif($software_ eq 'svim_asm') {
        if( $F[4] =~ '^<DUP') {
            $type = 'INS';
            &Update_ref_alt(\@F, 'DUP');
            if($F[3] eq 'N') {
                &Update_ref_alt(\@F, 'INS');
            }
        } elsif($F[2]=~/^svim_asm\.INV\.\d+$/) {
            $type = 'INV';
            if (defined $out_inv) {
                $inv_bed->add_inv_bed(\@F, $svid); next LINE;
            }
        } elsif($F[2]=~/BND/) {
            next LINE;
        }
    } elsif($software_ eq 'cutesv') {
        $dp = &get_DR_DV(\@F);
        if($svid =~ /cuteSV.DEL\.\w+$/) {
            $type = 'DEL';
            next line if $min_dp>0 and $dp < $min_dp;
            &Update_ref_alt(\@F, 'DEL');
            # do nothing
        } elsif($svid =~ /cuteSV.INS\.\w+$/) {
            $type = 'INS';
            next line if $min_dp>0 and $dp < $min_dp;
        } elsif( $F[4] eq '<DUP>') {
            $type = 'INS';
            next line if $min_dp>0 and $dp < $min_dp;
            &Update_ref_alt(\@F, 'DUP');
        } elsif( $F[4] eq '<INV>' ) {
            $type = 'INV';
            if (defined $out_inv) {
                $inv_bed->add_inv_bed(\@F, $svid); next LINE;
            } else {
                &Update_ref_alt(\@F, 'INV');
            }
        } else {
            next LINE;
        }
    } elsif($software_ eq 'svimmer') {
        /SVTYPE=(\w+)/ or die;
        my $svtype = $1;
        next if $svtype eq 'BND';
        my $alts = $F[4];
        if($alts =~ /^<DUP/) {
            $type = 'INS';
            &Update_ref_alt(\@F, 'DUP');
            &Update_ref_alt(\@F, 'INS') if $F[3] eq 'N';
        } elsif ($alts =~ /^<DEL/) {
            $type = 'DEL';
            &Update_ref_alt(\@F, 'DEL');
        } elsif ($alts =~ /^<INV/) {
            $type = 'INV';
            if (defined $out_inv) {
                $inv_bed->add_inv_bed(\@F, $svid); next LINE;
            } else {
                &Update_ref_alt(\@F, 'INV');
            }
        } elsif ($alts =~ /^<INS/) {
            $type = 'INS';
            &Update_ref_alt(\@F, 'INS') if $F[3] eq 'N';
            # next LINE;
        } elsif ($alts =~ /(\[|<)/) {
            die "@F";
        } else {
            # do nothing
        }
        $F[8] = 'GT';
        $F[9] = '1/1';
    } else {
        say STDERR "Warn: skip line for unknown software: $software: @F";
        next LINE;
    }
    #die Dumper \@F;
    if( $not_prase==0 and ($F[4]=~m/</ or $F[3]=~m/</) ) {
        unless($software eq 'svimmer' and $F[4] eq '<INS>') {
            next LINE;
        }
    }
    my $svlen = &max(length($F[3]), length($F[4]));
    next if $svlen > $max_len;
    if($clear_info eq '1') {
        # remove all infos
        $F[7] = '.' 
    } elsif($clear_info eq '0') {
        # do nothing
    } else {
        # remove all infos except SVLEN ; END ; SVTYPE
        $F[7] = &clear_info_cal($F[7], $clear_info);
    }
    if(defined $ins_fh and $type eq 'INS') {
        say $ins_fh join "\t", @F;
    } elsif(defined $del_fh and $type eq 'DEL') {
        say $del_fh join "\t", @F;
    } else {
        say $O join "\t", @F;
    }
}



if (defined $out_inv) {
    #die Dumper $inv_bed;
    $inv_bed->sort_bed();
    $inv_bed->merge_bed();
    $inv_bed->print_bed();
}


exit;

sub clear_info_cal {
    my ($info, $keeps) = @_;
    my @keeps = split /,/, $keeps;
    my @infos = split /;/, $info;
    my @infos_new;
    foreach my $info (@infos) {
        my ($key, $value) = split /=/, $info;
        if($key ~~ @keeps) {
            push @infos_new, $info;
        }
    }
    return join ';', @infos_new;
}


sub svim_ins_cal_seq {
    my ($F) = @_;
    my $info = $$F[7];
    # SEQS
    $info =~ /SEQS=([^;]+)/ or return undef;
    my ($seqs) = $1;
    my @seqs = split /,/, $seqs;
    my %seqs_count;
    for my $seq (@seqs) {
        $seqs_count{$seq}++;
    }
    #die Dumper \%seqs_count;
    my @seqs_sorted = sort {$seqs_count{$b} <=> $seqs_count{$a}} keys %seqs_count;
    my $max_count_seq = $seqs_sorted[0];
    #my @seqs_sorted = grep {$seqs_count{$_} == $max_count_seq} @seqs_sorted_i;
    return $max_count_seq;
}

sub guess_software_by_svid {
    my ($svid) = @_;
    $svid = lc($svid);
    if($svid =~ /pbsv/) {
        return 'pbsv';
    } elsif($svid =~ /svim_asm/) {
        return 'svim_asm';
    } elsif($svid =~ /svim/) {
        return 'svim';
    } elsif($svid =~ /cutesv/) {
        return 'cutesv';
    } elsif($svid =~ /assemblytics/) {
        return 'assemblytics';
    } elsif($svid =~ /sniffles2/) {
        return 'sniffles2';
    } else {
        say STDERR "Warn: skip line for unknown software: $svid";
        return undef;
    }
}

sub max {
    my (@datas) = @_;
    my $max = $datas[0];
    for my $data (@datas) {
        $max = $data if $data > $max;
    }
    return $max;
}


sub Update_ref_alt {
    my ($line, $type, $svlen_parm) = @_;
    if($not_prase==1) {
        return;
    }
    my $chr = $$line[0];
    my $pos = $$line[1];
    my $infos = $$line[7];
    my $svlen;
    if($svlen_parm) {
        $svlen = $svlen_parm;
    } else {
        if($infos =~ /SVLEN=([\-\de.+]+)(;|$)/) {
            $svlen = $1;
        } else {
            # END
            if($infos =~ /END=(\d+)(;|$)/) {
                $svlen = $1 - $pos + 1;
            }
        }
        die "no svlen: @$line" if ! defined $svlen;
    }
    $svlen = abs($svlen);
    die "chr not in fasta: $chr" if ! exists $$ref_fasta{$chr};
    my $ref_seq = \$$ref_fasta{$chr} // die "no ref seq for $chr";
    if($type eq 'DEL') {
        $$line[3] = substr($$ref_seq, $pos-1, $svlen);
        $$line[4] = substr($$line[3], 0, 1);
    } elsif($type eq 'INS') {
        $$line[3] = substr($$ref_seq, $pos-1, 1);
    } elsif($type eq 'INV') {
        $$line[3] = substr($$ref_seq, $pos-1, $svlen);
        $$line[4] = &revcomp($$line[3]);
    } elsif($type eq 'DUP') {
        $$line[4] = substr($$ref_seq, $pos-1, $svlen);
    } else {
        die;
    }
}

sub revcomp {
    my ($seq) = @_;
    my $rev = reverse $seq;
    $rev =~ tr/ACGTacgt/TGCAtgca/;
    $rev =~ s/[^ACGTacgt]/N/g;
    return $rev;
}

sub get_DP { # pbsv svim
    # extract DP from vcf line
    my ($line) = @_;
    my @format = split /:/, $$line[8];
    my @values = split /:/, $$line[9];
    my %hash;
    foreach my $i (0..$#format) {
        $hash{$format[$i]} = $values[$i];
    }
    my $DP = $hash{DP} // 10; # die "no DP in @$line";
    $DP = 0 if $DP eq '.';
    return($DP);
}

sub get_DR_DV { # cuteSV sniffles2
    # extract DR and DV from vcf line
    # ##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
    # ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
    my ($line) = @_;
    my @format = split /:/, $$line[8];
    my @values = split /:/, $$line[9];
    my %hash;
    foreach my $i (0..$#format) {
        $hash{$format[$i]} = $values[$i];
    }
    my $DR = $hash{DR} // 5; # die "no DR in @$line";
    my $DV = $hash{DV} // 5; # die "no DV in @$line";
    if($DR=~/,/) {
        my @DRs = split /,/, $DR;
        $DR = sum(@DRs);
    }
    $DR = 0 if $DR eq '.';
    $DV = 0 if $DV eq '.';
    my $DP = $DR + $DV;
    return($DP);
}


sub read_ref_fasta {
    my $ref_fasta_file = shift;
    my $I = open_in_fh($ref_fasta_file);
    my %ref;
    my $id;
    while(<$I>) {
        chomp;
        if(/^>(\S+)/) {
            $id = $1;
            next;
        }
        $ref{$id} .= $_;
    }
    return \%ref;
}

