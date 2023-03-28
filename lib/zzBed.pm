#!/usr/bin/perl
package zzBed;

require strict;
require warnings;
use v5.10;
use Data::Dumper;

use FindBin qw($Bin);
use lib "$Bin/../lib";
use zzIO;

no warnings qw( experimental::smartmatch );

sub new {
    my $class = shift;
    my $args = shift;
    my ($file, $max_len);
    $file = $$args{outfile} if defined $$args{outfile};
    #$max_len = defined $$args{max_len} ? $$args{max_len} : 50000;
    $max_len = $$args{max_len} if defined $$args{max_len};
    my $supp_min = defined $$args{supp_min} ? $$args{supp_min} : 2;
    my $fh;
    $fh = zzIO::open_out_fh($file) if defined $file;
    my $self = {
        array => [],
        fh => \$fh,
        max_len => $max_len,
        supp_min => $supp_min,
    };
    bless $self, $class;
    return $self;
}

sub read_add_bed {
    my ($self, $file) = @_;
    my $bed_t = new zzBed({max_len=>$self->{max_len}});
    my $fh = zzIO::open_in_fh($file);
    while(<$fh>) {
        chomp;
        my @line = split /\t/;
        $bed_t->add_inv_bed2(\@line);
    }
    close $fh;
    $bed_t->sort_bed();
    $bed_t->merge_bed();
    push $$self{array}->@*, $bed_t->{array}->@*;
}

sub add_inv_bed2 {
    my ($self, $line) = @_;
    my $svlen = $$line[2] - $$line[1];
    my $max_len = $self->{max_len};
    #return undef if $svlen > $max_len;
    if(defined $max_len and $svlen > $max_len) {
        return undef;
    }
    push $$self{array}->@*, $line;
}

sub add_inv_bed {
    my ($self, $line, $info) = @_;
    my $chr = $$line[0];
    my $pos = $$line[1];
    my $infos = $$line[7];
    my $gt_info = $$line[9];
    $gt_info=~m#^(\d)[/|](\d)# or return undef;
    return undef if($1==0 and $2==0); # no mut
    my $gt = join("/", sort {$a <=> $b} ($1, $2));
    my $svlen = -1;
    $infos =~ /SVLEN=(-?\d+)(;|$)/ and $svlen = abs($1);
    my $max_len = $self->{max_len};
    return undef if $svlen > $max_len;
    my $svend = $pos + $svlen;
    my @prints = ($chr, $pos-1, $svend, $gt);
    push @prints, $info if $info;
    push $$self{array}->@*, \@prints;
}

sub sort_bed {
    my ($self) = @_;
    my $array = $$self{array};
    my $maxi = scalar(@$array) - 1;
    say STDERR "now sort bed, length: $maxi";
    my @new_is = sort {$$array[$a][0] cmp $$array[$b][0] || $$array[$a][1] <=> $$array[$b][1]} 0..$maxi;
    my @new_array = $array->@[@new_is];
    $$self{array} = \@new_array;
}

sub cal_cov_count {
    # 统计每个区间的覆盖次数
    my ($self) = @_;
    my $array = $$self{array};
    my $maxi = scalar(@$array) - 1;
    say STDERR "now cal_cov_count, length: $maxi";
    my %start_ends = ();
    for (my $i=0; $i<=$maxi; $i++) {
        my ($chr, $poss, $pose) = $$array[$i]->@*;
        push $start_ends{$chr}->@*, $poss;
        push $start_ends{$chr}->@*, $pose;
    }
    # sort
    foreach my $chr (keys %start_ends) {
        my @new_is = sort {$a <=> $b} $start_ends{$chr}->@*;
        # unique
        my @new_is2 = ();
        my $old = -1;
        foreach my $new_is (@new_is) {
            if($new_is != $old) {
                push @new_is2, $new_is;
                $old = $new_is;
            }
        }
        $start_ends{$chr} = \@new_is2;
    }
    #die Dumper \%start_ends;
    # count
    my %cov_count;
    for(my $i=0; $i<=$maxi; $i++) {
        my ($chr, $poss, $pose, $gt, $info) = $$array[$i]->@*;
        my $start_ends = $start_ends{$chr};
        my $pos_old = $poss;
        foreach my $start_ends_now ($start_ends->@*) {
            if($start_ends_now <= $pose and $start_ends_now > $poss) {
                $cov_count{$chr}{$pos_old} //= [$start_ends_now, 0, []];
                $cov_count{$chr}{$pos_old}[1]++;
                push $cov_count{$chr}{$pos_old}[2]->@*, $gt;
                #say STDERR "$i: $chr:$poss-$pose.$info. $pos_old $start_ends_now " . join"-", $cov_count{$chr}{$pos_old}->@*;
            }
            $pos_old = $start_ends_now;
        }
    }
    $$self{cov_count} = \%cov_count;
    #die Dumper $array;
    #die Dumper \%cov_count;
    return(\%cov_count);
}

sub cov_count_filter {
    my ($self, $min_count) = @_;
    my $cov_count = $$self{cov_count};
    my $array = $$self{array};
    my $maxi = scalar(@$array) - 1;
    say STDERR "now cov_count_filter, length: $maxi";
    my %cov_count_flt;
    foreach my $chr (keys %$cov_count) {
        foreach my $pos (keys %{$$cov_count{$chr}}) {
            my ($end, $count, $gts) = $$cov_count{$chr}{$pos}->@*;
            my $supp_min = $self->{supp_min};
            if($count < $supp_min) {
                next;
            }
            if($count >= $min_count) {
                my $gt = &cal_most_gt($gts);
                $cov_count_flt{$chr}{$pos} = [$end, $count, $gt];
            }
        }
    }
    $$self{cov_count_flt} = \%cov_count_flt;
    return(\%cov_count_flt);
}

sub cov_count_filter_merge {
    # 合并相邻区间
    my ($self) = @_;
    my $cov_count_flt = $$self{cov_count_flt};
    say STDERR "now cov_count_filter_merge";
    my $first = 0;
    my ($chr_old, $pos_old, $info_old);
    my %merged;
    foreach my $chr(sort keys %$cov_count_flt) {
        POS:foreach my $pos(sort {$a <=> $b} keys %{$$cov_count_flt{$chr}}) {
            if($first==0) {
                $first=1;
                $chr_old = $chr;
                $pos_old = $pos;
                $info_old = $$cov_count_flt{$chr}{$pos};
                next POS;
            }
            my ($end, $count, $gt) = $$cov_count_flt{$chr}{$pos}->@*;
                #say STDERR "$info_old->[0] $pos";
            if($chr_old eq $chr and ($info_old->[0] == $pos)) {
                # to merge
                my $len_old = $info_old->[0] - $pos_old;
                my $len = $end - $pos;
                my $gt_new = $len_old > $len ? $info_old->[2] : $gt;
                $info_old = [$end, $count, $gt_new];
                $merged{$chr_old}{$pos_old} = $info_old;
            } else {
                # not to merge
                $merged{$chr_old}{$pos_old} = $info_old;
                $chr_old = $chr;
                $pos_old = $pos;
                $info_old = $$cov_count_flt{$chr}{$pos};
            }
        }
    }
    $$self{cov_count_flt_merge} = \%merged;
    return(\%merged);
}

sub cal_most_gt {
    my ($gts) = @_;
    my $hetero_count = 0;
    my $all_count = 0;
    foreach my $gt ($gts->@*) {
        $all_count++;
        $gt=~m#^(\d)[/|](\d)$# or die "GT error: $gt";
        $hetero_count++ if $1 != $2;
    }
    die if $all_count==0;
    if($hetero_count/$all_count > 0.5) {
        return "0/1";
    } else {
        return "1/1";
    }
}

sub merge_bed {
    my ($self) = @_;
    my $array = $$self{array};
    my $maxi = scalar(@$array) - 1;
    say STDERR "now merge_bed, length: $maxi";
    my ($chr_old, $poss_old, $pose_old) = $$array[0]->@*;
    for (my $i=1; $i<=$maxi; $i++){
        my ($chr, $poss, $pose) = $$array[$i]->@*;
        if($chr eq $chr_old and $poss < $pose_old) { # overlap
            $$array[$i-1][2] = $pose;
            $$array[$i-1][4] = 'merged' if exists $$array[$i-1][4];
            splice(@$array, $i, 1);
            $maxi--;
            redo;
        }
        $chr_old = $chr;
        $poss_old = $poss;
        $pose_old = $pose;
    }
    say STDERR "after merge_bed, length: $maxi";
}

sub print_bed {
    my ($self) = @_;
    my $array = $$self{array};
    my $maxi = scalar(@$array) - 1;
    for (my $i=0; $i<=$maxi; $i++) {
        my $line = join("\t", $$array[$i]->@*);
        $line .= "\n";
        print {$$self{fh}->$*} $line;
    }
}
1;
