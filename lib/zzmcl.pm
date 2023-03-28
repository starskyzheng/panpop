#!/usr/local/bin/perl -w

#  (C) Copyright 2006, 2007, 2008, 2009 Stijn van Dongen
#
#  This file is part of MCL.  You can redistribute and/or modify MCL under the
#  terms of the GNU General Public License; either version 3 of the License or
#  (at your option) any later version.  You should have received a copy of the
#  GPL along with MCL, in the file COPYING.

# mod from: http://www.micans.org/mcl/scripts/minimcl

package zzmcl;

require strict;
require warnings;
use Data::Dumper;
use v5.10;

sub explain {
    print <<EOH;
purpose:
   A small mcl implementation for educational purposes.  It is written
   in moderately terse perl.

implementation:
   It is hash based, which implies that we get sparse matrices easily but at
   the cost of using hashes.  The hash-based matrices only store non-zero
   entries.

   The code is pretty straightforward. The interpretation routine implements
   the mapping as described in the publications referenced in the (maxi) mcl
   manual.

bonus:
   Since the implementation is hash based you can use any type of labels, not
   necessarily numbers.

Usage:
   minimcl [--I=<num>] [--verbose] LABEL-INPUT

   This means --I=<num> is optional (with 2.0 the default) and so is --verbose.
   LABEL-INPUT should be a file name or stream (STDIN) where each line is of
   the form
            LABEL1 LABEL2 NUMBER
   or
            LABEL1 LABEL2
EOH
}


sub new {
    my ($class, $args) = @_;
    my %args = %$args;
    my $self={
        verbose => $args{verbose},
        I => $args{I},
        mx => {},
        #pair2indentities => {},
    };
    bless($self, $class || ref($class));
    return $self;
}

sub addEdge {
   my ($self, $x, $y, $val) = @_;
   if(!defined $val) {
      $val = 1.0;
   } elsif ($val !~ /^[0-9]/) {
      $val = 1.0;
   } elsif($val < 0) {
      $val = 0;
   } elsif($val > 1) {
      $val = 1;
   }
   $val = $val + 0;
   $self->{mx}{$x}{$y} = $val;
   $self->{mx}{$y}{$x} = $val;
}

sub run {
   my ($self) = @_;
   my $mx = $self->{mx};
   my $I = $self->{I};
   my $verbose = $self->{verbose};
   #die Dumper $mx;
   matrix_add_loops($mx);
   matrix_make_stochastic($mx);
   #die Dumper $mx;
   matrix_dump( $mx, 3, "start" ) if $verbose;
   my ( $cl, $limit ) = $self->mcl( $mx, $I );
   my $group2arrays = group2arrays($cl);
   #die Dumper $group2arrays;
   return $group2arrays;
   #matrix_dump( $limit, 1, "limit" ) if $verbose;
   #matrix_dump( $cl, 0, "clustering" );
}

sub group2arrays {
   my ($cl) = @_;
   my @group2arrays;
   my %ids;
   for my $n ( sort {$a<=>$b} keys %$cl ) {
      my $group = $cl->{$n};
      my @ids_now_all = sort keys %$group;
      my @ids_now;
      foreach my $id (@ids_now_all) {
         if(! exists $ids{$id}) {
            push @ids_now, $id;
            $ids{$id} = 1;
         }
      }
      $group2arrays[$n] = \@ids_now;
   }
   return(\@group2arrays);
}

sub mcl {
    my ( $self, $mx, $I ) = @_;
    my $chaos = 1;
    my $ite = 1;
    my $verbose = $self->{verbose};
    while ( $chaos > 0.001 ) {
        my $sq = matrix_square($mx);
        #my $progress = sprintf "chaos %.5f ite %d", $chaos, $ite;
        #matrix_dump( $sq, 3, "X $progress" ) if $verbose;
        $chaos = matrix_inflate( $sq, $I );
        #matrix_dump( $sq, 3, sprintf "I $progress" ) if $verbose;
        #print STDERR "$progress\n" if !$verbose;
        $mx = $sq;
        $ite++;
    }
    my $cl = matrix_interpret($mx);
    return ( $cl, $mx );
}

# dangersign:
# can this yield a < b < c < a ?

sub cmpany { local $^W = 0; $a <=> $b || $a cmp $b }


sub matrix_dump {
    my ( $mx, $modes, $msg ) = @_;
    print "($msg\n";
    for my $n ( sort cmpany keys %$mx ) {
        my @nb
            = $modes & 2
            ? map { sprintf "%s:%.3f", $_, $mx->{$n}{$_}; }
            sort cmpany keys %{ $mx->{$n} }
            : map { sprintf "%s", $_; } sort cmpany keys %{ $mx->{$n} };
        local $" = "\t";
        if ( $modes & 1 ) {
            printf "%-20s%s\n", $n, "@nb";
        } else {
            print "@nb\n";
        }
    }
    print ")\n";
}


sub matrix_square {
    my ($mx) = @_;
    my $sq = {};
    my @nodes = keys %$mx;
    for my $n (@nodes) {
        $sq->{$n} = matrix_multiply_vector( $mx, $mx->{$n} );
    }
    return $sq;
}


sub matrix_multiply_vector {
    my ( $mx, $v ) = @_;
    my $w = {};
    for my $e ( keys %$v ) {
        my $val = $v->{$e};
        for my $f ( keys %{ $mx->{$e} } ) {
            $w->{$f} += $val * $mx->{$e}{$f};
        }
    }
    return $w;
}


sub matrix_make_stochastic {
    my ($mx) = @_;
    matrix_inflate( $mx, 1 );  # return value chaos is meaningless for

    # non stochastic input.
}


sub matrix_add_loops {
    my ($mx) = @_;
    for my $n ( keys %$mx ) {
        $mx->{$n}{$n} = 0 if defined( $mx->{$n}{$n} );
        my $max = vector_max( $mx->{$n} );
        $mx->{$n}{$n} = $max ? $max : 1;
    }
}


sub vector_max {
    my ($v) = (@_);
    my $max = 0;
    for my $n ( keys %$v ) {
        $max = $v->{$n} if $v->{$n} > $max;
    }
    return $max;
}


sub vector_sum {
    my ( $v, $p ) = (@_);
    my $sum = 0;
    for my $n ( keys %$v ) {
        $sum += $v->{$n}**$p;
    }
    return $sum;
}


sub matrix_inflate {  # prunes small elements as well.
    my ( $mx, $I ) = @_;
    my @nodes = keys %$mx;
    my $chaos = 0;
    for my $n (@nodes) {
        my $sum = 0;
        my $sumsq = 0;
        my $max = 0;
        for my $nb ( keys %{ $mx->{$n} } ) {
            if ( $mx->{$n}{$nb} < 0.00001 ) {
                delete( $mx->{$n}{$nb} );
                next;
            }
            $mx->{$n}{$nb}**= $I;
            $sum += $mx->{$n}{$nb};
        }
        if ($sum) {
            for my $nb ( keys %{ $mx->{$n} } ) {
                $mx->{$n}{$nb} /= $sum;
                $sumsq += $mx->{$n}{$nb}**2;  # sum x_i^2 over stochastic vector x
                $max = $mx->{$n}{$nb} if $max < $mx->{$n}{$nb};
            }
        }
        $chaos = $max - $sumsq if $max - $sumsq > $chaos;
    }
    return $chaos;  # only meaningful if input is stochastic
}

# assumes but does not check doubly idempotent matrix.
# can handle attractor systems of size < 10.
sub matrix_interpret {  # recognizes/preserves overlap.
    my ($limit) = @_;
    my $clusters = {};  # hash of arrayrefs.
    my $attrid = {};
    my $clid = 0;
    my $count = scalar(keys %$limit);
    $count = 1 if $count < 1;
    my $min = 0.5/$count ; # default 0.1;
    $min = 0.1 if $min > 0.1;
    for my $n ( keys %$limit ) {  # crude removal of small elements.
        for my $nb ( keys %{ $limit->{$n} } ) {
            delete $limit->{$n}{$nb} if $limit->{$n}{$nb} < $min;
            #################
        }
    }
    my $attr = { map { ( $_, 1 ) } grep { $limit->{$_}{$_} } keys %$limit };

    # _ contract 'connected attractors', assign cluster id.
    for my $a ( keys %$attr ) {
        next if defined( $attrid->{$a} );
        my @aa = ($a);
        while (@aa) {
            my @bb = ();
            for my $aa (@aa) {
                $attrid->{$aa} = $clid;
                push @bb, grep { defined( $attr->{$_} ) } keys %{ $limit->{$aa} };
            }
            @aa = grep { !defined( $attrid->{$_} ) } @bb;
        }
        $clid++;
    }
    for my $n ( keys %$limit ) {
        if ( !defined( $attr->{$n} ) ) {  # look at attractors
            for my $a ( grep { defined( $attr->{$_} ) } keys %{ $limit->{$n} } ) {
                $clusters->{ $attrid->{$a} }{$n}++;
            }
        } else {
            $clusters->{ $attrid->{$n} }{$n}++;
        }
    }
    return $clusters;
}

1;
