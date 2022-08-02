#!/usr/bin/env perl
#!/usr/bin/perl
package zzIO;

#use strict;
#use warnings;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use read_config qw/read_config_yaml/;
use Carp;

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
            open_in_fh
            open_out_fh
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw(open_in_fh open_out_fh);


#BEGIN {}

my $config = read_config_yaml("$Bin/../config.yaml");
my $bgzip = $$config{bgzip} or die "bgzip not defined in config file!";
my $pigztest = `which pigz 2>/dev/null`; chomp $pigztest;
my $pigz = $pigztest ? $pigztest : $bgzip;

sub open_in_fh($;$) {
    my $in = $_[0] or confess "zzIO::open_in_fh: given no filename";
    my $parallel = $_[1] // 1;
    my $fh;
    if ($in eq '-') {
        $fh = \*STDIN;
        return $fh;
    } elsif ($in eq *main::STDIN) {
        return \*STDIN;
    } elsif (ref $in eq 'GLOB' ) {
        return $in;
    } elsif ($in =~ /.gz$/) {
        confess "$in not exists! " unless -e $in;
        open ($fh,"zcat $in|");
        if ($parallel>1) {
            open ($fh,"$pigz --stdout --decompress --processes $parallel $in|");
        } else {
            open ($fh,"zcat $in|");
        }
        return $fh;
    }else {
        open ($fh,"<$in") or confess "cannot open_in_fh with $in: $!";
        return $fh;
    }
    confess;
}


sub open_out_fh($;$) {
    my $out = $_[0] or confess "zzIO::open_out_fh: given no filename";
    my $parallel = $_[1] // 2;
    confess "'$parallel' is not a Integer" unless $parallel=~/^\d+$/;
    my $fh;
    if ($out =~ /\.gz$/) {
    #if ($out =~ /\.vcf\.gz$/) {
        #$parallel ? open ($fh,"| $pigz -p $parallel -c - > $out") : open ($fh,"| $bgzip -c > $out");
        open ($fh,"| $bgzip -c -@ $parallel > $out");
        return $fh;
    #}elsif ($out =~ /\.gz$/) {
    #    $parallel > 1 ? open ($fh,"| pigz -p $parallel -c - > $out") : open ($fh,"| gzip -c - > $out");
    #    return $fh;
    }elsif (ref $out eq 'GLOB' ) {
        return $out;
    }elsif ($out eq '-') {
        $fh = \*STDOUT;
        return $fh;
    }else {
        open ($fh,">$out") or confess  "cannot open fh with $out: $!";
        return $fh;
    }
    confess;
}

1;
