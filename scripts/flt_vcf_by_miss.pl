#!/usr/bin/env perl
#===============================================================================
#
#         FILE: 04-05.MISS_SNP_ALL.pl
#
#        USAGE: ./04-05.MISS_SNP_ALL.pl
#
#  DESCRIPTION:  
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zeyu Zheng (Lanzhou University), zhengzy2014@lzu.edu.cn
# ORGANIZATION:
#      VERSION: 1.0
#      CREATED: 10/24/2021 07:23:37 PM
#     REVISION: ---
#===============================================================================


use strict;
use warnings;
use utf8;
use v5.10;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;
use zzIO;
use MCE::Flow;
use MCE::Candy;

my($in_vcf, $opt_help, $out_vcf);
my $threshold_miss = 0;
my $thread=1;
GetOptions (
	'help|h!' => \$opt_help,
	'vcf=s' => \$in_vcf,
	'out=s' => \$out_vcf,
	't=i' => \$thread,
	'threshold_miss=s' => \$threshold_miss,
);


die "usage: perl xxx.pl [--threshold_miss 0(no miss allowed)] --vcf raw_snp.vcf.gz [-t 1] --out out.vcf\n" if ($opt_help || ! $in_vcf || ! $out_vcf);


print STDERR "\n** Now in node: ";
print STDERR`hostname`;


my $INVCF = open_in_fh($in_vcf);
my $O = open_out_fh($out_vcf);
my @idline;
print STDERR "Processing snp_vcf: $in_vcf\n";
while ( <$INVCF> ){
	if ( /^##/ ){
		print $O $_;
		next;
	}
	if (/^#/){
		chomp;
		print $O $_."\n";
		if (/^#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+\w+/){
			@idline=split(/\s+/,$_);
			last;
		}else{
			print STDERR "$_\n";
			last;
		}
	}
}
die "NO idline: \n" if ! @idline;

if ( $thread == 1) {
	while ( <$INVCF> ){
		my $ret = &cal($_);
		say $O $ret if defined $ret;
	}
	exit;
}


sub flt_mce {
	my ( $mce, $chunk_ref, $chunk_id ) = @_;
	my $result = &cal($$chunk_ref[0]);
	if (defined $result) {
		$mce->gather($chunk_id, $result, "\n") 
	} else {
		$mce->gather($chunk_id);
	}
}

mce_flow_f {
    chunk_size => 1,
    max_workers => $thread,  
    #max_workers => 1,
	#use_slurpio => 0,
    gather => MCE::Candy::out_iter_fh($O)
    },
    \&flt_mce,
    $INVCF;


say STDERR 'run ok';

exit;


sub cal {
	my $in = $_;
	chomp $in;
	my @a = split(/\s+/,$in);
	my $miss=0;
	foreach my $i (9..@a-1) {
		my $id = $idline[$i];
		$miss++ if $a[$i] =~/^\.\/\./;
	}
	my $percent = $miss/(@a-9);
	#print "$miss\n";
	return join("\t",@a) if $percent <= $threshold_miss;
	return undef;
}

