#
#===============================================================================
#
#         FILE: read_config.pm
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zeyu Zheng (LZU), zhengzy2014@lzu.edu.cn
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 12/01/2021 10:11:37 PM
#     REVISION: ---
#===============================================================================

package read_config;
use strict;
use warnings;
use v5.24;
use File::Spec;
no warnings 'experimental::smartmatch';

use FindBin qw($Bin);
use lib "$Bin/../lib";

#use zzIO;
use YAML qw(LoadFile);

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
            read_config_yaml
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw(read_config_yaml);


my @config_files = (
    "$Bin/../configs/software.yaml",
    "$Bin/../configs/base.yaml",
);

my @software_all = qw/stmsa bcftools vg tabix muscle3 bgzip famsa  mafft minigraph/;

my @software_musthave = qw/stmsa bcftools vg tabix muscle3 bgzip famsa/;

sub read_config_yaml {
    #my (@files) = @_;
    my @files = @config_files;
    my %config;
    foreach my $file(@files) {
        my $config_tmp = LoadFile($file);
        # remove empty item
        foreach my $key (keys %$config_tmp) {
            next if ref($config_tmp->{$key}) eq '' and $config_tmp->{$key} eq '';
            $config{$key} = $config_tmp->{$key};
        }
    }

    foreach my $software_name (@software_all) {
        my $path = $config{$software_name};
        next unless defined $path and $path ne '';
        my $bin_path_new = "$Bin/../$path";
        my $which_path = `which $path 2>/dev/null`;
        chomp $which_path;
        if ($which_path and -e $which_path and -x $which_path) {
            $config{$software_name} = File::Spec->rel2abs($which_path);
            next;
        } elsif (-e $path and -x $path) {
            $config{$software_name} = File::Spec->rel2abs($path);
            next;
        } elsif (-e $bin_path_new and -x $bin_path_new) {
            $config{$software_name} = File::Spec->rel2abs($bin_path_new);
            next;
        } elsif($software_name ~~ @software_musthave) {
            die "Error: $software_name not found ( $path ) or $bin_path_new or is not excutable\n";
        }
    }
    return \%config;
}

1;
