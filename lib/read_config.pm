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

use FindBin qw($Bin);
use lib "$Bin/../lib";

use zzIO;
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


sub read_config_yaml {
    my ($file) = @_;
    my $config = LoadFile($file);
    # remove empty item
    foreach my $key (keys %$config) {
        delete $config->{$key} if ref($config->{$key}) eq '' and $config->{$key} eq '';
    }
    foreach my $bin (qw/stmsa bcftools vg tabix muscle3 
                            mafft minigraph stretcher seqs2identity/) {
        my $path = $config->{$bin};
        next unless defined $path and $path ne '';
        my $bin_path_new = "$Bin/../$path";
        if (-e $path and -x $path) {
            $config->{$bin} = File::Spec->rel2abs($path);
            next;
        }
        if (-e $bin_path_new and -x $bin_path_new) {
            $config->{$bin} = File::Spec->rel2abs($bin_path_new);
            next;
        } else {
            die "Error: $bin not found in $path or $bin_path_new or is not excutable\n";
        }
    }
    return $config;
}

1;
