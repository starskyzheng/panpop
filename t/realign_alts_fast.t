use strict;
use warnings;

use Test::More;

BEGIN {
    eval { require YAML; 1 } or plan skip_all => 'YAML.pm not available';
}

use FindBin;
use lib "$FindBin::Bin/../lib";

use realign_alts;

my @cases = (
    {
        name => 'simple_substitution',
        alts => [ 'ACGT', 'ATGT', 'ACGT' ],
    },
    {
        name => 'with_internal_gap',
        alts => [ 'ACGT', 'A-GT', 'AC-T' ],
    },
    {
        name => 'leading_gap',
        alts => [ '-CGT', 'ACGT', '-C-T' ],
    },
    {
        name => 'trailing_gap_with_ext',
        alts => [ 'ACG-', 'A-G-', 'ACGT' ],
        ext  => [ 'T', 'A' ],
    },
);

realign_alts::_set_fast_helpers_enabled(0);
my %baseline;
foreach my $case (@cases) {
    my $alts = $case->{alts};
    my $max_alts = @$alts - 1;
    my ($muts, $len) = realign_alts::alt_alts_to_muts($alts, $max_alts, $case->{ext});
    ok(defined $muts, "slow path produced result for $case->{name}");
    $baseline{$case->{name}} = [$muts, $len];
}

if (realign_alts::_fast_helpers_available()) {
    ok(realign_alts::_set_fast_helpers_enabled(1), 'enabled native helpers');
    foreach my $case (@cases) {
        my $alts = $case->{alts};
        my $max_alts = @$alts - 1;
        my ($muts_fast, $len_fast) = realign_alts::alt_alts_to_muts($alts, $max_alts, $case->{ext});
        my ($muts_slow, $len_slow) = @{ $baseline{$case->{name}} };
        is_deeply($muts_fast, $muts_slow, "fast matches slow for $case->{name}");
        is($len_fast, $len_slow, "alignment length matches for $case->{name}");
    }
} else {
    diag('Native helpers not available; fast-path comparisons skipped');
}

done_testing();
