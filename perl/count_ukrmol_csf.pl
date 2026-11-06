#!/usr/bin/env perl

# -- cont target and l2  CSFs given a UKRmol+ model file using the MAS approach.
#    Use this script as
#
#      perl -I<UKRMOL_SCRIPTS>/lib <THIS_SCRIPT> model_file.pl
#

use strict;
use warnings;
use MultiSpace;
use File::Basename qw(basename);
use File::Spec;

my $program = basename($0);

our %model;
my $file = shift or die "usage: $program model.pl\n";
do File::Spec -> rel2abs($file) or die die "cannot load $file: ", ($@ || $!), "\n";

die "model file is missing a 'MAS' entry" unless $model{'MAS'};
my $nelec = $model{'nelectrons'};
my $sym = $model{'symmetry'};
my $gsm = $model{'groundStateSpinMultiplicity'};
my $type = $model{'use_MAS'} || 1;

my %irreps = (
  C1  => [ qw(A) ],
  Ci  => [ qw(Ag Au) ],
  C2  => [ qw(A B) ],
  Cs  => [ qw(Ap App) ],
  C2v => [ qw(A1 B1 B2 A2) ],
  C2h => [ qw(Ag Au Bu Bg) ],
  D2  => [ qw(A  B3 B2 B1) ],
  D2h => [ qw(Ag B3u B2u B1g B1u B2g B3g Au) ],
);


sub report_csf_count {

  my ($label, $space, $ne, $spinmult) = @_;
  my %n = $space -> dimension($ne, $spinmult);
  my $total = 0;

  printf "%-20s", "$label (spinmult $spinmult):";

  for my $i (0 .. $#{$irreps{$sym}}) {
    my $key = (grep { exists $n{$_} } ($i, $i+1, $irreps{$sym}[$i]))[0];
    my $val = defined $key ? $n{$key} : 0;
    printf "  %-5s %8d", $irreps{$sym}[$i], $val;
    $total += $val;
  }

  printf "    total %10d\n", $total;

}

print "\nmodel: $file\nsymmetry: $sym\nnelectrons: $nelec\nuse_MAS: $type\n\n";

my $target = MultiSpace -> new({ symmetry => $sym, MAS => $model{'MAS'}, constraints => $model{'constraints'}, type => $type });

report_csf_count("target", $target, $nelec, $gsm);

my $l2mas = $model{'l2_MAS'} || $model{'MAS'};
my $l2 = MultiSpace -> new({ symmetry => $sym, MAS => $l2mas, constraints => $model{'l2_constraints'}, type => $type });

print "\n", ($model{'l2_MAS'} ? "L2 from l2_MAS" : "L2 = target MAS + 1 electron"), "\n";

for my $spinmult (grep  { $_ > 0 } ($gsm - 1, $gsm + 1)) {
  report_csf_count("L2", $l2, $nelec + 1, $spinmult);
}

print "\nMOCSF per symmetry ~ L2 count + [ contiuum orbitals x target states ]\n";
