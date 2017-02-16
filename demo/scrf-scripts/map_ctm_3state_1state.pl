#!/usr/bin/env perl

use strict;
use warnings;

my $last_file = "";
my $last_chan = "";
my $last_start = "";
my $last_phn = "";
my $last_state = 0;
my $last_phn_len = 0;
while (<>) {
  chomp;
  if ($_ =~ m/^(\S+) (\S+) (\d+\.*\d*) (\d+\.\d*) (\S+)$/) {
    my $file = $1;
    my $chan = $2;
    my $start = $3;
    my $dur = $4;
    my $lab = $5;
    my $state = chop $lab;
    my $phn = $lab;
    if ($state !~ m/\d/) {
      die "Error: unrecognized state ($state) for the label: $lab\n";
    }
    if ($file ne $last_file) {
      if ($last_file ne "") {
        print "$last_file $last_chan $last_start $last_phn_len $last_phn\n";
      }
      $last_file = $file;
      $last_chan = $chan;
      $last_start = $start;
      $last_phn = $phn;
      $last_state = $state;
      $last_phn_len = 0;
    }
    if ($phn ne $last_phn || $state < $last_state) {
      print "$last_file $last_chan $last_start $last_phn_len $last_phn\n";
      $last_start = $start;
      $last_phn = $phn;
      $last_state = $state;
      $last_phn_len = 0;
    }
    $last_phn_len += $dur;
  } else {
    die "Error: unrecognized pattern: $_\n";
  }
}
if ($last_file ne "") {
  print "$last_file $last_chan $last_start $last_phn_len $last_phn\n";
}
