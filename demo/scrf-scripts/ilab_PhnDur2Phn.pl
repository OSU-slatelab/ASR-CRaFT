#!/usr/bin/perl

#use Getopt::Std;
#getopt('dg');
#
#$usage="$0 -d maximum_duration -g group_size(default=1)";
#
#if ("$opt_d" eq "") {
#  die "$usage\n";
#}
#
#if ("$opt_d" ne "") {
#  #number of durations = |{1,2,3,...,$maxDur,>$maxDur}| = $maxDur + 1
#  $maxDur = $opt_d;
#  $numDur = $maxDur + 1;
#}
#
#$groupSize = 1;
#if ("$opt_g" ne "") {
#  $groupSize = $opt_g;
#}
#{
#use integer;
#$maxDurGroup = $maxDur / $groupSize;
#}
#$numDurGroup = $maxDurGroup + 1;
#
#while (<>) {
#  chomp;
#  ($utt, $seg, $bgn, $end, $phn) = split('\s+',$_);
#  $dur = $end - $bgn;
#  {
#  use integer;
#  $segDurGroup = ($dur - 1) / $groupSize + 1; 
#  }
#  if ($segDurGroup > $maxDurGroup) {
#      $segDurGroup = $maxDurGroup + 1;
#  }
#  $lab = $phn * $numDurGroup + $segDurGroup - 1;
#  print "$utt $seg $lab\n";
#}


use Getopt::Std;
getopt('m');

$usage="$0 -m Phn_Dur_to_Lab_mapping_file ";

if ("$opt_m" eq "") {
  die "$usage\n";
}

open MAPFILE,$opt_m;

$header=<MAPFILE>;
chomp($header);
($tag, $NUML)=split('\s+',$header);
if ($tag ne "NUML") {die "Error: The header format of the mapping file is incorrect: missing NUML.";}

$header=<MAPFILE>;
chomp($header);
($tag, $splitLongSeg)=split('\s+',$header);
if ($tag ne "SPLITLONGSEG") {die "Error: The header format of the mapping file is incorrect: missing SPLITLONGSEG.";}

$header=<MAPFILE>;
chomp($header);
($tag, $maxDur)=split('\s+',$header);
if ($tag ne "MAXDUR") {die "Error: The header format of the mapping file is incorrect: missing MAXDUR.";}

$header=<MAPFILE>;
chomp($header);
($tag, $groupSize)=split('\s+',$header);
if ($tag ne "GROUPSIZE") {die "Error: The header format of the mapping file is incorrect: missing GROUPSIZE.";}

$header=<MAPFILE>;
chomp($header);
($tag, $CUTOFF)=split('\s+',$header);
if ($tag ne "CUTOFF") {die "Error: The header format of the mapping file is incorrect: missing CUTOFF.";}

$counter=0;
%lab2Phn = (); 
while (<MAPFILE>) {
  chomp;
  ($phn,$dur,$lab)=split('\s+',$_);
  $lab2Phn{$lab} = $phn;
  $counter++;
}
close MAPFILE;
if ($counter != $NUML) {
  die "Error: The PhnDurLab mapping file claims $NUML elements; $counter elements found\n";
}

while (<>) {
  chomp;
  ($utt, $seg, $lab) = split('\s+',$_);
  if (exists $lab2Phn{$lab}) {
    $phn = $lab2Phn{$lab};
  } else {
    die "Error: Cannot find a phone for this label \"$_\"";
#    $phn = "<UNK>";
  }
  print "$utt $seg $phn\n";
}
