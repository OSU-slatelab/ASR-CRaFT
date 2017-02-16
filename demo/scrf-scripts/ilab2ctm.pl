#!/usr/bin/perl

use Getopt::Std;
#getopts('ts:o:d:');
getopts('s:o:d:');

#$usage="$0 -s symlist_file -o olist_file -d dur_ilab_ascii_file -t";
$usage="$0 -s symlist_file -o olist_file -d dur_ilab_ascii_file";
# dur_ilab_ascii_file should have exactly the same number of lines as the input phone_ilab_ascii_file, but indicating the duration of each corresponding phone segment

if ("$opt_s" eq "" || "$opt_o" eq "") {
  die "$usage\n";
}

#$showt=0;
#if ("$opt_t" ne "") {
#  $showt=1;
#}

open SYMS,$opt_s;
$symc=<SYMS>;
chomp($symc);
$counter=0;
undef @symList;
while (<SYMS>) {
  chomp;
  ($symList[$counter],$junk)=split('\s+',$_);
  $counter++;
}
close SYMS;
if ($counter != $symc) {
  die "Symbol List File claims $symc elements; $counter elements found\n";
}

open OLIST,$opt_o;
undef @olist;
$c=0;
while (<OLIST>) {
  chomp;
  $olist[$c]=$_;
  $c++;
}
close OLIST;

if ("$opt_d" ne "") {
  open(DUR, "$opt_d") || die "Error: cannot open the duration ilab ascii file: $opt_d\n";
}

$cur_utt=-1;
$last_lab=-1;
$start=0; #inclusive
$end=0;   #exclusive
$phn_len=0;
while (<>) {
  chomp;
  ($utt, $frm, $lab) = split('\s+',$_);
  if ("$opt_d" ne "") {
    $dur_line = <DUR>;
    ($d_utt, $d_frm, $dur) = split('\s+',$dur_line);
    $d_utt eq $utt && $d_frm eq $frm || die "Error: phone ilab file and duration ilab file does not match:\n$_\n$dur_line\n";
  } else {
    $dur=1;
  }
  if ($utt != $cur_utt) {
    if ($cur_utt >=0) {
      #print_label($utt_name, $start, $phn_len, $last_lab, $showt);
      print_label($utt_name, $start, $phn_len, $last_lab);
    }
    $start = 0;
    $last_lab=$lab;
    $phn_len=0;
    $cur_utt=$utt;
    $utt_name=$olist[$utt];
    $utt_name =~ s/^dr\d_//;
    $utt_name =~ s/.lat$//;
  } elsif ($lab != $last_lab) {
    #print_label($utt_name, $start, $phn_len, $last_lab, $showt);
    print_label($utt_name, $start, $phn_len, $last_lab);
    $start += $phn_len;
    $last_lab=$lab;
    $phn_len=0;
  }
  $phn_len += $dur;
}
#print_label($utt_name, $start, $phn_len, $last_lab, $showt);
print_label($utt_name, $start, $phn_len, $last_lab);

if ("$opt_d" ne "") {
  close DUR;
}

sub print_label {
  #my ($utt_name, $start, $phn_len, $last_lab, $showt) = @_;
  my ($utt_name, $start, $phn_len, $last_lab) = @_;
  print "$utt_name 1 ";
  #if ($showt == 1) {
    printf("%.2f %.2f ", $start/100, $phn_len/100);
  #}
  print "$symList[$last_lab]\n";
}
