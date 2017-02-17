#!/usr/bin/env bash

function HELP {
  echo "Usage: $0 [options] <dataset [cv|dt|core|enh]> <phn_ilab_ascii>"
  echo "options: -s symlist (required)"
  echo "         -o olist (required)"
  echo "         -n num_state (optional, [1(default)|3|...]"
  echo "         -d dur_ilab_ascii (required for segmental models; not needed for frame-based models)"
  exit 1
}

if [ $# -eq 0 ]; then
  HELP
fi

echo "$0 $@"

numstate=1

### Start getopts code ###

#Parse command line flags
#If an option should be followed by an argument, it should be followed by a ":".
#Notice there is no ":" after "h". The leading ":" suppresses error messages from
#getopts. This is required to get my unrecognized option code to work.

while getopts :s:o:n:d:h FLAG; do
  case $FLAG in
    s)  #set option "s"
      symlist=$OPTARG
      echo "-s used: $OPTARG"
      echo "symlist = $symlist"
      ;;
    o)  #set option "o"
      olist=$OPTARG
      echo "-o used: $OPTARG"
      echo "olist = $olist"
      ;;
    n)  #set option "n"
      numstate=$OPTARG
      echo "-n used: $OPTARG"
      echo "numstate = $numstate"
      ;;
    d)  #set option "d"
      dur=$OPTARG
      echo "-d used: $OPTARG"
      echo "dur = $dur"
      ;;
    h)  #show help
      HELP
      ;;
    \?) #unrecognized option - show help
      echo -e \\n"Option -$OPTARG not allowed."
      HELP
      #If you just want to display a simple error message instead of the full
      #help, remove the 2 lines above and uncomment the 2 lines below.
      #echo -e "Use $0 -h to see the help documentation."\\n
      #exit 2
      ;;
  esac
done

shift $((OPTIND-1))  #This tells getopts to move on to the next argument.

### End getopts code ###

if [ $# -ne 2 ]; then
  echo "Error: no input file."
  HELP
fi

dataset=$1
ilab=$2

if [ ! -f $ilab ]; then
  echo "Error: file not exit: $ilab"
  exit 1
fi

script_dir=`dirname $0`

dir=`dirname $ilab`
prefix=`basename $ilab .phn.lab.ascii`

if [ ! -z $dur ]; then
  if [ ! -f $dur ]; then
    echo "Error: file does not exist: $dur"
    exit 1
  fi
  dur_opt="-d $dur"
else
  dur_opt=
fi

echo "ilab2ctm.pl -s $symlist -o $olist $dur_opt $ilab > $dir/$prefix.ctm"
ilab2ctm.pl -s $symlist -o $olist $dur_opt $ilab > $dir/$prefix.ctm
out_ctm=$prefix.ctm
if [ "$numstate" == 3 ]; then
  $script_dir/map_ctm_3state_1state.pl < $dir/$prefix.ctm > $dir/$prefix.1state.ctm
  out_ctm=$prefix.1state.ctm
fi

pushd $dir

$script_dir/timit_norm_trans.pl -i $out_ctm -m $script_dir/timit_score_new/phones.60-48-39.map -from 48 -to 39 > $prefix.map39.ctm
ln -fsn $prefix.map39.ctm ctm_39phn

cp -p $script_dir/timit_score_new/glm_39phn ./

ref=$script_dir/timit_score_new/$dataset/stm_39phn
if [ ! -f $ref ]; then
  echo "$0 Error: reference file does not exist: $ref" 1>&2
  exit 1
fi
cp -p $ref ./

hubscr_dir=$script_dir/sctk/bin
$hubscr_dir/hubscr.pl -p $hubscr_dir -V -l english -h hub5 -g glm_39phn -r stm_39phn ctm_39phn

tail -n 10 ctm_39phn.filt.sys

popd
