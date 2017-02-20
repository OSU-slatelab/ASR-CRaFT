#!/bin/bash

USAGE="Usage: $0 lrstart=[start learning rate] lrend=[end learning rate] lrstep=[learning rate step] wt_pre=[weight directory, optional]"
#exit;

lrstart=
lrend=
lrstep=
lr_range=
etastart=
etaend=
etastep=
eta_range=
wt_pre=
dataset=

parsed_opts= `getopt -o h --long lrstart::,lrend::,lrstep::,lr_range::,\
etastart::,etaend::,etastep::,eta_range::,wt_pre::,dataset:: \
-n "$0" -- "$@"`
eval set -- "$parsed_opts"

while true; do
    case "$1" in
	-h)
	    echo "$USAGE"
	    exit ;;
	--lrstart)
	    case "$2" in
		"") shift 2 ;;
		*) lrstart=$2; shift 2 ;;
	    esac ;;
	--lrend)
	    case "$2" in
		"") shift 2 ;;
		*) lrend=$2; shift 2 ;;
	    esac ;;
	--lrstep)
	    case "$2" in
		"") shift 2 ;;
		*) lrstep=$2; shift 2 ;;
	    esac ;;
	--lr_range)
	    case "$2" in
		"") shift 2 ;;
		*) lr_range=$2; shift 2 ;;
	    esac ;;
	--etastart)
	    case "$2" in
		"") shift 2 ;;
		*) etastart=$2; shift 2 ;;
	    esac ;;
	--etaend)
	    case "$2" in
		"") shift 2 ;;
		*) etaend=$2; shift 2 ;;
	    esac ;;
	--etastep)
	    case "$2" in
		"") shift 2 ;;
		*) etastep=$2; shift 2 ;;
	    esac ;;
	--eta_range)
	    case "$2" in
		"") shift 2 ;;
		*) eta_range=$2; shift 2 ;;
	    esac ;;
	--wt_pre)
	    case "$2" in
		"") shift 2 ;;
		*) wt_pre=$2; shift 2 ;;
	    esac ;;
	--dataset)
	    case "$2" in
		"") shift 2 ;;
		*) dataset=$2; shift 2 ;;
	    esac ;;
	--)
	    shift; break;;
	*)
	    echo "Unknown option $1"; exit 1 ;;
    esac
done
	

if [ -z "$lrstart" ] && [ -z "$lr_range" ]; then
  lr_range="0.0005 `seq 0.001 0.001 0.009 | tr '\n' ' '` `seq 0.01 0.01 0.09 | tr '\n' ' '` `seq 0.1 0.1 0.9 | tr '\n' ' '`"
fi

if [ -z "$lrend" ]; then
  lrend=$lrstart
fi

if [ -z "$lrstep" ]; then
  lrstep=0.1
fi

if [ -z "$lr_range" ]; then
  lr_range=`seq $lrstart $lrstep $lrend`
fi

if [ -z "$etastart" ] && [ -z "$eta_range" ]; then
  eta_range="`seq 0.001 0.001 0.009 | tr '\n' ' '` `seq 0.01 0.01 0.09 | tr '\n' ' '` `seq 0.1 0.1 0.9 | tr '\n' ' '` `seq 1.0 1.0 9.0 | tr '\n' ' '` `seq 10.0 10.0 90.0 | tr '\n' ' '`"
fi

if [ -z "$etaend" ]; then
  etaend=$etastart
fi

if [ -z "$etastep" ]; then
  etastep=0.1
fi

if [ -z "$eta_range" ]; then
  eta_range=`seq $etastart $etastep $etaend`
fi

if [ -z "$wt_pre" ]; then
  wt_pre=./weights
  echo "wt_pre not defined, using default of $wt_pre"
fi

if [ -z "$dataset" ]; then
  dataset="cv core enh"
  echo "dataset not defined, using default of $dataset"
fi

for ds in $dataset; do

  allresult="wt_pre=$wt_pre\nlr_range=$lr_range\n\n"
  for lr in $lr_range; do
    allresult="$allresult\nlr=$lr\n"
    decode_dir=$wt_pre.lr$lr/${ds}_out.lr$lr
    if [ -d $decode_dir ]; then
      result=`grep -r "Sum/Avg" $decode_dir/*/*.sys | sort`
      printf "%b\n" "$result" > $wt_pre.lr$lr/${ds}_results.all.txt
      allresult="$allresult$result\n"
    fi
  done

  allresult=$allresult"\neta_range=$eta_range\n\n"
  for eta in $eta_range; do
    allresult="$allresult\neta=$eta\n"
    decode_dir=$wt_pre.ag_eta$eta/${ds}_out.ag_eta$eta
    if [ -d $decode_dir ]; then
      result=`grep -r "Sum/Avg" $decode_dir/*/*.sys | sort`
      printf "%b\n" "$result" > $wt_pre.ag_eta$eta/${ds}_results.all.txt
      allresult="$allresult$result\n"
    fi
  done

  if [ -z "$lrstart" ] && [ -z "$etastart" ]; then
    outfile=$wt_pre.all_lr_eta.${ds}_results.all.txt
  elif [ ! -z "$lrstart" ]; then
    outfile=$wt_pre.lr$lrstart-$lrstep-$lrend.${ds}_results.all.txt
  elif [ ! -z "$etastart" ]; then
    outfile=$wt_pre.eta$etastart-$etastep-$etaend.${ds}_results.all.txt
  else
    echo "Error: cannot have both lrstart and etastart"
    exit 1
  fi

  printf "%b\n" "$allresult" > $outfile

  echo outfile: $outfile

done
