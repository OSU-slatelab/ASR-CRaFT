#!/usr/bin/env bash

set -e
set -o pipefail

if [ $# -ne 2 ]; then
  echo "Usage: $0 <config file> <exp dir>"
  exit 1
fi

script_dir=`dirname $0`
script_dir=`readlink -f $script_dir`

conf=$1
exp_dir=$2

if [ ! -f $conf ]; then
  echo "Error: cannot access the config file: $conf" 1>&2
  exit 1
fi
conf=`readlink -f $conf`

mkdir -p $exp_dir

cd $exp_dir
ln -s $conf conf
ln -s $script_dir scripts

cmd_path=/u/hey/src/kaldi/myrecipes/gen
ln -s $cmd_path/cmd_slurm.sh
ln -s $cmd_path/cmd_runpl.sh
ln -s cmd_slurm.sh cmd.sh
ln -s $cmd_path/utils myutils
