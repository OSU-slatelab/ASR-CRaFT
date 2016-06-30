#!/usr/bin/env bash

if [ $# -ne 3 ]; then
  echo "Usage: <dir> <train_outfile> <test_outfile>"
  echo "       This script assumes the posterior or bottleneck features were already generated through mylocal/prep_bn_feats.sh."
  echo "       The input <dir> to this script should be equal to the output directory in mylocal/prep_bn_feats.sh defined by --data_bn."
  exit 1
fi

dir=$1
train_outfile=$2
test_outfile=$3

script_dir=`dirname $0`
if [[ ! $script_dir =~ ^/ ]]; then
  script_dir=`pwd`/$script_dir;  # relative path -> absolute path
fi

if [ ! -d $dir ]; then
  echo "Error: directory not exist: $dir"
  exit 1
fi

pushd $dir
feacat -ip htk -op pf -lists -i $script_dir/timit_sisx_train.olist.kaldi_to_htk_to_pfile -o $train_outfile
feacat -ip htk -op pf -lists -i $script_dir/timit_sisx_test.neworder.olist.kaldi_to_htk_to_pfile -o $test_outfile
popd
