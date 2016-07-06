#!/bin/bash

# Copyright 2014  The Ohio State University (Author: Yanzhang He)
# Apache 2.0

# Extract and store DNN bottleneck features. 

. ./cmd.sh ## You'll want to change cmd.sh to something that will work on your system.
           ## This relates to the queue.

. ./path.sh ## Source the tools/utils (import the queue.pl)

# Config:
srcdata=data  # set it to data-fmllr-tri3 for speaker-adapted triphone state input features
dnndir=exp/dnn_mono_pretrain-dbn_dnn   # set it to dnn4_pretrain-dbn_dnn for speaker-adapted triphone state dnn
data_bn=data-bn-dnn-mono
use_prior=false
htk_save=true
remove_last_components=0
stage=0 # resume training with --stage=N
# End of config.
. utils/parse_options.sh || exit 1;
#

if [ $stage -le 0 ]; then
  # Store bottleneck features, so we can train on them easily,
  # train
  dir=$data_bn/train
  if [ ! -f $dir/.done ]; then
    mysteps/nnet/make_bn_feats.sh --nj 10 --cmd "$train_cmd" \
       --remove-last-components $remove_last_components --htk-save $htk_save --use-prior $use_prior \
       $dir $srcdata/train $dnndir $dir/log $dir/data && touch $dir/.done &
    # split the data : 90% train 10% cross-validation (held-out)
    #utils/subset_data_dir_tr_cv.sh $dir ${dir}_tr90 ${dir}_cv10 || exit 1
    #|| exit 1
  fi
  # dev
  dir=$data_bn/dev
  if [ ! -f $dir/.done ]; then
    mysteps/nnet/make_bn_feats.sh --nj 10 --cmd "$train_cmd" \
       --remove-last-components $remove_last_components --htk-save $htk_save --use-prior $use_prior \
       $dir $srcdata/dev $dnndir $dir/log $dir/data && touch $dir/.done &
    #|| exit 1
  fi
  # test
  dir=$data_bn/test
  if [ ! -f $dir/.done ]; then
    mysteps/nnet/make_bn_feats.sh --nj 10 --cmd "$train_cmd" \
       --remove-last-components $remove_last_components --htk-save $htk_save --use-prior $use_prior \
       $dir $srcdata/test $dnndir $dir/log $dir/data && touch $dir/.done &
    #|| exit 1
  fi
  # enhanced
  dir=$data_bn/enh
  if [ ! -f $dir/.done ]; then
    mysteps/nnet/make_bn_feats.sh --nj 10 --cmd "$train_cmd" \
       --remove-last-components $remove_last_components --htk-save $htk_save --use-prior $use_prior \
       $dir $srcdata/enh $dnndir $dir/log $dir/data && touch $dir/.done &
    #|| exit 1
  fi
fi

wait
echo Success
exit 0

