#!/bin/bash

# Copyright 2012-2014  Brno University of Technology (Author: Karel Vesely)
# Apache 2.0

# This example script trains a DNN on top of fMLLR features. 
# The training is done in 3 stages,
#
# 1) RBM pre-training:
#    in this unsupervised stage we train stack of RBMs, 
#    a good starting point for frame cross-entropy trainig.
# 2) frame cross-entropy training:
#    the objective is to classify frames to correct pdfs.
# 3) sequence-training optimizing sMBR: 
#    the objective is to emphasize state-sequences with better 
#    frame accuracy w.r.t. reference alignment.

. ./cmd.sh ## You'll want to change cmd.sh to something that will work on your system.
           ## This relates to the queue.

. ./path.sh ## Source the tools/utils (import the queue.pl)

# Config:
#gmmdir=exp/tri3
gmmdir=exp/mono
data_fmllr=data-fmllr-tri3
data=$data_fmllr
#pretraindir=exp/dnn4_pretrain-dbn
pretraindir=exp/dnn_fmllr_mono_pretrain-dbn
stage=0 # resume training with --stage=N
skip_smbr=true
# End of config.
. utils/parse_options.sh || exit 1;
#

#if [ $stage -le 0 ]; then
  ## Store fMLLR features, so we can train on them easily,
  ## test
  #dir=$data_fmllr/test
  #steps/nnet/make_fmllr_feats.sh --nj 10 --cmd "$train_cmd" \
  #   --transform-dir $gmmdir/decode_test \
  #   $dir data/test $gmmdir $dir/log $dir/data || exit 1
  ## dev
  #dir=$data_fmllr/dev
  #steps/nnet/make_fmllr_feats.sh --nj 10 --cmd "$train_cmd" \
  #   --transform-dir $gmmdir/decode_dev \
  #   $dir data/dev $gmmdir $dir/log $dir/data || exit 1
  # train
  #dir=$data_fmllr/train
  #steps/nnet/make_fmllr_feats.sh --nj 10 --cmd "$train_cmd" \
  #   --transform-dir ${gmmdir}_ali \
  #   $dir data/train $gmmdir $dir/log $dir/data || exit 1

  ## split the data : 90% train 10% cross-validation (held-out)
  #dir=$data/train
  #utils/subset_data_dir_tr_cv.sh $dir ${dir}_tr90 ${dir}_cv10 || exit 1
#fi

if [ $stage -le 1 ]; then
  # Pre-train DBN, i.e. a stack of RBMs (small database, smaller DNN)
  dir=$pretraindir
  if [ ! -f $dir/.done ]; then
    (tail --pid=$$ -F $dir/log/pretrain_dbn.log 2>/dev/null)& # forward log
    $cuda_cmd $dir/log/pretrain_dbn.log \
      steps/nnet/pretrain_dbn.sh \
        --hid-dim 1024 --rbm-iter 20 $data/train $dir || exit 1;
    touch $dir/.done
  fi
fi

if [ $stage -le 2 ]; then
  # Train the DNN optimizing per-frame cross-entropy.
  dir=${pretraindir}_dnn
  ali=${gmmdir}_ali
  feature_transform=$pretraindir/final.feature_transform
  dbn=$pretraindir/6.dbn
  (tail --pid=$$ -F $dir/log/train_nnet.log 2>/dev/null)& # forward log
  # Train
  if [ ! -f $dir/.done ]; then
    $cuda_cmd $dir/log/train_nnet.log \
      steps/nnet/train.sh --feature-transform $feature_transform --dbn $dbn --hid-layers 0 --learn-rate 0.008 \
      $data/train_tr90 $data/train_cv10 data/lang $ali $ali $dir || exit 1;
    touch $dir/.done
  fi
  # Decode (reuse HCLG graph)
  if [ ! -f $dir/decode_test/.done ]; then
    steps/nnet/decode.sh --nj 20 --cmd "$decode_cmd" --acwt 0.2 \
      $gmmdir/graph $data/test $dir/decode_test && touch $dir/decode_test/.done &
  fi
  if [ ! -f $dir/decode_dev/.done ]; then
    steps/nnet/decode.sh --nj 20 --cmd "$decode_cmd" --acwt 0.2 \
      $gmmdir/graph $data/dev $dir/decode_dev && touch $dir/decode_dev/.done &
  fi
  if [ ! -f $dir/decode_enh/.done ]; then
    steps/nnet/decode.sh --nj 20 --cmd "$decode_cmd" --acwt 0.2 \
      $gmmdir/graph $data/enh $dir/decode_enh && touch $dir/decode_enh/.done &
  fi
fi

if [ "$skip_smbr" == "false" ]; then
  # Sequence training using sMBR criterion, we do Stochastic-GD 
  # with per-utterance updates. We use usually good acwt 0.1
  dir=${pretraindir}_dnn_smbr
  srcdir=${pretraindir}_dnn
  acwt=0.2

  if [ $stage -le 3 ]; then
    # First we generate lattices and alignments:
    if [ ! -f ${srcdir}_ali/.done ]; then
      steps/nnet/align.sh --nj 20 --cmd "$train_cmd" \
	$data/train data/lang $srcdir ${srcdir}_ali || exit 1;
      touch ${srcdir}_ali/.done
    fi
    if [ ! -f ${srcdir}_denlats/.done ]; then
      steps/nnet/make_denlats.sh --nj 20 --cmd "$decode_cmd" --acwt $acwt --mono true \
	--lattice-beam 10.0 --beam 18.0 \
	$data/train data/lang $srcdir ${srcdir}_denlats || exit 1;
      touch ${srcdir}_denlats/.done
    fi
  fi

  if [ $stage -le 4 ]; then
    # Re-train the DNN by 6 iterations of sMBR 
    if [ ! -f $dir/.done ]; then
      steps/nnet/train_mpe.sh --cmd "$cuda_cmd" --num-iters 6 --acwt $acwt \
	--do-smbr true \
	$data/train data/lang $srcdir ${srcdir}_ali ${srcdir}_denlats $dir || exit 1
      touch $dir/.done
    fi
    # Decode
    for ITER in 1 6; do
      if [ ! -f $dir/decode_test_it${ITER}/.done ]; then
	steps/nnet/decode.sh --nj 20 --cmd "$decode_cmd" \
	  --nnet $dir/${ITER}.nnet --acwt $acwt \
	  $gmmdir/graph $data/test $dir/decode_test_it${ITER} && touch $dir/decode_test_it${ITER}/.done &
      fi
      if [ ! -f $dir/decode_dev_it${ITER}/.done ]; then
	steps/nnet/decode.sh --nj 20 --cmd "$decode_cmd" \
	  --nnet $dir/${ITER}.nnet --acwt $acwt \
	  $gmmdir/graph $data/dev $dir/decode_dev_it${ITER} && touch $dir/decode_dev_it${ITER}/.done &
      fi
      if [ ! -f $dir/decode_enh_it${ITER}/.done ]; then
	steps/nnet/decode.sh --nj 20 --cmd "$decode_cmd" \
	  --nnet $dir/${ITER}.nnet --acwt $acwt \
	  $gmmdir/graph $data/enh $dir/decode_enh_it${ITER} && touch $dir/decode_enh_it${ITER}/.done &
      fi
    done 
  fi
fi # if [ "$skip_smbr" == "false" ]

echo Success
exit 0

# Getting results [see RESULTS file]
# for x in exp/*/decode*; do [ -d $x ] && grep WER $x/wer_* | utils/best_wer.sh; done
