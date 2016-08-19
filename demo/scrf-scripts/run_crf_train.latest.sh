#!/bin/bash

input_arg=`echo $0 $@`

. ./cmd.sh  # for $train_cmd

. /u/drspeech/share/lib/icsiargs.sh #TODO: change to getopt
ASRCRAFTBASE=/u/stiff/ASR-CRaFT-built/build
#PATH=$ASRCRAFTBBASE/v0.01h/bin/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01i/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01i/paramTying/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01j_NewDurRep/paramTying/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01k_segmental/paramTying/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01l_segmental_WithoutDurLab/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01m_all_segmental_model/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01n_all_segmental_model_with_boundary_feature/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01l_all_segmental_model_with_boundary_feature_fast_model/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01l_all_segmental_model_with_boundary_feature_fast_model/fast_model_no_mem_leak_with_init_iter/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01q_all_segmental_model_with_boundary_feature_fast_model_with_boundary_delta_ftr/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01r_1state_segmental_all_one_file/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01s_1state_segmental_all_one_file_fixed_multistate_weight_indexing/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01t_nstate_segmental_all_one_file/:$ASRCRAFTBASE/helper/:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01o_individual_segmental_features/sample-avg-max-min-dur_debugNStateSegDecodeFst:$ASRCRAFTBASE/helper/:/data/data2/hey/SegmentalCRF/misc/:$ICSIPATH:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01v_nstate_segmental_correct_maxv_all_one_file:$ASRCRAFTBASE/helper/:/data/data2/hey/SegmentalCRF/misc/:$ICSIPATH:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01w_nstate_segmental_done_file:$ASRCRAFTBASE/helper/:/data/data2/hey/SegmentalCRF/misc/:$ICSIPATH:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01.latest:$ASRCRAFTBASE/helper/:/data/data2/hey/SegmentalCRF/misc/:$ICSIPATH:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01x_nstate_segmental_parellel.detail2:$ASRCRAFTBASE/helper/:/data/data2/hey/SegmentalCRF/misc/:$ICSIPATH:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01x_nstate_segmental_parellel.detail4.fixmem:$ASRCRAFTBASE/helper/:/data/data2/hey/SegmentalCRF/misc/:$ICSIPATH:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01x_nstate_segmental_parellel.detail5.adagrad_resume:$ASRCRAFTBASE/helper/:/data/data2/hey/SegmentalCRF/misc/:$ICSIPATH:$PATH; export PATH
#PATH=/data/data2/hey/SegmentalCRF/bin/ASR-CRaFT_v0.01x_nstate_segmental_parellel.detail6.avg_minibatch.non_verbose:$ASRCRAFTBASE/helper/:/data/data2/hey/SegmentalCRF/misc/:$ICSIPATH:$PATH; export PATH
PATH=$ASRCRAFTBASE/CRFTrain:$PATH; export PATH
LD_LIBRARY_PATH=/u/drspeech/opt/OpenFst-1.3.4/lib/:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH

USAGE="Usage: 

$0 
  argfile=[arg file, default=./conf]
  lr=[learning rate]
  lr_decay=[learning rate decaying factor, default=1.0]
  eta=[AdaGrad scaling factor, default=1.0. Only one of eta and lr can be set, based on whether to use AdaGrad.]
  wt_pre=[weight dir prefix, default=./weights]
  SEGMODEL=[stdframe(default)|stdseg|stdseg_no_dur|stdseg_no_dur_no_transftr|stdseg_no_dur_no_segtransftr]
  autoresume=[1(default)|0]
  mb=[minibatch size, default is 1]
  nj=[number of jobs, default is 1]
  mem=[memory requirement, default=nj*938 (MB), optional]
  nodes=[machines to allocate the jobs, e.g. feldspar, optional]
  init_weight=[init weight file, optional]
  avg_weight=[avg_weight_file, optional]
  grad_sqr_acc=[grad_sqr_acc_file, optional]
  present=[avg_weight_presentation, optional]
  start_iter=[starting training iteration, optional]
"

if [ ! -z "$1" ] && [ "$1" == "-h" ]; then
  echo "$USAGE"
  exit 1
fi

if [ -z "$argfile" ]; then
  argfile=conf
  echo "The arg file argument is not set, using default: $argfile"
fi

if [ ! -f "$argfile" ]; then
  echo "Error: cannot access the argfile: $argfile" 1>&2
  exit 1
fi

. $argfile

if [ ! -z "$lr" ] && [ -z "$eta" ]; then
  grad_opt="crf_lr=$lr"
elif [ -z "$lr" ] && [ ! -z "$eta" ]; then
  grad_opt="crf_use_adagrad=true crf_adagrad_eta=$eta"
else
  echo "Either lr or eta has to be set, but cannot be both set."
  echo "$USAGE"
  exit 1
fi

if [ -z "$lr_decay" ]; then
  lr_decay=1.0
  echo "lr_decay is set to default: 1.0"
fi

if [ -z "$wt_pre" ]; then
  wt_pre=./weights
  echo "No output weight directory found, using default: $wt_pre"
fi

if [ -z "$autoresume" ]; then
  autoresume=1
  echo "autoresume is set to default: 1"
fi

if [ -z "$mb" ]; then
  mb=1
  echo "mb is set to default: 1"
fi

if [ -z "$nj" ]; then
  nj=1
  echo "nj is set to default: 1"
fi

if [ -z "$mem" ]; then
  mem=$((nj * 938)) # MB
  echo "mem is set to default: $mem"
fi

[ ! -z "$nodes" ] && nodes_opt="-w $nodes" || nodes_opt=

if [ ! -z "$lr" ]; then
  WBASE=$wt_pre.lr$lr
else
  WBASE=$wt_pre.ag_eta$eta
fi
echo "Using output weight directory $WBASE"

if [ ! -d "$WBASE" ]; then
  mkdir $WBASE
fi

WF=$WBASE/$WFPREFIX
LOGF=$WBASE/crf.sgtrain.log

# resume from previous stopped iteration
TrainingResumeParam=
if [ -z $init_weight ] && [ -z $avg_weight ]; then
  if [ -z $start_iter ] && [ $autoresume -eq 1 ]; then
    #start_iter=`/data/data2/hey/SegmentalCRF/misc/getStartIterForTrain.sh WBASE=$WBASE WFPREFIX=$WFPREFIX`
    iter=0
    while [ 1 ]; do
      if [ ! -f $WBASE/.done.train.i$iter ]; then
        start_iter=$iter
        break
      fi
      iter=$((iter + 1))
    done
  fi
  echo "start_iter=$start_iter"
  if [ ! -z $start_iter ]; then
    if [ $start_iter -lt 0 ]; then
      echo "start_iter should be greater than or equal to 0." 1>&2
      exit 1
    elif [ $start_iter -gt 0 ]; then
      if [ $END -ge $START ]; then
        utts_per_iter=$((END - START + 1))
        seed_iter=$((start_iter - 1))
        init_weight=${WF}.i${seed_iter}.out
        avg_weight=${WF}.i${seed_iter}.avg.out
        present=$((start_iter * utts_per_iter))

        if [ -z $grad_sqr_acc ]; then
          tmp=${WF}.i${seed_iter}.gradSqrAcc.out
          if [ -f $tmp ]; then
            grad_sqr_acc=$tmp
          fi
        fi

        #TrainingResumeParam=`/data/data2/hey/SegmentalCRF/misc/getTrainingResumeParam.pl -w $WF -s $start_iter -n $utts_per_iter`
        TrainingResumeParam="init_weight_file=$init_weight avg_weight_file=$avg_weight avg_weight_present=$present grad_sqr_acc_file=$grad_sqr_acc init_iter=$start_iter"
      fi
    fi
  fi
else
  TrainingResumeParam="init_weight_file=$init_weight avg_weight_file=$avg_weight avg_weight_present=$present grad_sqr_acc_file=$grad_sqr_acc init_iter=$start_iter"
fi
TRAIN_ARGS=$TRAIN_ARGS" $TrainingResumeParam"

echo "TRAIN_ARGS=$TRAIN_ARGS"

which CRFTrain | tee -a $LOGF

echo `date` 2>&1 | tee -a $LOGF
echo 2>&1 | tee -a $LOGF
echo $input_arg 2>&1 | tee -a $LOGF
echo 2>&1 | tee -a $LOGF
echo "CRFTrain
  out_weight_file=$WF
  crf_label_size=$NUML
  ftr1_file=$FEAF1
  ftr2_file=$FEAF2
  hardtarget_file=$LABF
  train_sent_range=$START-$END
  cv_sent_range=0
  $grad_opt
  crf_lr_decay_rate=$lr_decay
  crf_utt_rpt=$CHECKPOINT
  crf_states=$NSTATES
  crf_epochs=$ITERS
  use_broken_class_label=$USEBROKENCLASSLABEL
  crf_bunch_size=$mb
  threads=$nj
  $TRAIN_ARGS" 2>&1 | tee -a $LOGF
#  ftr3_file=$FEAF3

#$train_cmd -pe smp $nj -l mem_free=${nj}G $LOGF \
srun -c $nj --mem=$mem $nodes_opt \
CRFTrain \
  out_weight_file=$WF \
  crf_label_size=$NUML \
  ftr1_file=$FEAF1 \
  ftr2_file=$FEAF2 \
  hardtarget_file=$LABF \
  train_sent_range=$START-$END \
  cv_sent_range=0 \
  $grad_opt \
  crf_lr_decay_rate=$lr_decay \
  crf_utt_rpt=$CHECKPOINT \
  crf_states=$NSTATES \
  crf_epochs=$ITERS \
  use_broken_class_label=$USEBROKENCLASSLABEL \
  crf_bunch_size=$mb \
  threads=$nj \
  $TRAIN_ARGS 2>&1 | tee -a $LOGF
#  removed for demo
#  ftr3_file=$FEAF3 \

echo `date` 2>&1 | tee -a $LOGF
