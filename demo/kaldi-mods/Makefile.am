timitdir = $(TIMITDIR)
nobase_dist_timit_DATA = \
	conf/all_spk.list \
	conf/dev_spk.list \
	conf/enh_spk.list \
	conf/test_spk.list \
	mylocal/timit_sisx_test.neworder.olist.kaldi_to_htk_to_pfile \
	mylocal/timit_sisx_train.olist.kaldi_to_htk_to_pfile 
nobase_dist_timit_SCRIPTS = \
	run.sh \
	mylocal/htk_to_pf.sh \
	mylocal/prep_bn_feats.sh \
	mylocal/run_dnn_fmllr_mono.sh \
	mylocal/timit_data_prep.sh \
	mylocal/timit_format_data.sh \
	mysteps/decode_fmllr.sh \
	mysteps/make_bn_feats.sh

*.sh : $@
	chmod +x $@
