SUBDIRS = kaldi-mods
bin_SCRIPTS = \
	scrf-timit-demo.sh
nobase_dist_pkglibexec_SCRIPTS = \
	scrf-scripts/prep_exp.sh \
	scrf-scripts/sortResult.kaldi.sh \
	util/ilab_PhnDur2Dur.pl \
	util/ilab_PhnDur2Phn.pl \
	util/ilab2ctm.pl \
	util/map_ctm_3state_1state.pl
nobase_pkglibexec_SCRIPTS = \
	scrf-scripts/run_crf_fstdecode.sh \
	scrf-scripts/run_crf_train.sh \
	util/timit_score.sh 
nobase_dist_pkgdata_DATA = \
	timit-aux/phn-dur-lab-map-timit.txt \
	timit-aux/symlist-timit.txt \
	timit-aux/timit_sisx_test.core.neworder.olist \
	timit-aux/timit_sisx_test.core+rest.neworder.olist \
	timit-aux/timit_sisx_test.dt_set.neworder.olist \
	timit-aux/timit_sisx_train.olist \
	timit-aux/timit_test39.mlf \
	timit-aux/timit_train.48labs.ilab \
	timit-aux/timit_train.mlf
pkgdata_DATA = segmental-timit-demo.cfg
EXTRA_DIST = \
	scrf-scripts/run_crf_fstdecode.sh.in \
	scrf-scripts/run_crf_train.sh.in \
	util/timit_score.sh.in \
	segmental-timit-demo.cfg.in \
	scrf-timit-demo.sh.in
CLEANFILES = $(nobase_pkglibexec_SCRIPTS) $(pkgdata_DATA) $(bin_SCRIPTS)

subs= sed -e 's,[@]timitdir[@],$(TIMITDIR),g' \
	-e 's,[@]pkglibexecdir[@],$(pkglibexecdir),g' \
	-e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

scrf-scripts/run_crf_fstdecode.sh : scrf-scripts/run_crf_fstdecode.sh.in
	mkdir -p ./scrf-scripts
	$(subs) < $(srcdir)/scrf-scripts/run_crf_fstdecode.sh.in > scrf-scripts/run_crf_fstdecode.sh

scrf-scripts/run_crf_train.sh : scrf-scripts/run_crf_train.sh.in
	mkdir -p ./scrf-scripts
	$(subs) < $(srcdir)/scrf-scripts/run_crf_train.sh.in > scrf-scripts/run_crf_train.sh

util/timit_score.sh : util/timit_score.sh.in
	mkdir -p ./util
	$(subs) < $(srcdir)/util/timit_score.sh.in > util/timit_score.sh

segmental-timit-demo.cfg : segmental-timit-demo.cfg.in
	$(subs) < $(srcdir)/segmental-timit-demo.cfg.in > segmental-timit-demo.cfg

scrf-timit-demo.sh : scrf-timit-demo.sh.in
	$(subs) < $(srcdir)/scrf-timit-demo.sh.in > scrf-timit-demo.sh

*.sh : $@
	chmod +x $@

*.pl : $@
	chmod +x $@
