#include "QN_config.h"

#include "CRF.h"

#include "fst/fstlib.h"
#include "CRF_Model.h"
#include "io/CRF_FeatureStreamManager.h"
#include "trainers/CRF_SGTrainer.h"
#include "trainers/CRF_LBFGSTrainer.h"
#include "trainers/CRF_AISTrainer.h"
#include "ftrmaps/CRF_StdFeatureMap.h"
#include "ftrmaps/CRF_StdSparseFeatureMap.h"
//#include "CRF_StdRange.h"
//#include "CRF_StdTransRange.h"
//#include "CRF_StdTransFeatureMap.h"


using namespace std;

static struct {
	char* ftr1_file;
	char* ftr1_format;
	int ftr1_width;
	int ftr1_ftr_start;
	int ftr1_ftr_count;
	int ftr1_window_offset;
	int ftr1_window_len;
	int ftr1_delta_order;
	int ftr1_delta_win;
	char* ftr1_norm_file;
	char* ftr2_file;
	char* ftr2_format;
	int ftr2_width;
	int ftr2_ftr_start;
	int ftr2_ftr_count;
	int ftr2_window_offset;
	int ftr2_window_len;
	int ftr2_delta_order;
	int ftr2_delta_win;
	char* ftr2_norm_file;
	char* ftr3_file;
	char* ftr3_format;
	int ftr3_width;
	int ftr3_ftr_start;
	int ftr3_ftr_count;
	int ftr3_window_offset;
	int ftr3_window_len;
	int ftr3_delta_order;
	int ftr3_delta_win;
	char* ftr3_norm_file;
	int train_sent_start;
	int train_sent_count;
	int window_extent;
	char* train_sent_range;
	char* cv_sent_range;
	char* hardtarget_file;
	int hardtarget_window_offset;
	char* out_weight_file;
	char* init_weight_file;
	char* avg_weight_file;
	int avg_weight_present;
	int crf_label_size;
	int crf_bunch_size;
	int crf_epochs;
	int crf_utt_rpt;
	char* crf_train_method;
	char* crf_objective_function;
	char* crf_train_order;
	float crf_lr;
	int crf_random_seed;
	char* crf_featuremap;
	char* crf_featuremap_file;
	int crf_stateftr_start;
	int crf_stateftr_end;
	int crf_transftr_start;
	int crf_transftr_end;
	int crf_states;
	float crf_gauss_var;
	int crf_use_state_bias;
	int crf_use_trans_bias;
	float crf_state_bias_value;
	float crf_trans_bias_value;
	float crf_ais_l1alpha;
	int threads;
	int verbose;
	int dummy;
} config;

QN_ArgEntry argtab[] =
{
	{ NULL, "ASR CRaFT CRF training program version " CRF_VERSION, QN_ARG_DESC },
	{ "ftr1_file", "Input feature file", QN_ARG_STR, &(config.ftr1_file), QN_ARG_REQ },
	{ "ftr1_format", "Input feature file format [pfile]", QN_ARG_STR, &(config.ftr1_format) },
	{ "ftr1_width", "Input feature file columns", QN_ARG_INT, &(config.ftr1_width) },
	{ "ftr1_ftr_start", "First feature used", QN_ARG_INT, &(config.ftr1_ftr_start) },
	{ "ftr1_ftr_count", "Number of features used", QN_ARG_INT, &(config.ftr1_ftr_count) },
	{ "ftr1_window_offset", "Offset of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr1_window_offset) },
	{ "ftr1_window_len", "Length of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr1_window_len) },
	{ "ftr1_delta_order", "Delta order", QN_ARG_INT, &(config.ftr1_delta_order) },
	{ "ftr1_delta_win", "Delta window", QN_ARG_INT, &(config.ftr1_delta_win) },
	{ "ftr1_norm_file", "Normalization parameters for ftr1_file", QN_ARG_STR, &(config.ftr1_norm_file) },
	{ "ftr2_file", "Input feature file", QN_ARG_STR, &(config.ftr2_file) },
	{ "ftr2_format", "Input feature file format [pfile]", QN_ARG_STR, &(config.ftr2_format) },
	{ "ftr2_width", "Input feature file columns", QN_ARG_INT, &(config.ftr2_width) },
	{ "ftr2_ftr_start", "First feature used", QN_ARG_INT, &(config.ftr2_ftr_start) },
	{ "ftr2_ftr_count", "Number of features used", QN_ARG_INT, &(config.ftr2_ftr_count) },
	{ "ftr2_window_offset", "Offset of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr2_window_offset) },
	{ "ftr2_window_len", "Length of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr2_window_len) },
	{ "ftr2_delta_order", "Delta order", QN_ARG_INT, &(config.ftr2_delta_order) },
	{ "ftr2_delta_win", "Delta window", QN_ARG_INT, &(config.ftr2_delta_win) },
	{ "ftr2_norm_file", "Normalization parameters for ftr2_file", QN_ARG_STR, &(config.ftr2_norm_file) },
	{ "ftr3_file", "Input feature file", QN_ARG_STR, &(config.ftr3_file) },
	{ "ftr3_format", "Input feature file format [pfile]", QN_ARG_STR, &(config.ftr3_format) },
	{ "ftr3_width", "Input feature file columns", QN_ARG_INT, &(config.ftr3_width) },
	{ "ftr3_ftr_start", "First feature used", QN_ARG_INT, &(config.ftr3_ftr_start) },
	{ "ftr3_ftr_count", "Number of features used", QN_ARG_INT, &(config.ftr3_ftr_count) },
	{ "ftr3_window_offset", "Offset of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr3_window_offset) },
	{ "ftr3_window_len", "Length of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr3_window_len) },
	{ "ftr3_delta_order", "Delta order", QN_ARG_INT, &(config.ftr3_delta_order) },
	{ "ftr3_delta_win", "Delta window", QN_ARG_INT, &(config.ftr3_delta_win) },
	{ "ftr3_norm_file", "Normalization parameters for ftr3_file", QN_ARG_STR, &(config.ftr3_norm_file) },
	{ "window_extent", "Extent of all windows (frames)", QN_ARG_INT, &(config.window_extent) },
	{ "train_sent_range", "Training sentence indices in QN_Range(3) format", QN_ARG_STR, &(config.train_sent_range), QN_ARG_REQ },
	{ "cv_sent_range", "CV sentence indices in QN_Range(3) format", QN_ARG_STR, &(config.cv_sent_range) },
	{ "hardtarget_file", "Target Label File", QN_ARG_STR, &(config.hardtarget_file), QN_ARG_REQ },
	{ "hardtarget_window_offset", "Offset of hardtarget file (frames)", QN_ARG_INT, &(config.hardtarget_window_offset) },
	{ "out_weight_file", "Output Weight File", QN_ARG_STR, &(config.out_weight_file), QN_ARG_REQ },
	{ "init_weight_file", "Input Seed Weight File", QN_ARG_STR, &(config.init_weight_file) },
	{ "avg_weight_file", "Input Average Weight File", QN_ARG_STR, &(config.avg_weight_file) },
	{ "avg_weight_present", "No. of presentations used to build average weight file", QN_ARG_INT, &(config.avg_weight_present) },
	{ "crf_label_size", "Number of CRF output labels", QN_ARG_INT, &(config.crf_label_size), QN_ARG_REQ },
	{ "crf_bunch_size", "Bunch size for CRF processing", QN_ARG_INT, &(config.crf_bunch_size) },
	{ "crf_epochs", "Maximum number of epochs", QN_ARG_INT, &(config.crf_epochs) },
	{ "crf_utt_rpt", "Block size to report progress on", QN_ARG_INT, &(config.crf_utt_rpt) },
	{ "crf_lr", "Learning rate", QN_ARG_FLOAT, &(config.crf_lr) },
	{ "crf_random_seed", "Presentation order random seed", QN_ARG_INT, &(config.crf_random_seed) },
	{ "crf_train_method", "CRF training method (sg|lbfgs)", QN_ARG_STR, &(config.crf_train_method) },
	{ "crf_objective_function", "Objective function for gradient training (expf|softexpf|ferr)", QN_ARG_STR, &(config.crf_objective_function) },
	{ "crf_train_order", "Presentation order of samples for training (seq|random|noreplace)", QN_ARG_STR, &(config.crf_train_order) },
	{ "crf_featuremap", "Association of inputs to feature functions (stdstate|stdtrans|stdsparse|stdsparsetrans|file)", QN_ARG_STR, &(config.crf_featuremap) },
	{ "crf_featuremap_file", "File containing map of inputs to feature functions", QN_ARG_STR, &(config.crf_featuremap_file) },
	{ "crf_stateftr_start", "Feature index to start computing state feature funcs", QN_ARG_INT, &(config.crf_stateftr_start) },
	{ "crf_stateftr_end", "Feature index to end computing state feature funcs", QN_ARG_INT, &(config.crf_stateftr_end) },
	{ "crf_transftr_start", "Feature index to start computing trans feature funcs", QN_ARG_INT, &(config.crf_transftr_start) },
	{ "crf_transftr_end", "Feature index to end computing trans feature funcs", QN_ARG_INT, &(config.crf_transftr_end) },
	{ "crf_states", "Number of states per label", QN_ARG_INT, &(config.crf_states) },
	{ "crf_gauss_var", "Gaussian variance", QN_ARG_FLOAT, &(config.crf_gauss_var) },
	{ "crf_use_state_bias", "Use state bias functions", QN_ARG_BOOL, &(config.crf_use_state_bias) },
	{ "crf_use_trans_bias", "Use transition bias functions", QN_ARG_BOOL, &(config.crf_use_trans_bias) },
	{ "crf_state_bias_value", "Function value for state bias functions", QN_ARG_FLOAT, &(config.crf_state_bias_value) },
	{ "crf_trans_bias_value", "Function value for transition bias functions", QN_ARG_FLOAT, &(config.crf_trans_bias_value) },
	{ "threads", "Number of threads to use for multithreaded trainers", QN_ARG_INT, &(config.threads) },
	{ "crf_ais_l1alpha", "l1 alpha threshold for AIS training", QN_ARG_FLOAT, &(config.crf_ais_l1alpha) },
	//	{ "dummy", "Output status messages", QN_ARG_INT, &(config.dummy) },
	{ "verbose", "Output status messages", QN_ARG_INT, &(config.verbose) },
	{ NULL, NULL, QN_ARG_NOMOREARGS }
};

static struct CRF_FeatureMap_config fmap_config;

static void set_defaults(void) {
	config.ftr1_file="";
	config.ftr1_format="pfile";
	config.ftr1_width=0;
	config.ftr1_ftr_start=0;
	config.ftr1_ftr_count=0;
	config.ftr1_window_offset=0;
	config.ftr1_window_len=1;
	config.ftr1_delta_order=0;
	config.ftr1_delta_win=0;
	config.ftr1_norm_file=NULL;
	config.ftr2_file="";
	config.ftr2_format="pfile";
	config.ftr2_width=0;
	config.ftr2_ftr_start=0;
	config.ftr2_ftr_count=0;
	config.ftr2_window_offset=0;
	config.ftr2_window_len=1;
	config.ftr2_delta_order=0;
	config.ftr2_delta_win=0;
	config.ftr2_norm_file=NULL;
	config.ftr3_file="";
	config.ftr3_format="pfile";
	config.ftr3_width=0;
	config.ftr3_ftr_start=0;
	config.ftr3_ftr_count=0;
	config.ftr3_window_offset=0;
	config.ftr3_window_len=1;
	config.ftr3_delta_order=0;
	config.ftr3_delta_win=0;
	config.ftr3_norm_file=NULL;
	config.window_extent=1;
	config.train_sent_range="";
	config.cv_sent_range=0;
	config.hardtarget_file="";
	config.hardtarget_window_offset=0;
	config.out_weight_file="out.weights";
	config.init_weight_file=NULL;
	config.avg_weight_file=NULL;
	config.avg_weight_present=0;
	config.crf_bunch_size=1;
	config.crf_epochs=10;
	config.crf_utt_rpt=100;
	config.crf_lr=0.008;
	config.crf_random_seed=0;
	config.crf_train_method="sg";
	config.crf_objective_function="expf";
	config.crf_train_order="random";
	config.crf_featuremap="stdstate";
	config.crf_featuremap_file=NULL;
	config.crf_stateftr_start=0;
	config.crf_stateftr_end=-1;
	config.crf_transftr_start=0;
	config.crf_transftr_end=-1;
	config.crf_states=1;
	config.crf_gauss_var=0.0;
	config.crf_use_state_bias=1;
	config.crf_use_trans_bias=1;
	config.crf_state_bias_value=1.0;
	config.crf_trans_bias_value=1.0;
	config.crf_ais_l1alpha=0.0;
	config.threads=1;
	config.verbose=0;
};

static void set_fmap_config(QNUInt32 nfeas) {
	fmap_config.map_type=STDSTATE;
	if (strcmp(config.crf_featuremap,"stdtrans")==0) { fmap_config.map_type=STDTRANS;}
	if (strcmp(config.crf_featuremap,"stdsparse")==0) { fmap_config.map_type=STDSPARSE;}
	if (strcmp(config.crf_featuremap,"stdsparsetrans")==0) { fmap_config.map_type=STDSPARSETRANS;}
	if (strcmp(config.crf_featuremap,"file")==0) { fmap_config.map_type=INFILE;}
	fmap_config.numLabs=config.crf_label_size;
	fmap_config.numFeas=nfeas;
	fmap_config.numStates=config.crf_states;

	fmap_config.useStateFtrs=true;
	fmap_config.stateFidxStart=config.crf_stateftr_start;
	if (config.crf_stateftr_end>=0) {
		fmap_config.stateFidxEnd=config.crf_stateftr_end;
	}
	else {
		fmap_config.stateFidxEnd=nfeas-1;
	}
	if (fmap_config.map_type==STDTRANS || fmap_config.map_type==STDSPARSETRANS) {
		fmap_config.useTransFtrs=true;
		fmap_config.transFidxStart=config.crf_transftr_start;
		if (config.crf_transftr_end>=0) {
			fmap_config.transFidxEnd=config.crf_transftr_end;
		}
		else {
			fmap_config.transFidxEnd=nfeas-1;
		}
	}
	else {
		fmap_config.useTransFtrs=false;
	}
	if (config.crf_use_state_bias == 0) {
		fmap_config.useStateBias=false;
	}
	else {
		fmap_config.useStateBias=true;
	}
	if (config.crf_use_trans_bias == 0) {
		fmap_config.useTransBias=false;
	}
	else {
		fmap_config.useTransBias=true;
	}
	fmap_config.stateBiasVal=config.crf_state_bias_value;
	fmap_config.transBiasVal=config.crf_trans_bias_value;

};



int main(int argc, const char* argv[]) {
	char* progname;

	cout << "SIZEOF_SIZE_T: " << QN_SIZEOF_SIZE_T;
	cout << " QN_SIZET_BAD: " << QN_SIZET_BAD;
	cout << " QN_SIZET_MAX: " << QN_SIZET_MAX;
	cout << " QN_SEGID_BAD: " << QN_SEGID_BAD;
	cout << " QN_ALL: " << QN_ALL << endl;
	set_defaults();
	QN_initargs(&argtab[0], &argc, &argv, &progname);
	QN_printargs(NULL, progname, &argtab[0]);
	cout << "IN LABELS: " << config.crf_label_size << endl;

	seqtype trn_seq = RANDOM_REPLACE;
	ftrmaptype trn_ftrmap = STDSTATE;
	trntype trn_type = SGTRAIN;
	objfunctype ofunc_type = EXPF;

	if (strcmp(config.crf_train_order,"seq")==0) { trn_seq=SEQUENTIAL;}
	if (strcmp(config.crf_train_order,"noreplace")==0) { trn_seq=RANDOM_NO_REPLACE;}

	if (strcmp(config.crf_featuremap,"stdtrans")==0) { trn_ftrmap=STDTRANS;}
	if (strcmp(config.crf_featuremap,"stdsparse")==0) { trn_ftrmap=STDSPARSE;}
	if (strcmp(config.crf_featuremap,"stdsparsetrans")==0) { trn_ftrmap=STDSPARSETRANS;}
	if (strcmp(config.crf_featuremap,"file")==0) { trn_ftrmap=INFILE;}

	if (strcmp(config.crf_train_method,"sg")==0) {trn_type=SGTRAIN;}
	if (strcmp(config.crf_train_method,"lbfgs")==0) {trn_type=LBFGSTRAIN;}
	if (strcmp(config.crf_train_method,"ais")==0) {trn_type=AISTRAIN;}

	if (strcmp(config.crf_objective_function,"expf")==0) { ofunc_type=EXPF; }
	if (strcmp(config.crf_objective_function,"softexpf")==0) { ofunc_type=EXPFSOFT; }
	if (strcmp(config.crf_objective_function,"ferr")==0) { ofunc_type=FERR; }

	if ((trn_ftrmap == INFILE) && (config.crf_featuremap_file==NULL)) {
		cerr << "ERROR: crf_featuremap_file must be non-NULL when crf_featuremap set to 'file'" << endl;
		exit(-1);
	}
	CRF_FeatureStreamManager* str2=NULL;
	CRF_FeatureStreamManager* str3=NULL;

	CRF_FeatureStreamManager str1(1,"ftr1_file",config.ftr1_file,config.ftr1_format,config.hardtarget_file, config.hardtarget_window_offset,
							(size_t) config.ftr1_width, (size_t) config.ftr1_ftr_start, (size_t) config.ftr1_ftr_count,
							config.window_extent, config.ftr1_window_offset, config.ftr1_window_len,
							config.ftr1_delta_order, config.ftr1_delta_win,
							config.train_sent_range, config.cv_sent_range,
							NULL,0,0,0,trn_seq,config.crf_random_seed,config.threads);
	if (strcmp(config.ftr2_file,"") != 0) {

		/*CRF_FeatureStreamManager str2(1,"ftr2_file",config.ftr2_file,config.ftr2_format,config.hardtarget_file, config.hardtarget_window_offset,
							(size_t) config.ftr2_width, (size_t) config.ftr2_ftr_start, (size_t) config.ftr2_ftr_count,
							config.window_extent, config.ftr2_window_offset, config.ftr2_window_len,
							config.ftr2_delta_order, config.ftr2_delta_win,
							config.train_sent_range, config.cv_sent_range,
							NULL,0,0,0,trn_seq,config.crf_random_seed,config.threads);
		str1.join(&str2);*/
		str2=new CRF_FeatureStreamManager(1,"ftr2_file",config.ftr2_file,config.ftr2_format,config.hardtarget_file, config.hardtarget_window_offset,
				(size_t) config.ftr2_width, (size_t) config.ftr2_ftr_start, (size_t) config.ftr2_ftr_count,
				config.window_extent, config.ftr2_window_offset, config.ftr2_window_len,
				config.ftr2_delta_order, config.ftr2_delta_win,
				config.train_sent_range, config.cv_sent_range,
				NULL,0,0,0,trn_seq,config.crf_random_seed,config.threads);
		str1.join(str2);
	}
	if (strcmp(config.ftr3_file,"") != 0) {

		/*CRF_FeatureStreamManager str3(1,"ftr3_file",config.ftr3_file,config.ftr3_format,config.hardtarget_file, config.hardtarget_window_offset,
							(size_t) config.ftr3_width, (size_t) config.ftr3_ftr_start, (size_t) config.ftr3_ftr_count,
							config.window_extent, config.ftr3_window_offset, config.ftr3_window_len,
							config.ftr3_delta_order, config.ftr3_delta_win,
							config.train_sent_range, config.cv_sent_range,
							NULL,0,0,0,trn_seq,config.crf_random_seed,config.threads);
		str1.join(&str3);*/
		str3=new CRF_FeatureStreamManager(1,"ftr3_file",config.ftr3_file,config.ftr3_format,config.hardtarget_file, config.hardtarget_window_offset,
									(size_t) config.ftr3_width, (size_t) config.ftr3_ftr_start, (size_t) config.ftr3_ftr_count,
									config.window_extent, config.ftr3_window_offset, config.ftr3_window_len,
									config.ftr3_delta_order, config.ftr3_delta_win,
									config.train_sent_range, config.cv_sent_range,
									NULL,0,0,0,trn_seq,config.crf_random_seed,config.threads);
		str1.join(str3);
	}

	CRF_Model my_crf(config.crf_label_size);
	cout << "LABELS: " << my_crf.getNLabs() << endl;

	set_fmap_config(str1.getNumFtrs());
	my_crf.setFeatureMap(CRF_FeatureMap::createFeatureMap(&fmap_config));
	/*CRF_FeatureMap* my_map=NULL;
	if (trn_ftrmap == STDSPARSE || trn_ftrmap == STDSPARSETRANS) {
		CRF_StdSparseFeatureMap* tmp_map=new CRF_StdSparseFeatureMap(config.crf_label_size,str1.getNumFtrs());
		if (trn_ftrmap == STDSPARSE) {
			if (config.crf_stateftr_end>=0) {
				tmp_map->setStateFtrRange(config.crf_stateftr_start,config.crf_stateftr_end);
			}
			tmp_map->recalc();
		}
		else {
			if (config.crf_stateftr_end>=0) {
				tmp_map->setStateFtrRange(config.crf_stateftr_start,config.crf_stateftr_end);
			}
			tmp_map->setUseTransFtrs(true);
			if (config.crf_transftr_end>=0) {
				tmp_map->setTransFtrRange(config.crf_transftr_start,config.crf_transftr_end);
			}
			tmp_map->recalc();
		}
		tmp_map->setStateBiasVal(config.crf_state_bias_value);
		tmp_map->setTransBiasVal(config.crf_trans_bias_value);
		if (config.crf_use_state_bias != 1) {
			tmp_map->setUseStateBias(false);
		}
		if (config.crf_use_trans_bias != 1) {
			tmp_map->setUseTransBias(false);
		}
		my_map=tmp_map;
	}
	else if (trn_ftrmap == STDSTATE) {
		CRF_StdFeatureMap* tmp_map=new CRF_StdFeatureMap(config.crf_label_size,str1.getNumFtrs());
		if (config.crf_stateftr_end>=0) {
			tmp_map->setStateFtrRange(config.crf_stateftr_start,config.crf_stateftr_end);
		}
		tmp_map->setStateBiasVal(config.crf_state_bias_value);
		tmp_map->setTransBiasVal(config.crf_trans_bias_value);
		if (config.crf_use_state_bias != 1) {
			tmp_map->setUseStateBias(false);
		}
		if (config.crf_use_trans_bias != 1) {
			tmp_map->setUseTransBias(false);
		}
		tmp_map->recalc();
		my_map=tmp_map;
	}
	else if (trn_ftrmap == STDTRANS) {
		CRF_StdFeatureMap* tmp_map=new CRF_StdFeatureMap(config.crf_label_size,str1.getNumFtrs());
		if (config.crf_stateftr_end>=0) {
			tmp_map->setStateFtrRange(config.crf_stateftr_start,config.crf_stateftr_end);
		}
		tmp_map->setUseTransFtrs(true);
		if (config.crf_transftr_end>=0) {
			tmp_map->setTransFtrRange(config.crf_transftr_start,config.crf_transftr_end);
		}
		tmp_map->setStateBiasVal(config.crf_state_bias_value);
		tmp_map->setTransBiasVal(config.crf_trans_bias_value);
		if (config.crf_use_state_bias != 1) {
			tmp_map->setUseStateBias(false);
		}
		if (config.crf_use_trans_bias != 1) {
			tmp_map->setUseTransBias(false);
		}
		tmp_map->recalc();
		my_map=tmp_map;
	}
	else if (trn_ftrmap== INFILE) {
		cerr << "Reading featuremaps from file not yet implemented" << endl;
		exit(-1);
	}
	my_map->setNumStates(config.crf_states);
	my_map->recalc();
	my_crf.setFeatureMap(my_map);*/
	cout << "FEATURES: " << my_crf.getLambdaLen() << endl;
	cout << "LABELS: " << my_crf.getNLabs() << endl;
	if (config.init_weight_file != NULL) {
		bool chkfile=my_crf.readFromFile(config.init_weight_file);
		if (!chkfile) {
			cerr << "ERROR! File " << config.init_weight_file << " unable to be opened for reading.  ABORT!" << endl;
			exit(-1);
		}
		if ((config.avg_weight_present >0) && (config.avg_weight_file != NULL)) {
			bool chkfile=my_crf.readAverageFromFile(config.avg_weight_file,config.avg_weight_present);
			if (!chkfile) {
				cerr << "ERROR! File " << config.avg_weight_file << " unable to be opened for reading.  ABORT!" << endl;
				exit(-1);
			}
		}
		else {

		}
	}
	CRF_Trainer* my_trainer;
	switch (trn_type) {
	case LBFGSTRAIN :
		my_trainer = new CRF_LBFGSTrainer(&my_crf,&str1,config.out_weight_file);
		((CRF_LBFGSTrainer *)my_trainer)->setObjectiveFunction(ofunc_type);
		break;
	case AISTRAIN :
		my_trainer = new CRF_AISTrainer(&my_crf,&str1,config.out_weight_file);
		((CRF_AISTrainer *)my_trainer)->setl1alpha(config.crf_ais_l1alpha);
		break;
	case SGTRAIN :
	default:
		my_trainer = new CRF_SGTrainer(&my_crf,&str1,config.out_weight_file);
		((CRF_SGTrainer *)my_trainer)->setObjectiveFunction(ofunc_type);
		break;
	}

	my_trainer->setMaxIters(config.crf_epochs);
	my_trainer->setLR(config.crf_lr);
	my_trainer->setUttRpt(config.crf_utt_rpt);
	if (config.crf_gauss_var != 0.0) {
		my_trainer->setGaussVar(config.crf_gauss_var);
	}
	try {
		my_trainer->train();
	}
	catch (exception &e) {
		cerr << "Exception: " << e.what() << endl;
		exit(-1);
	}
}
