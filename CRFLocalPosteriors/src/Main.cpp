#include "QN_config.h"

#include "CRF.h"
#include "CRF_Model.h"
#include "CRF_FeatureStream.h"
#include "CRF_FeatureStreamManager.h"
#include "CRF_LocalPosteriorBuilder.h"
#include "CRF_NewLocalPosteriorBuilder.h"
#include "CRF_StdFeatureMap.h"
#include "CRF_StdTransFeatureMap.h"
#include "CRF_StdSparseFeatureMap.h"
#include "CRF_LatticeBuilder.h"

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
	char* ftr2_file;
	char* ftr2_format;
	int ftr2_width;
	int ftr2_ftr_start;
	int ftr2_ftr_count;
	int ftr2_window_offset;
	int ftr2_window_len;
	int ftr2_delta_order;
	int ftr2_delta_win;
	char* ftr3_file;
	char* ftr3_format;
	int ftr3_width;
	int ftr3_ftr_start;
	int ftr3_ftr_count;
	int ftr3_window_offset;
	int ftr3_window_len;
	int ftr3_delta_order;
	int ftr3_delta_win;
	int train_sent_start;
	int train_sent_count;
	int window_extent;
	char* hardtarget_file;
	int hardtarget_window_offset;
	char* weight_file;
	int crf_label_size;
	int crf_bunch_size;
	int crf_epochs;
	char* crf_output_ftrfile;
	char* crf_eval_range;
	char* crf_ftr_type;
	char* crf_featuremap;
	char* crf_featuremap_file;
	int crf_stateftr_start;
	int crf_stateftr_end;
	int crf_transftr_start;
	int crf_transftr_end;
	int crf_states;
	int crf_use_state_bias;
	int crf_use_trans_bias;
	float crf_state_bias_value;
	float crf_trans_bias_value;
	int verbose;
	int dummy;
} config;

QN_ArgEntry argtab[] =
{
	{ NULL, "ASR CRaFT CRF local posterior program version " CRF_VERSION, QN_ARG_DESC },
	{ "ftr1_file", "Input feature file", QN_ARG_STR, &(config.ftr1_file), QN_ARG_REQ },
	{ "ftr1_format", "Input feature file format [pfile]", QN_ARG_STR, &(config.ftr1_format) },
	{ "ftr1_width", "Input feature file columns", QN_ARG_INT, &(config.ftr1_width) },
	{ "ftr1_ftr_start", "First feature used", QN_ARG_INT, &(config.ftr1_ftr_start) },
	{ "ftr1_ftr_count", "Number of features used", QN_ARG_INT, &(config.ftr1_ftr_count) },
	{ "ftr1_window_offset", "Offset of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr1_window_offset) },
	{ "ftr1_window_len", "Length of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr1_window_len) },
	{ "ftr1_delta_order", "Delta order", QN_ARG_INT, &(config.ftr1_delta_order) },
	{ "ftr1_delta_win", "Delta window", QN_ARG_INT, &(config.ftr1_delta_win) },
	{ "ftr2_file", "Input feature file", QN_ARG_STR, &(config.ftr2_file) },
	{ "ftr2_format", "Input feature file format [pfile]", QN_ARG_STR, &(config.ftr2_format) },
	{ "ftr2_width", "Input feature file columns", QN_ARG_INT, &(config.ftr2_width) },
	{ "ftr2_ftr_start", "First feature used", QN_ARG_INT, &(config.ftr2_ftr_start) },
	{ "ftr2_ftr_count", "Number of features used", QN_ARG_INT, &(config.ftr2_ftr_count) },
	{ "ftr2_window_offset", "Offset of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr2_window_offset) },
	{ "ftr2_window_len", "Length of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr2_window_len) },
	{ "ftr2_delta_order", "Delta order", QN_ARG_INT, &(config.ftr2_delta_order) },
	{ "ftr2_delta_win", "Delta window", QN_ARG_INT, &(config.ftr2_delta_win) },
	{ "ftr3_file", "Input feature file", QN_ARG_STR, &(config.ftr3_file) },
	{ "ftr3_format", "Input feature file format [pfile]", QN_ARG_STR, &(config.ftr3_format) },
	{ "ftr3_width", "Input feature file columns", QN_ARG_INT, &(config.ftr3_width) },
	{ "ftr3_ftr_start", "First feature used", QN_ARG_INT, &(config.ftr3_ftr_start) },
	{ "ftr3_ftr_count", "Number of features used", QN_ARG_INT, &(config.ftr3_ftr_count) },
	{ "ftr3_window_offset", "Offset of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr3_window_offset) },
	{ "ftr3_window_len", "Length of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr3_window_len) },
	{ "ftr3_delta_order", "Delta order", QN_ARG_INT, &(config.ftr3_delta_order) },
	{ "ftr3_delta_win", "Delta window", QN_ARG_INT, &(config.ftr3_delta_win) },
	{ "window_extent", "Extent of all windows (frames)", QN_ARG_INT, &(config.window_extent) },
	{ "hardtarget_file", "Target Label File", QN_ARG_STR, &(config.hardtarget_file) },
	{ "hardtarget_window_offset", "Offset of hardtarget file (frames)", QN_ARG_INT, &(config.hardtarget_window_offset) },
	{ "weight_file", "Input Weight File", QN_ARG_STR, &(config.weight_file), QN_ARG_REQ },
	{ "crf_label_size", "Number of CRF output labels", QN_ARG_INT, &(config.crf_label_size), QN_ARG_REQ },
	{ "crf_bunch_size", "Bunch size for CRF processing", QN_ARG_INT, &(config.crf_bunch_size) },
	{ "crf_epochs", "Maximum number of epochs", QN_ARG_INT, &(config.crf_epochs) },
	{ "crf_output_ftrfile", "Output feature file name", QN_ARG_STR, &(config.crf_output_ftrfile), QN_ARG_REQ },
	{ "crf_eval_range", "Range of utterances to evaluate", QN_ARG_STR, &(config.crf_eval_range), QN_ARG_REQ },
	{ "crf_ftr_type", "Type of feature to output (log|posterior|linear)",QN_ARG_STR, &(config.crf_ftr_type)},
	{ "crf_featuremap", "Association of inputs to feature functions (stdstate|stdtrans|stdparse|stdparsetrans|file)", QN_ARG_STR, &(config.crf_featuremap) },
	{ "crf_featuremap_file", "File containing map of inputs to feature functions", QN_ARG_STR, &(config.crf_featuremap_file) },
	{ "crf_stateftr_start", "Feature index to start computing state feature funcs", QN_ARG_INT, &(config.crf_stateftr_start) },
	{ "crf_stateftr_end", "Feature index to end computing state feature funcs", QN_ARG_INT, &(config.crf_stateftr_end) },
	{ "crf_transftr_start", "Feature index to start computing trans feature funcs", QN_ARG_INT, &(config.crf_transftr_start) },
	{ "crf_transftr_end", "Feature index to end computing trans feature funcs", QN_ARG_INT, &(config.crf_transftr_end) },
	{ "crf_states", "Number of states per label", QN_ARG_INT, &(config.crf_states) },
	{ "crf_use_state_bias", "Use state bias functions", QN_ARG_BOOL, &(config.crf_use_state_bias) },
	{ "crf_use_trans_bias", "Use transition bias functions", QN_ARG_BOOL, &(config.crf_use_trans_bias) },
	{ "crf_state_bias_value", "Function value for state bias functions", QN_ARG_FLOAT, &(config.crf_state_bias_value) },
	{ "crf_trans_bias_value", "Function value for transition bias functions", QN_ARG_FLOAT, &(config.crf_trans_bias_value) },
	//{ "dummy", "Output status messages", QN_ARG_INT, &(config.dummy) },
	{ "verbose", "Output status messages", QN_ARG_INT, &(config.verbose) }
};

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
	config.ftr2_file="";
	config.ftr2_format="pfile";
	config.ftr2_width=0;
	config.ftr2_ftr_start=0;
	config.ftr2_ftr_count=0;
	config.ftr2_window_offset=0;
	config.ftr2_window_len=1;
	config.ftr2_delta_order=0;
	config.ftr2_delta_win=0;
	config.ftr3_file="";
	config.ftr3_format="pfile";
	config.ftr3_width=0;
	config.ftr3_ftr_start=0;
	config.ftr3_ftr_count=0;
	config.ftr3_window_offset=0;
	config.ftr3_window_len=1;
	config.ftr3_delta_order=0;
	config.ftr3_delta_win=0;
	config.window_extent=1;
	config.hardtarget_file=NULL;
	config.hardtarget_window_offset=0;
	config.crf_bunch_size=1;
	config.crf_epochs=10;
	config.crf_output_ftrfile="";
	config.crf_eval_range=NULL;
	config.crf_ftr_type="log";
	config.crf_featuremap="stdstate";
	config.crf_featuremap_file=NULL;
	config.crf_stateftr_start=0;
	config.crf_stateftr_end=-1;
	config.crf_transftr_start=0;
	config.crf_transftr_end=-1;
	config.crf_states=1;
	config.crf_use_state_bias=1;
	config.crf_use_trans_bias=1;
	config.crf_state_bias_value=1.0;
	config.crf_trans_bias_value=1.0;
	config.verbose=0;
};


int main(int argc, const char* argv[]) {
	char* progname;

	bool logftrs=true;
	bool normfeas=true;
	ftrmaptype trn_ftrmap = STDSTATE;

	set_defaults();
	QN_initargs(&argtab[0], &argc, &argv, &progname);
	QN_printargs(NULL, progname, &argtab[0]);
	cout << "IN LABELS: " << config.crf_label_size << endl;

    if (strcmp(config.crf_ftr_type,"posterior")==0) {
    	logftrs=false;
    }
    if (strcmp(config.crf_ftr_type,"linear")==0) {
    	logftrs=true;
    	normfeas=false;
    }

	if (strcmp(config.crf_featuremap,"stdtrans")==0) { trn_ftrmap=STDTRANS;}
	if (strcmp(config.crf_featuremap,"stdsparse")==0) { trn_ftrmap=STDSPARSE;}
	if (strcmp(config.crf_featuremap,"stdsparsetrans")==0) { trn_ftrmap=STDSPARSETRANS;}
	if (strcmp(config.crf_featuremap,"file")==0) { trn_ftrmap=INFILE;}


    CRF_FeatureStreamManager str1(1,"ftr1_file",config.ftr1_file,config.ftr1_format,config.hardtarget_file,config.hardtarget_window_offset,
							(size_t) config.ftr1_width, (size_t) config.ftr1_ftr_start, (size_t) config.ftr1_ftr_count,
							config.window_extent, config.ftr1_window_offset, config.ftr1_window_len,
							config.ftr1_delta_order, config.ftr1_delta_win,
							config.crf_eval_range, 0,
							NULL,0,0,0,SEQUENTIAL);

	if (strcmp(config.ftr2_file,"") != 0) {

		CRF_FeatureStreamManager str2(1,"ftr2_file",config.ftr2_file,config.ftr2_format,config.hardtarget_file,config.hardtarget_window_offset,
							(size_t) config.ftr2_width, (size_t) config.ftr2_ftr_start, (size_t) config.ftr2_ftr_count,
							config.window_extent, config.ftr2_window_offset, config.ftr2_window_len,
							config.ftr2_delta_order, config.ftr2_delta_win,
							config.crf_eval_range, 0,
							NULL,0,0,0,SEQUENTIAL);
		str1.join(&str2);
	}
    cout << "Feature File created" << endl;

	if (strcmp(config.ftr3_file,"") != 0) {

		CRF_FeatureStreamManager str3(1,"ftr3_file",config.ftr3_file,config.ftr3_format,config.hardtarget_file,config.hardtarget_window_offset,
							(size_t) config.ftr3_width, (size_t) config.ftr3_ftr_start, (size_t) config.ftr3_ftr_count,
							config.window_extent, config.ftr3_window_offset, config.ftr3_window_len,
							config.ftr3_delta_order, config.ftr3_delta_win,
							config.crf_eval_range, 0,
							NULL,0,0,0,SEQUENTIAL);
		cout << "Joining files" << endl;
		str1.join(&str3);
		cout << "Files joined" << endl;
	}
    cout << "Feature File created" << endl;

	CRF_FeatureStream* crf_ftr_str = str1.trn_stream;
	CRF_Model my_crf(config.crf_label_size);
	cout << "LABELS: " << my_crf.getNLabs() << endl;
	CRF_FeatureMap* my_map;
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
	my_crf.setFeatureMap(my_map);
	bool openchk=my_crf.readFromFile(config.weight_file);
	if (!openchk) {
		cerr << "ERROR: Failed opening file: " << config.weight_file << endl;
		exit(-1);
	}

#ifndef USE_NEW_POSTERIOR_CODE

	//CRF_LocalPosteriorBuilder* lpb = new CRF_LocalPosteriorBuilder(&my_crf);
	CRF_NewLocalPosteriorBuilder* lpb = new CRF_NewLocalPosteriorBuilder(&my_crf,normfeas);

	crf_ftr_str->rewind();
	FILE* outl=fopen(config.crf_output_ftrfile,"w+");
	QNUInt32 len=config.crf_label_size;
	QN_OutFtrStream* ftrout = new QN_OutFtrLabStream_PFile(1, "localposterior", outl, len, 0, 1);
	QN_SegID segid = crf_ftr_str->nextseg();
	int count=0;
	CRF_Seq* seq_head;
	CRF_StateVector* posteriorList;
	float ab_f[len];
	while (segid != QN_SEGID_BAD) {
		cout << "Processing segment " << count;
		try {
			//seq_head=lpb->buildFtrSeq(crf_ftr_str);
			posteriorList=lpb->buildFtrSeq(crf_ftr_str);
		}
		catch (exception &e) {
			cerr << "Exception: " << e.what() << endl;
			exit(-1);
		}
		cout << " ... Segment processed" << endl;
		CRF_Seq* cur_seq=seq_head;
		//while (cur_seq != NULL) {
		//for (QNUInt32 j=0; j<posteriorList->size(); j++) {
		for (QNUInt32 j=0; j<posteriorList->getNodeCount(); j++) {
			//double* ab = cur_seq->getAlphaBeta();
			double* ab = posteriorList->at(j)->getAlphaBeta();
			for (QNUInt32 i=0; i<len; i++) {
				if (!logftrs) {
					try {
						ab[i]=expE(ab[i]);
					}
					catch (exception& e) {
						cerr << "Exception: " << e.what() << endl;
						exit(-1);
					}
				}
				ab_f[i]=ab[i]; // Convert array from double to float for QuickNet interface
			}
			//cout << "Starting Write" << endl;
			ftrout->write_ftrs(1,ab_f);
			//cout << "Ending Write" << endl << "Starting Delete" << endl;
			//CRF_Seq* last_seq=cur_seq;
			//cur_seq=cur_seq->getNext();
			//delete last_seq;
			//cout << "Ending Delete" << endl;
		}
		ftrout->doneseg(segid);
		segid=crf_ftr_str->nextseg();
		count++;
		//delete posteriorList;
	}
	delete ftrout; // explicitly delete the labelstream to flush contents to disk.
	fclose(outl);

#else // USE_NEW_POSTERIOR_CALC
	//CRF_LocalPosteriorBuilder* lpb = new CRF_LocalPosteriorBuilder(&my_crf);
	//CRF_NewLocalPosteriorBuilder* lpb = new CRF_NewLocalPosteriorBuilder(&my_crf,normfeas);
	cout << "hello1" << endl;
	CRF_LatticeBuilder *lb=new CRF_LatticeBuilder(crf_ftr_str,&my_crf);
	cout << "hello2" << endl;
	crf_ftr_str->rewind();
	cout << "hello3 " << config.crf_output_ftrfile << endl;
	FILE* outl=fopen(config.crf_output_ftrfile,"w+");
	QNUInt32 len=config.crf_label_size;
	QN_OutFtrStream* ftrout = new QN_OutFtrLabStream_PFile(1, "localposterior", outl, len, 0, 1);
	QN_SegID segid = crf_ftr_str->nextseg();
	int count=0;
	cout << "hello4" << endl;
	vector<double> denominatorGamma;
	float ab_f[len];
	while (segid != QN_SEGID_BAD) {
		cout << "Processing segment " << count;
		try {
			//seq_head=lpb->buildFtrSeq(crf_ftr_str);
			//posteriorList=lpb->buildFtrSeq(crf_ftr_str);
			lb->getAlignmentGammas(&denominatorGamma,NULL,NULL,NULL);
		}
		catch (exception &e) {
			cerr << "Exception: " << e.what() << endl;
			exit(-1);
		}
		cout << " ... Segment processed" << endl;
		//CRF_Seq* cur_seq=seq_head;
		//while (cur_seq != NULL) {
		//for (QNUInt32 j=0; j<posteriorList->size(); j++) {
		for (QNUInt32 j=0; j<denominatorGamma.size(); j++) {
			//double* ab = cur_seq->getAlphaBeta();
			//double* ab = posteriorList->at(j)->getAlphaBeta();
			for (QNUInt32 i=0; i<len; i++) {
				if (!logftrs) {
					try {
						ab_f[i]=(float)(expE(denominatorGamma.at(j*len+i)));
					}
					catch (exception& e) {
						cerr << "Exception: " << e.what() << endl;
						exit(-1);
					}
				}
			}
			ftrout->write_ftrs(1,ab_f);
		}
		ftrout->doneseg(segid);
		segid=crf_ftr_str->nextseg();
		count++;
		//delete posteriorList;
		denominatorGamma.clear();
	}
	delete ftrout; // explicitly delete the labelstream to flush contents to disk.
	fclose(outl);

#endif
}
