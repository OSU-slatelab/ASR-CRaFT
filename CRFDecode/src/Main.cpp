#include "QN_config.h"
#include "CRF.h"

#include "CRF_FeatureStream.h"
#include "CRF_FeatureStreamManager.h"
#include "CRF_Model.h"
//#include "CRF_StdRange.h"
//#include "CRF_StdTransRange.h"
#include "CRF_StdFeatureMap.h"
#include "CRF_StdSparseFeatureMap.h"
#include "CRF_StdTransFeatureMap.h"
//#include "CRF_Viterbi.h"
#include "CRF_NewViterbi.h"
#include "CRF_LabelPath.h"

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
	char* crf_output_labelfile;
	char* crf_eval_range;
	char* crf_decode_mode;
	char* crf_featuremap;
	char* crf_featuremap_file;
	int crf_stateftr_start;
	int crf_stateftr_end;
	int crf_transftr_start;
	int crf_transftr_end;
	int crf_states;
	int verbose;
	int dummy;
} config;

QN_ArgEntry argtab[] =
{
	{ NULL, "ASR CRaFT CRF decoding program version " CRF_VERSION, QN_ARG_DESC },
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
	{ "crf_output_labelfile", "Output label file name", QN_ARG_STR, &(config.crf_output_labelfile), QN_ARG_REQ },
	{ "crf_eval_range", "Range of utterances to evaluate", QN_ARG_STR, &(config.crf_eval_range), QN_ARG_REQ },
	{ "crf_decode_mode", "Mode of decode (decode|align) ", QN_ARG_STR, &(config.crf_decode_mode) },
	{ "crf_featuremap", "Association of inputs to feature functions (stdstate|stdtrans|stdsparse|stdsparsetrans|file)", QN_ARG_STR, &(config.crf_featuremap) },
	{ "crf_featuremap_file", "File containing map of inputs to feature functions", QN_ARG_STR, &(config.crf_featuremap_file) },
	{ "crf_stateftr_start", "Feature index to start computing state feature funcs", QN_ARG_INT, &(config.crf_stateftr_start) },
	{ "crf_stateftr_end", "Feature index to end computing state feature funcs", QN_ARG_INT, &(config.crf_stateftr_end) },
	{ "crf_transftr_start", "Feature index to start computing trans feature funcs", QN_ARG_INT, &(config.crf_transftr_start) },
	{ "crf_transftr_end", "Feature index to end computing trans feature funcs", QN_ARG_INT, &(config.crf_transftr_end) },
	{ "crf_states", "Number of states per label", QN_ARG_INT, &(config.crf_states) },
	{ "dummy", "Output status messages", QN_ARG_INT, &(config.dummy) },
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
	config.hardtarget_file="";
	config.hardtarget_window_offset=0;
	config.crf_bunch_size=1;
	config.crf_epochs=10;
	config.crf_output_labelfile="";
	config.crf_eval_range=NULL;
	config.crf_decode_mode="decode";
	config.crf_featuremap="stdstate";
	config.crf_featuremap_file=NULL;
	config.crf_stateftr_start=0;
	config.crf_stateftr_end=-1;
	config.crf_transftr_start=0;
	config.crf_transftr_end=-1;
	config.crf_states=1;
	config.verbose=0;
};
 

int main(int argc, const char* argv[]) {
	char* progname;
	
	bool alignMode=false;
	
	set_defaults();
	QN_initargs(&argtab[0], &argc, &argv, &progname);
	QN_printargs(NULL, progname, &argtab[0]);
	cout << "IN LABELS: " << config.crf_label_size << endl;
	
	if (strcmp(config.crf_decode_mode,"align")==0) { alignMode=true; }
	    
    if (strcmp(config.hardtarget_file,"")==0) { config.hardtarget_file=NULL; }
    
    if ((config.hardtarget_file==NULL) && (alignMode)) {
    	cerr << "hardtarget_file required when crf_decode_mode=align" << endl;
    	exit(-1);
    }
    
    ftrmaptype dec_ftrmap = STDSTATE;

	if (strcmp(config.crf_featuremap,"stdtrans")==0) { dec_ftrmap=STDTRANS;}
	if (strcmp(config.crf_featuremap,"stdsparse")==0) { dec_ftrmap=STDSPARSE;}
	if (strcmp(config.crf_featuremap,"stdsparsetrans")==0) { dec_ftrmap=STDSPARSETRANS;}
		
	if ((dec_ftrmap == INFILE) && (config.crf_featuremap_file==NULL)) {
		cerr << "ERROR: crf_featuremap_file must be non-NULL when crf_featuremap set to 'file'" << endl;
		exit(-1);
	} 
    
    
    
    CRF_FeatureStreamManager str1(1,"ftr1_file",config.ftr1_file,config.ftr1_format,config.hardtarget_file, config.hardtarget_window_offset,
							(size_t) config.ftr1_width, (size_t) config.ftr1_ftr_start, (size_t) config.ftr1_ftr_count,
							config.window_extent, config.ftr1_window_offset, config.ftr1_window_len,
							config.ftr1_delta_order, config.ftr1_delta_win,
							config.crf_eval_range, 0,
							NULL,0,0,0,SEQUENTIAL);
	if (strcmp(config.ftr2_file,"") != 0) {
	
		CRF_FeatureStreamManager str2(1,"ftr2_file",config.ftr2_file,config.ftr2_format,config.hardtarget_file, config.hardtarget_window_offset,
							(size_t) config.ftr2_width, (size_t) config.ftr2_ftr_start, (size_t) config.ftr2_ftr_count,
							config.window_extent, config.ftr2_window_offset, config.ftr2_window_len,
							config.ftr2_delta_order, config.ftr2_delta_win,
							config.crf_eval_range, 0,
							NULL,0,0,0,SEQUENTIAL);
		cout << "Joining files" << endl;
		str1.join(&str2);
		cout << "Files joined" << endl;
	}
	if (strcmp(config.ftr3_file,"") != 0) {
	
		CRF_FeatureStreamManager str3(1,"ftr3_file",config.ftr3_file,config.ftr3_format,config.hardtarget_file, config.hardtarget_window_offset,
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
	if (dec_ftrmap == STDSPARSE || dec_ftrmap == STDSPARSETRANS) {
		CRF_StdSparseFeatureMap* tmp_map=new CRF_StdSparseFeatureMap(config.crf_label_size,str1.getNumFtrs());
		if (dec_ftrmap == STDSPARSE) {
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
		my_map=tmp_map;
	}
	else if (dec_ftrmap == STDSTATE) {
		CRF_StdFeatureMap* tmp_map=new CRF_StdFeatureMap(config.crf_label_size,str1.getNumFtrs());
		if (config.crf_stateftr_end>=0) {
			tmp_map->setStateFtrRange(config.crf_stateftr_start,config.crf_stateftr_end);
		}
		tmp_map->recalc();
		my_map=tmp_map;
	}
	else if (dec_ftrmap == STDTRANS) {
		CRF_StdFeatureMap* tmp_map=new CRF_StdFeatureMap(config.crf_label_size,str1.getNumFtrs());
		if (config.crf_stateftr_end>=0) {
			tmp_map->setStateFtrRange(config.crf_stateftr_start,config.crf_stateftr_end);
		}
		tmp_map->setUseTransFtrs(true);
		if (config.crf_transftr_end>=0) {
			tmp_map->setTransFtrRange(config.crf_transftr_start,config.crf_transftr_end);
		}
		tmp_map->recalc();
		my_map=tmp_map;
	}
	else if (dec_ftrmap== INFILE) {
		cerr << "Reading featuremaps from file not yet implemented" << endl;
		exit(-1);
	}
	my_map->setNumStates(config.crf_states);
	my_crf.setFeatureMap(my_map);
	bool openchk=my_crf.readFromFile(config.weight_file);
	if (!openchk) {
		cerr << "ERROR: Failed opening file: " << config.weight_file << endl;
		exit(-1);
	}
	//CRF_Viterbi* vit = new CRF_Viterbi(crf_ftr_str,&my_crf);
	CRF_NewViterbi* vit = new CRF_NewViterbi(crf_ftr_str,&my_crf);
	crf_ftr_str->rewind();
	FILE* outl=fopen(config.crf_output_labelfile,"w+");
	QN_OutLabStream* labout = new QN_OutLabStream_ILab(1, "out_labfile", outl, 255, 1);
	QN_SegID segid = crf_ftr_str->nextseg();
	int count=0;
	CRF_LabelPath* best;
	while (segid != QN_SEGID_BAD) {
		cout << "Processing segment " << count << endl;
		try {
			if (alignMode) {
				best=vit->alignPath();
			}
			else {
				best=vit->bestPathVec();
			}
		}
		catch (exception &e) {
			cerr << "Exception: " << e.what() << endl;
		}
		size_t* bp=best->getLabelPath();
		size_t len=best->getLength();
		labout->write_labs(len,bp);
		labout->doneseg(segid);
		segid=crf_ftr_str->nextseg();
		count++;
		delete best;
	}
	delete labout; // explicitly delete the labelstream to flush contents to disk.
	fclose(outl);
}
