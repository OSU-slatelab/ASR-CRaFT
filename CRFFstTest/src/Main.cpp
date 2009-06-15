#include "QN_config.h"
#include "fst/lib/fstlib.h"

#include "CRF.h"
#include "CRF_FeatureStream.h"
#include "CRF_FeatureStreamManager.h"
#include "CRF_Model.h"
#include "CRF_StdFeatureMap.h"
#include "CRF_StdTransFeatureMap.h"
#include "CRF_LatticeBuilder.h"
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
	char* crf_lm;
	char* crf_lat_outdir;
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
	{ "crf_featuremap", "Association of inputs to feature functions (stdstate|stdtrans|file)", QN_ARG_STR, &(config.crf_featuremap) },
	{ "crf_featuremap_file", "File containing map of inputs to feature functions", QN_ARG_STR, &(config.crf_featuremap_file) },
	{ "crf_stateftr_start", "Feature index to start computing state feature funcs", QN_ARG_INT, &(config.crf_stateftr_start) },
	{ "crf_stateftr_end", "Feature index to end computing state feature funcs", QN_ARG_INT, &(config.crf_stateftr_end) },
	{ "crf_transftr_start", "Feature index to start computing trans feature funcs", QN_ARG_INT, &(config.crf_transftr_start) },
	{ "crf_transftr_end", "Feature index to end computing trans feature funcs", QN_ARG_INT, &(config.crf_transftr_end) },
	{ "crf_states", "Number of states per label", QN_ARG_INT, &(config.crf_states) },
	{ "crf_lm", "Language model file name (in OpenFST binary format", QN_ARG_STR, &(config.crf_lm) },
	{ "crf_lat_outdir", "Output directory for lattice files (in OpenFST binary format)", QN_ARG_STR, &(config.crf_lat_outdir) },
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
	config.crf_lm=NULL;
	config.crf_lat_outdir=NULL;
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
	if (strcmp(config.crf_featuremap,"file")==0) { dec_ftrmap=INFILE;}
	
	if ((dec_ftrmap == INFILE) && (config.crf_featuremap_file==NULL)) {
		cerr << "ERROR: crf_featuremap_file must be non-NULL when crf_featuremap set to 'file'" << endl;
		exit(-1);
	} 
    
    
    
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
		cout << "Joining files" << endl;
		str1.join(&str2);
		cout << "Files joined" << endl;
	}
    cout << "Feature File created" << endl;
	
	
	CRF_FeatureStream* crf_ftr_str = str1.trn_stream;
	CRF_Model my_crf(config.crf_label_size);
	cout << "LABELS: " << my_crf.getNLabs() << endl;
	
	CRF_FeatureMap* my_map;
	if (dec_ftrmap == STDSTATE) {
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
	else {
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
	StdFst* lm_fst = NULL;
	bool* test;
	test = new bool;
	*test=false;
	
	if (config.crf_lm != NULL) {
		cout << "Yo!" << endl;
		
		lm_fst=StdVectorFst::Read(config.crf_lm);
		if (lm_fst == NULL) {
			cerr << "ERROR: Failed opening file: " << config.crf_lm << endl;
			exit(-1);
		}
		cout << "Yo! Yo!" << endl;
	}
	//CRF_Viterbi* vit = new CRF_Viterbi(crf_ftr_str,&my_crf);
	//CRF_NewViterbi* vit = new CRF_NewViterbi(crf_ftr_str,&my_crf);
	CRF_LatticeBuilder lb(crf_ftr_str,&my_crf);
	crf_ftr_str->rewind();
	StdVectorFst* fst;
	QN_SegID segid = crf_ftr_str->nextseg();
	int count=0;
	typedef StdArc::StateId StateId;
	StdVectorFst final_result;
	FILE* outl=fopen(config.crf_output_labelfile,"w+");
	QN_OutLabStream* labout = new QN_OutLabStream_ILab(1, "out_labfile", outl, 255, 1);
	while (segid != QN_SEGID_BAD) {
		cout << "Processing segment " << count << endl;
		
		try {
			if (config.crf_states == 1) {
				fst=lb.buildLattice();
			}
			else {
				cout << "**Processing as nState Lattice" << endl; 
				fst=lb.nStateBuildLattice();
			}
			if (lm_fst != NULL) {
				cout << "**Applying language model fst" << endl;
				//result=new StdVectorFst();
				StdComposeFst* result=new StdComposeFst(*fst,*lm_fst);
				//Compose(*fst, *lm_fst, result);
				cout << "**Projecting language model fst" << endl;
				//Project(result,PROJECT_OUTPUT);
				cout << "**Removing epsilonds from fst" << endl;
				//RmEpsilon(result);
				//TopSort(result);
				delete fst;
				if (config.crf_lat_outdir != NULL) {
					cout << "**Writing Result lattice to " << string(config.crf_lat_outdir) << endl;
					string fst_fname=string(config.crf_lat_outdir)+"/fst."+stringify(count)+".fst";
					result->Write(fst_fname);
				}				
				cout << "**Finding Shortest Path " << endl;
				ShortestPath(*result,&final_result,1);
				Project(&final_result,PROJECT_OUTPUT);
				RmEpsilon(&final_result);
				TopSort(&final_result);
				delete result;
				
			}
			else {
				if (config.crf_lat_outdir != NULL) {
					cout << "**Writing Result lattice to " << string(config.crf_lat_outdir) << endl;
					string fst_fname=string(config.crf_lat_outdir)+"/fst."+stringify(count)+".fst";
					fst->Write(fst_fname);
				}				
				cout << "**Finding Shortest Path " << endl;
				ShortestPath(*fst,&final_result,1);
				RmEpsilon(&final_result);
				TopSort(&final_result);
			}
			// This portion constructs an ilab file from the above
			if (config.crf_lat_outdir != NULL ) {
				cout << "**Writing Shortest Path lattice to " << string(config.crf_lat_outdir) << endl;
				string fst_fname=string(config.crf_lat_outdir)+"/fst."+stringify(count)+".final.fst";
				final_result.Write(fst_fname);
			}
			//Get total number of labels - with a single path this is number of states - 1
			if (config.crf_output_labelfile != NULL) {
				QNUInt32 num_labels = final_result.NumStates() - 1;
				QNUInt32* labarr = new QNUInt32[num_labels];
				//Gets the initial state; if kNoState => empty FST.
				StateId initial_state = final_result.Start();
				//	Iterates over the FSTs states.
				QNUInt32 frm_no=0;
				cout << "**Writing Shortest Path Sequence to ilabfile " << endl;
				for (StateIterator<StdFst> siter(final_result); !siter.Done(); siter.Next())
				{ 
  					StateId state_id = siter.Value();
  					//Iterates over state state_id's arcs.
  					for (ArcIterator<StdFst> aiter(final_result, state_id); !aiter.Done(); aiter.Next())
					{
  						const StdArc &arc = aiter.Value();
  						//cout << count << "\t" << frm_no << "\t" << arc.olabel-1 << endl;
  						labarr[frm_no]=arc.olabel-1;
  						frm_no++;
					}
				}
				labout->write_labs(num_labels,labarr);
				labout->doneseg(segid);
				delete labarr;
			}
			//delete result;
		}
		catch (exception &e) {
			cerr << "Exception: " << e.what() << endl;
		}
		segid=crf_ftr_str->nextseg();
		count++;
	}
	//final_result.Write("binary.fst");
	delete labout; // explicitly delete the labelstream to flush contents to disk.
	fclose(outl);
}
