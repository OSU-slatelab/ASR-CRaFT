/*
 * CRFDecode.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Command line interface for CRF Viterbi decoding using a language model
 *  in OpenFst format.
 * Follows command line interface model for ICSI Quicknet.
 * Requires the OpenFst finite state library.
 */
#include "quicknet3/QN_config.h"
#include "quicknet3/QN_Range.h"
#include "fst/fstlib.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include "CRF.h"
#include "CRF_Model.h"
#include "io/CRF_FeatureStream.h"
#include "io/CRF_FeatureStreamManager.h"
#include "ftrmaps/CRF_StdFeatureMap.h"
#include "ftrmaps/CRF_StdSparseFeatureMap.h"
#include "io/CRF_MLFManager.h"
#include "decoders/CRF_ViterbiDecoder.h"
#include "decoders/CRF_ViterbiNode_PruneTrans.h"
#include "decoders/CRF_ViterbiDecoder_StdSeg_NoSegTransFtr.h"

using namespace std;
typedef StdArc::StateId StateId;
typedef StdArc::Weight Weight;

static struct CRF_FeatureMap_config fmap_config;

/*
 * command line options
 */
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

	// Added by Ryan, for context features and segmental models
	int ftr1_left_context_len;
	int ftr1_right_context_len;
	bool ftr1_extract_seg_ftr;
	bool ftr1_use_boundary_delta_ftr;

	char* ftr2_file;
	char* ftr2_format;
	int ftr2_width;
	int ftr2_ftr_start;
	int ftr2_ftr_count;
	int ftr2_window_offset;
	int ftr2_window_len;
	int ftr2_delta_order;
	int ftr2_delta_win;

	// Added by Ryan, for context features and segmental models
	int ftr2_left_context_len;
	int ftr2_right_context_len;
	bool ftr2_extract_seg_ftr;
	bool ftr2_use_boundary_delta_ftr;

	char* ftr3_file;
	char* ftr3_format;
	int ftr3_width;
	int ftr3_ftr_start;
	int ftr3_ftr_count;
	int ftr3_window_offset;
	int ftr3_window_len;
	int ftr3_delta_order;
	int ftr3_delta_win;

	// Added by Ryan, for context features and segmental models
	int ftr3_left_context_len;
	int ftr3_right_context_len;
	bool ftr3_extract_seg_ftr;
	bool ftr3_use_boundary_delta_ftr;

	int train_sent_start;
	int train_sent_count;
	int window_extent;
	char* hardtarget_file;
	int hardtarget_window_offset;
	char* weight_file;
	int crf_label_size;
	int crf_bunch_size;
	int crf_epochs;
	char* crf_output_format;
	char* crf_output_labelfile;
	char* crf_output_mlffile;
	int crf_mlf_output_frames;
	int crf_mlf_output_states;
	char* crf_align_mlffile;
	char* crf_eval_range;
	char* crf_decode_mode;
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
	char* crf_lm_bin;

	// Added by Ryan
	char* crf_conf_bin;
	char* crf_disambig;

	char* crf_lm_arpa;
	char* crf_isymbols;
	char* crf_osymbols;
	char* crf_olist;
	char* crf_lat_outdir;

	// Added by Ryan
	bool crf_if_output_full_lat;
	char* htk_lat_outdir;

	float crf_decode_beam;
	int crf_decode_max_hyp;
	int crf_decode_min_hyp;
	float crf_decode_hyp_inc;
	int verbose;
	int dummy;

	//Added by Ryan, for segmental CRFs
	int label_maximum_duration;
	int num_actual_labs;
	char* crf_model_type;

} config;

/*
 * Command line options to be presented to the screen
 */
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

	// Added by Ryan, for context features and segmental models
	{ "ftr1_left_context_len", "Length of context features to the left of the window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr1_left_context_len) },
	{ "ftr1_right_context_len", "Length of context features to the right of the window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr1_right_context_len) },
	{ "ftr1_extract_seg_ftr", "Extract segment-level features (as opposed to frame-level features)", QN_ARG_BOOL, &(config.ftr1_extract_seg_ftr) },
	{ "ftr1_use_boundary_delta_ftr", "Use boundary delta features", QN_ARG_BOOL, &(config.ftr1_use_boundary_delta_ftr) },

	{ "ftr2_file", "Input feature file", QN_ARG_STR, &(config.ftr2_file) },
	{ "ftr2_format", "Input feature file format [pfile]", QN_ARG_STR, &(config.ftr2_format) },
	{ "ftr2_width", "Input feature file columns", QN_ARG_INT, &(config.ftr2_width) },
	{ "ftr2_ftr_start", "First feature used", QN_ARG_INT, &(config.ftr2_ftr_start) },
	{ "ftr2_ftr_count", "Number of features used", QN_ARG_INT, &(config.ftr2_ftr_count) },
	{ "ftr2_window_offset", "Offset of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr2_window_offset) },
	{ "ftr2_window_len", "Length of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr2_window_len) },
	{ "ftr2_delta_order", "Delta order", QN_ARG_INT, &(config.ftr2_delta_order) },
	{ "ftr2_delta_win", "Delta window", QN_ARG_INT, &(config.ftr2_delta_win) },

	// Added by Ryan, for context features and segmental models
	{ "ftr2_left_context_len", "Length of context features to the left of the window on ftr2_file (frames)", QN_ARG_INT, &(config.ftr2_left_context_len) },
	{ "ftr2_right_context_len", "Length of context features to the right of the window on ftr2_file (frames)", QN_ARG_INT, &(config.ftr2_right_context_len) },
	{ "ftr2_extract_seg_ftr", "Extract segment-level features (as opposed to frame-level features)", QN_ARG_BOOL, &(config.ftr2_extract_seg_ftr) },
	{ "ftr2_use_boundary_delta_ftr", "Use boundary delta features", QN_ARG_BOOL, &(config.ftr2_use_boundary_delta_ftr) },

	{ "ftr3_file", "Input feature file", QN_ARG_STR, &(config.ftr3_file) },
	{ "ftr3_format", "Input feature file format [pfile]", QN_ARG_STR, &(config.ftr3_format) },
	{ "ftr3_width", "Input feature file columns", QN_ARG_INT, &(config.ftr3_width) },
	{ "ftr3_ftr_start", "First feature used", QN_ARG_INT, &(config.ftr3_ftr_start) },
	{ "ftr3_ftr_count", "Number of features used", QN_ARG_INT, &(config.ftr3_ftr_count) },
	{ "ftr3_window_offset", "Offset of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr3_window_offset) },
	{ "ftr3_window_len", "Length of window on ftr1_file (frames)", QN_ARG_INT, &(config.ftr3_window_len) },
	{ "ftr3_delta_order", "Delta order", QN_ARG_INT, &(config.ftr3_delta_order) },
	{ "ftr3_delta_win", "Delta window", QN_ARG_INT, &(config.ftr3_delta_win) },

	// Added by Ryan, for context features and segmental models
	{ "ftr3_left_context_len", "Length of context features to the left of the window on ftr3_file (frames)", QN_ARG_INT, &(config.ftr3_left_context_len) },
	{ "ftr3_right_context_len", "Length of context features to the right of the window on ftr3_file (frames)", QN_ARG_INT, &(config.ftr3_right_context_len) },
	{ "ftr3_extract_seg_ftr", "Extract segment-level features (as opposed to frame-level features)", QN_ARG_BOOL, &(config.ftr3_extract_seg_ftr) },
	{ "ftr3_use_boundary_delta_ftr", "Use boundary delta features", QN_ARG_BOOL, &(config.ftr3_use_boundary_delta_ftr) },

	{ "window_extent", "Extent of all windows (frames)", QN_ARG_INT, &(config.window_extent) },
	{ "hardtarget_file", "Target Label File", QN_ARG_STR, &(config.hardtarget_file) },
	{ "hardtarget_window_offset", "Offset of hardtarget file (frames)", QN_ARG_INT, &(config.hardtarget_window_offset) },
	{ "weight_file", "Input Weight File", QN_ARG_STR, &(config.weight_file), QN_ARG_REQ },
	{ "crf_label_size", "Number of CRF output labels", QN_ARG_INT, &(config.crf_label_size), QN_ARG_REQ },
	{ "crf_bunch_size", "Bunch size for CRF processing", QN_ARG_INT, &(config.crf_bunch_size) },
	{ "crf_epochs", "Maximum number of epochs", QN_ARG_INT, &(config.crf_epochs) },
	{ "crf_output_format", "Output format type (ilab|mlf)", QN_ARG_STR, &(config.crf_output_format) },
	{ "crf_output_labelfile", "Output label file name", QN_ARG_STR, &(config.crf_output_labelfile) },
	{ "crf_output_mlffile", "Output MLF file name", QN_ARG_STR, &(config.crf_output_mlffile) },
	{ "crf_align_mlffile", "Target MLF filename for word alignment", QN_ARG_STR, &(config.crf_align_mlffile) },
	{ "crf_mlf_output_frames", "Output frame timings to MLF", QN_ARG_BOOL, &(config.crf_mlf_output_frames) },
	{ "crf_mlf_output_states", "Output state labels to MLF", QN_ARG_BOOL, &(config.crf_mlf_output_states) },
	{ "crf_eval_range", "Range of utterances to evaluate", QN_ARG_STR, &(config.crf_eval_range), QN_ARG_REQ },
	{ "crf_decode_mode", "Mode of decode (decode|align) ", QN_ARG_STR, &(config.crf_decode_mode) },
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
	{ "crf_lm_bin", "Language model file name (in OpenFST binary format)", QN_ARG_STR, &(config.crf_lm_bin) },
	{ "crf_lm_arpa", "Language model file name (in ARPA text format)", QN_ARG_STR, &(config.crf_lm_arpa) },

	// Added by Ryan
	{ "crf_conf_bin", "Confusion model file name (in OpenFST binary format)", QN_ARG_STR, &(config.crf_conf_bin) },
	{ "crf_disambig", "Disambiguation symbol IDs file name (in OpenFST binary format)", QN_ARG_STR, &(config.crf_disambig) },

	{ "crf_isymbols", "Input symbols file name (in OpenFST format)", QN_ARG_STR, &(config.crf_isymbols) },
	{ "crf_osymbols", "Output symbols file name (in OpenFST format)", QN_ARG_STR, &(config.crf_osymbols) },
    { "crf_olist", "Ordered list of output labels (for MLF)", QN_ARG_STR, &(config.crf_olist) },
	{ "crf_lat_outdir", "Output directory for lattice files (in OpenFST binary format)", QN_ARG_STR, &(config.crf_lat_outdir) },

	// Added by Ryan
	{ "crf_if_output_full_lat", "If output the full lattice (comparing to only a few best paths)", QN_ARG_BOOL, &(config.crf_if_output_full_lat) },
	{ "htk_lat_outdir", "Output directory for lattice files in HTK SLF format", QN_ARG_STR, &(config.htk_lat_outdir) },

	{ "crf_decode_beam", "Beam width for pruning", QN_ARG_FLOAT, &(config.crf_decode_beam) },
	{ "crf_decode_max_hyp", "Maximum hypotheses to keep in beam", QN_ARG_INT, &(config.crf_decode_max_hyp) },
	{ "crf_decode_min_hyp", "Minimum hypotheses to keep in beam", QN_ARG_INT, &(config.crf_decode_min_hyp) },
	{ "crf_decode_hyp_inc", "Increment for hypothesis beam pruning", QN_ARG_FLOAT, &(config.crf_decode_hyp_inc) },
	{ "verbose", "Output status messages", QN_ARG_INT, &(config.verbose) },

	//Added by Ryan, for segmental CRFs
	{ "label_maximum_duration", "The maximum duration if labels are phone-duration combination", QN_ARG_INT, &(config.label_maximum_duration) },
	{ "num_actual_labs", "The number of actual labels (without duration)", QN_ARG_INT, &(config.num_actual_labs) },
	{ "crf_model_type", "CRF model structure (stdframe|stdseg|stdseg_no_dur|stdseg_no_dur_no_transftr|stdseg_no_dur_no_segtransftr)", QN_ARG_STR, &(config.crf_model_type) },

	{ NULL, NULL, QN_ARG_NOMOREARGS }
};

/*
 * Default values for command line options
 */
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

	// Added by Ryan, for context features and segmental models
	config.ftr1_left_context_len=0;
	config.ftr1_right_context_len=0;
	config.ftr1_extract_seg_ftr=false;
	config.ftr1_use_boundary_delta_ftr=false;

	config.ftr2_file="";
	config.ftr2_format="pfile";
	config.ftr2_width=0;
	config.ftr2_ftr_start=0;
	config.ftr2_ftr_count=0;
	config.ftr2_window_offset=0;
	config.ftr2_window_len=1;
	config.ftr2_delta_order=0;
	config.ftr2_delta_win=0;

	// Added by Ryan, for context features and segmental models
	config.ftr2_left_context_len=0;
	config.ftr2_right_context_len=0;
	config.ftr2_extract_seg_ftr=false;
	config.ftr2_use_boundary_delta_ftr=false;

	config.ftr3_file="";
	config.ftr3_format="pfile";
	config.ftr3_width=0;
	config.ftr3_ftr_start=0;
	config.ftr3_ftr_count=0;
	config.ftr3_window_offset=0;
	config.ftr3_window_len=1;
	config.ftr3_delta_order=0;
	config.ftr3_delta_win=0;

	// Added by Ryan, for context features and segmental models
	config.ftr3_left_context_len=0;
	config.ftr3_right_context_len=0;
	config.ftr3_extract_seg_ftr=false;
	config.ftr3_use_boundary_delta_ftr=false;

	config.window_extent=1;
	config.hardtarget_file="";
	config.hardtarget_window_offset=0;
	config.crf_bunch_size=1;
	config.crf_epochs=10;
	config.crf_output_format="ilab";
	config.crf_output_labelfile=NULL;
	config.crf_output_mlffile=NULL;
	config.crf_align_mlffile=NULL;
	config.crf_mlf_output_frames=0;
	config.crf_mlf_output_states=0;
	config.crf_eval_range=NULL;
	config.crf_decode_mode="decode";
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
	config.crf_lm_bin=NULL;
	config.crf_lm_arpa=NULL;

	// Added by Ryan
	config.crf_conf_bin=NULL;
	config.crf_disambig=NULL;

	config.crf_isymbols=NULL;
	config.crf_osymbols=NULL;
	config.crf_olist=NULL;
	config.crf_lat_outdir=NULL;

	// Added by Ryan
	config.crf_if_output_full_lat=false;
	config.htk_lat_outdir=NULL;

	config.crf_decode_beam=0.0;
	config.crf_decode_max_hyp=0;
	config.crf_decode_min_hyp=0;
	config.crf_decode_hyp_inc=0.05;
	config.verbose=0;

	//Added by Ryan, for segmental CRFs
	config.label_maximum_duration=1;
	config.num_actual_labs=0;
	config.crf_model_type="stdframe";
};

/*
 * Sets initial values for the feature map object based on command line config object
 */
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

	//Added by Ryan, for segmental CRFs
	fmap_config.maxDur=config.label_maximum_duration;
	fmap_config.nActualLabs=config.num_actual_labs;
};

/*
 * logs an error message to standard output
 */
void log_msg(string outstr) {
	if (config.verbose) {
		cout << outstr << endl;
	}
}

/*
 * logs an error message to standard error
 */
void log_errmsg(string errstr) {
	if (config.verbose) {
		cerr << errstr << endl;
	}
}


// The following FST2HTK_lat related stuff are added by Ryan
/*
 * The class for storing the start and the end of an word arc with accumulated weights.
 * If it stores for a partial arc, just save the start node and leave the end node as -1.
 */
class HTK_word_arc {
public:
	int start_htk_node;
	int end_htk_node;
	double ac_weight;
	double lm_weight;

	HTK_word_arc(int s, double ac_w, double lm_w) : start_htk_node(s), end_htk_node(-1), ac_weight(ac_w), lm_weight(lm_w) {}
	HTK_word_arc(int s, int e, double ac_w, double lm_w) : start_htk_node(s), end_htk_node(e), ac_weight(ac_w), lm_weight(lm_w) {}
	HTK_word_arc(const HTK_word_arc& a) {
		start_htk_node = a.start_htk_node;
		end_htk_node = a.end_htk_node;
		ac_weight = a.ac_weight;
		lm_weight = a.lm_weight;
	}
	void add_weight(double ac, double lm) {
		ac_weight += ac;
		lm_weight += lm;
	}
};

/*
 * The FST lattice node that are constructed along the lattice expansion.
 * The in_paths collect all (partial) word arcs that go into this node.
 *
 */
class FST_lat_node {
public:
	StateId fst_node;
	int word_lab;
	int time_frame;
	int mapped_htk_node;
	vector<HTK_word_arc*> incoming_word_arcs;
	bool finish_word_arcs;  // true if all the word arcs that end at this node have already been added to the htk lattice

	FST_lat_node(StateId s, int l = kNoLabel, int t = BAD_TIME_FRAME) :
		fst_node(s), word_lab(l), time_frame(t), mapped_htk_node(kNoStateId), finish_word_arcs(false) {}

	~FST_lat_node() {
		for (vector<HTK_word_arc*>::iterator it = incoming_word_arcs.begin(); it != incoming_word_arcs.end(); ++it) {
			delete *it;
		}
	}

	void add_word_beginning_arc(vector<FST_lat_node*>* htk_nodes, FST_lat_node* prev_node, double arc_ac_weight, double arc_lm_weight) {
		// create an htk node if needed when adding an fst word beginning arc
		if (prev_node->mapped_htk_node == kNoStateId) {
			prev_node->mapped_htk_node = htk_nodes->size();
			htk_nodes->push_back(prev_node);
			//log_msg("FST_lat_node::add_word_beginning_arc() LOG: htk_nodes->push_back("
			//		+ stringify(prev_node->fst_node) + ") = " + stringify(prev_node->mapped_htk_node));
		}
		// add a partial htk word arc
		HTK_word_arc* htk_word_arc = new HTK_word_arc(prev_node->mapped_htk_node, arc_ac_weight, arc_lm_weight);
		incoming_word_arcs.push_back(htk_word_arc);
	}

	void add_word_internal_arc(const FST_lat_node& prev_node, double arc_ac_weight, double arc_lm_weight) {
		// no need to create an htk node for adding an fst word internal arc
		// accumulate all partial word arcs from previous node and extend them by adding the current weights
		for (uint i = 0; i < prev_node.incoming_word_arcs.size(); i++) {
			HTK_word_arc* htk_word_arc = new HTK_word_arc(*(prev_node.incoming_word_arcs[i]));  // copy from previous node
			htk_word_arc->add_weight(arc_ac_weight, arc_lm_weight);  // add weights
			incoming_word_arcs.push_back(htk_word_arc);  // add to current node
		}
	}

	void add_word_arcs_to_htk_lat(vector<FST_lat_node*>* htk_nodes, vector<HTK_word_arc*>* htk_arcs) {
		if (!finish_word_arcs) {
			if (mapped_htk_node == kNoStateId) {
				mapped_htk_node = htk_nodes->size();
				htk_nodes->push_back(this);
				//log_msg("FST_lat_node::add_word_arcs_to_htk_lat() LOG: htk_nodes->push_back("
				//		+ stringify(fst_node) + ") = " + stringify(mapped_htk_node));
			}
			for (uint i = 0; i < incoming_word_arcs.size(); ++i) {
				HTK_word_arc* htk_word_arc = new HTK_word_arc(*incoming_word_arcs[i]);
				htk_word_arc->end_htk_node = mapped_htk_node;
				htk_arcs->push_back(htk_word_arc);
				//log_msg("FST_lat_node::add_word_arcs_to_htk_lat() LOG: htk_arcs->push_back("
				//		+ stringify(htk_word_arc->start_htk_node) + ", "
				//		+ stringify(htk_word_arc->end_htk_node) + ", "
				//		+ stringify(htk_word_arc->ac_weight) + ", "
				//		+ stringify(htk_word_arc->lm_weight) + ") = "
				//		+ stringify(htk_arcs->size() - 1));
			}
			finish_word_arcs = true;
		}
	}
};

/*
 * The class for converting an FST lattice into an HTK lattice
 */
class FST2HTK_lat {

protected:
	static constexpr double sec_per_frame = 0.01;
	vector<FST_lat_node*> htk_nodes;
	vector<HTK_word_arc*> htk_arcs;
	map<StateId, FST_lat_node*> fst_node_map;
	set<StateId> queued_fst_nodes;

public:
//	~FST2HTK_lat();
//	FST_lat_node* findFstNode(StateId state_id);
//	FST_lat_node* insertFstNode(StateId state_id);
//	void convert(const VectorFst<StdArc>& fst_lat);
//	void Write(const string &filename, const SymbolTable& oSymTab) const;

	~FST2HTK_lat() {
		// the FST_lat_node objects in htk_nodes are a subject of those in fst_node_map.
		// so just need to delete those in fst_node_map.
		for (map<StateId, FST_lat_node*>::iterator nodeIt = fst_node_map.begin(); nodeIt != fst_node_map.end(); ++nodeIt) {
			delete nodeIt->second;
		}
		for (uint i = 0; i < htk_arcs.size(); ++i) {
			delete htk_arcs[i];
		}
	}

	FST_lat_node* findFstNode(StateId state_id) {
		map<StateId, FST_lat_node*>::iterator nodeIt = fst_node_map.find(state_id);
		if (nodeIt == fst_node_map.end()) return NULL;
		return nodeIt->second;
	}

	FST_lat_node* insertFstNode(StateId state_id) {
		if (findFstNode(state_id)) {
			cerr << "FST2HTK_lat::insertFstNode() ERROR: state " << state_id << " is already in the fst_node_map." << endl;
			exit(-1);
		}
		fst_node_map[state_id] = new FST_lat_node(state_id);
		return fst_node_map[state_id];
	}

	FST_lat_node* findOrInsertFstNode(StateId state_id, int label, int time_frame) {
		FST_lat_node* fst_node = findFstNode(state_id);
		if (!fst_node) {
			fst_node = insertFstNode(state_id);
			fst_node->word_lab = label;
			fst_node->time_frame = time_frame;
		} else {
			if (fst_node->word_lab != label) {
				cerr << "FST2HTK_lat::findOrInsertFstNode() ERROR: two incoming arcs going through the fst state " << state_id
						<< " with different word labels: " << fst_node->word_lab << " and " << label << endl;
				exit(-1);
			}
			if (fst_node->time_frame != time_frame) {
				cerr << "FST2HTK_lat::findOrInsertFstNode() ERROR: two paths reach the fst state " << state_id
						<< " at different time frame: " << fst_node->time_frame << " and " << time_frame << endl;
				exit(-1);
			}
		}
		return fst_node;
	}

	void convert(const VectorFst<StdArc>& fst_lat) {
		deque<StateId> state_list;
		StateId state_id = fst_lat.Start();
		if (state_id == kNoStateId) // empty fst
			return;
		state_list.push_back(state_id);
		// insert the start fst node to the map with the time frame as BAD_TIME_FRAME = -1
		FST_lat_node* start_node = insertFstNode(state_id);
		//start_node->in_queue = true;
		while (!state_list.empty()) {
			StateId state_id = state_list.front(); state_list.pop_front();
			//log_msg("FST2HTK_lat::convert() LOG: pop_front state_id " + stringify(state_id));
			FST_lat_node* fst_node = findFstNode(state_id);
			assert(fst_node);
			bool has_subsequent_non_eps_ilabels = false;
			for (ArcIterator<StdFst> aiter(fst_lat, state_id); !aiter.Done(); aiter.Next())
			{
				const StdArc &arc = aiter.Value();
				//log_msg("FST2HTK_lat::convert() LOG: state_id=" + stringify(state_id) +
				//		", arc.ilabel=" + stringify(arc.ilabel) +
				//		", arc.olabel=" + stringify(arc.olabel) +
				//		", arc.weight=" + stringify(arc.weight.Value()) +
				//		", arc.nextstate=" + stringify(arc.nextstate));

				// TODO: expand all nodes with epsilon labels
				if (arc.ilabel == 0) {
					cerr << "FST2HTK_lat::convert() ERROR: currently doesn't support epsilon input labels on any fst arc." << endl;
					exit(-1);
				} else {
					has_subsequent_non_eps_ilabels = true;
					StateId next_state_id = arc.nextstate;
					// avoid pushing duplicate fst nodes into the queue
					if (queued_fst_nodes.find(next_state_id) == queued_fst_nodes.end()) {
						state_list.push_back(next_state_id);
						queued_fst_nodes.insert(next_state_id);
						//log_msg("FST2HTK_lat::convert() LOG: push_back state_id " + stringify(next_state_id));
					}
					FST_lat_node* next_fst_node;
					if (arc.olabel == 0) {
						// next_fst_node is a word internal state
						next_fst_node = findOrInsertFstNode(next_state_id, fst_node->word_lab, fst_node->time_frame + 1);

						// Currently the FST only supports one type of weight, so all the weights are counted as acoustic weights as a hack
						// FST uses negative log weight for the tropical semiring.
						// But we need to use the regular log weight in HTK lattices. So we multiply the weights by -1.
						// TODO: make the FST support acoustic weight and LM weight.
						next_fst_node->add_word_internal_arc(*fst_node, -1 * arc.weight.Value(), -1 * 0.0);
					} else {
						// next_fst_node is a word beginning state, create an HTK arc for the previous word
						next_fst_node = findOrInsertFstNode(next_state_id, arc.olabel, fst_node->time_frame + 1);

						// currently the FST only supports one type of weight, so all the weights are counted as acoustic weights as a hack
						// FST uses negative log weight for the tropical semiring.
						// But we need to use the regular log weight in HTK lattices. So we multiply the weights by -1.
						// TODO: make the FST support acoustic weight and LM weight.
						next_fst_node->add_word_beginning_arc(&htk_nodes, fst_node, -1 * arc.weight.Value(), -1 * 0.0);

						// since the new word already starts from next_fst_node, all the partial word arcs in
						// the previous node fst_node can be all finished as full arcs at the time of fst_node.
						fst_node->add_word_arcs_to_htk_lat(&htk_nodes, &htk_arcs);
					}
				}
			}
			if (!has_subsequent_non_eps_ilabels) {
				// That means this fst_node does not have subsequent non-epsilon input labels
				// Then all the word arcs in it end at this node. Add them to the HTK lattice.
				fst_node->add_word_arcs_to_htk_lat(&htk_nodes, &htk_arcs);
			}
		}
	}

	void Write(const string &filename, const string &uttname, SymbolTable* oSymTab) const {
		try {
			ofstream htk_lat_stream;
			htk_lat_stream.open(filename.c_str());

			// print the header
			htk_lat_stream << "VERSION=1.0" << endl;
			htk_lat_stream << "UTTERANCE=" << uttname << endl;
			htk_lat_stream << "lmscale=1.00  wdpenalty=0.00" << endl;
			htk_lat_stream << "prscale=1.00" << endl;
			htk_lat_stream << "acscale=1.00" << endl;
			htk_lat_stream << "N=" << htk_nodes.size() << " L=" << htk_arcs.size() << endl;

			// print the nodes
			// I=0    t=0.00  W=!NULL
			// I=1    t=0.05  W=<s>                 v=1
			for (uint i = 0; i < htk_nodes.size(); ++i) {
				string label;
				if (htk_nodes[i]->word_lab == kNoLabel) {
					label = "!NULL";
				} else if (oSymTab != NULL && htk_lat_stream.is_open()) {
					label = oSymTab->Find(htk_nodes[i]->word_lab);
				} else {
					cerr << "CRFDecode ERROR: output symbol table has not been set or htk_lat_stream is already closed." << endl;
					exit(-1);
				}
				htk_lat_stream << "I=" << i
						<< " t=" << sec_per_frame * (htk_nodes[i]->time_frame + 1)
						<< " W=" << label;
				if (htk_nodes[i]->word_lab != kNoLabel) {
					htk_lat_stream << " v=1";
				}
				htk_lat_stream << endl;
			}

			// print the arcs
			// J=0     S=0    E=1    a=-263.35   l=0.000   r=0.00
			// J=1     S=0    E=2    a=-321.15   l=0.000   r=0.00
			for (uint j = 0; j < htk_arcs.size(); ++j) {
				htk_lat_stream << "J=" << j
						<< " S=" << htk_arcs[j]->start_htk_node
						<< " E=" << htk_arcs[j]->end_htk_node
						<< " a=" << htk_arcs[j]->ac_weight
						<< " l=" << htk_arcs[j]->lm_weight
						<< " r=0.00" << endl;
			}
		} catch (exception &e) {
			cerr << "Exception: " << e.what() << endl;
			exit(0);
		}
	}
};


/*
 * Main decoding block
 *
 */
int main(int argc, const char* argv[]) {
	char* progname;

	bool alignMode=false;
	bool mlfMode=false;
	SymbolTable* oSymTab = NULL;
	SymbolTable* iSymTab = NULL;

	VectorFst<StdArc>* lm_fst = NULL;

	set_defaults();
	QN_initargs(&argtab[0], &argc, &argv, &progname);
	QN_printargs(NULL, progname, &argtab[0]);
	//QN_logger = new QN_Logger_Simple(stdout,stderr,progname);

	cout << "IN LABELS: " << config.crf_label_size << endl;

	// Added by Ryan
	if (config.crf_lm_bin != NULL && strcmp(config.crf_lm_bin, "") == 0)
		config.crf_lm_bin = NULL;

	// Added by Ryan
	if (config.crf_conf_bin != NULL && strcmp(config.crf_conf_bin, "") == 0)
		config.crf_conf_bin = NULL;

	// Added by Ryan
	if (config.crf_disambig != NULL && strcmp(config.crf_disambig, "") == 0)
		config.crf_disambig = NULL;

	// Added by Ryan, for segmental CRFs
	if (config.label_maximum_duration <= 0)
	{
		string errstr="main() in CRFFstDecode caught exception: the maximum duration of labels must be larger than 0.";
		throw runtime_error(errstr);
	}
	if (config.num_actual_labs == 0)
	{
		config.num_actual_labs = config.crf_label_size;
	}
	if (strcmp(config.crf_model_type,"stdseg_no_dur_no_transftr")==0 &&
				strcmp(config.crf_featuremap, "stdstate") != 0)
	{
		string errstr="main() in CRFFstDecode caught exception: crf_featuremap must be \"stdstate\" for \"stdseg_no_dur_no_transftr\" CRF model.";
		throw runtime_error(errstr);
	}

	// commented out by Ryan for the class CRF_StdSegNode_WithoutDurLab
//	if (config.crf_label_size != config.label_maximum_duration * config.num_actual_labs)
//	{
//		string errstr="main() in CRFFstDecode caught exception: It should be crf_label_size == label_maximum_duration * num_actual_labs.";
//		throw runtime_error(errstr);
//	}


	if (strcmp(config.crf_decode_mode,"align")==0) { alignMode=true; }

    if (strcmp(config.hardtarget_file,"")==0) { config.hardtarget_file=NULL; }

    if ((config.hardtarget_file==NULL) && (alignMode)) {
    	cerr << "hardtarget_file required when crf_decode_mode=align" << endl;
    	exit(-1);
    }

    if ((config.crf_output_labelfile==NULL) && (config.crf_output_mlffile==NULL)) {
    	cerr << "At least one of crf_output_labelfile or crf_output_mlffile must be assigned" << endl;
    	exit(-1);
    }


	FILE* outl=NULL;
	QN_OutLabStream* labout=NULL;
	if (config.crf_output_labelfile != NULL ) {
		outl=fopen(config.crf_output_labelfile,"w+");

		// Changed by Ryan
		//labout = new QN_OutLabStream_ILab(1, "out_labfile", outl, 255, 1);
		labout = new QN_OutLabStream_ILab(1, "out_labfile", outl, 8192, 1);
	}

	CRF_MLFManager* mlfManager = NULL;
    ofstream mlfstream;
    vector<string> olist;
    if (config.crf_output_mlffile != NULL) {
    	mlfMode=true;
    	if (config.crf_olist == NULL) {
    		cerr << "crf_olist required when crf_output mlffile defined" << endl;
    		exit(-1);
    	}
    	ifstream olistfile;
    	olistfile.open(config.crf_olist);
    	if (olistfile.is_open()) {
    		while (!olistfile.eof()) {
    			string s;
    			getline(olistfile,s);
    			if (s!="") {
    			olist.push_back(s);
    			}}
    	}
    	olistfile.close();
    	if (config.crf_osymbols==NULL) {
    		cerr << "crf_osymbols required when crf_output_mlffile defined" << endl;
    		exit(-1);
    	}
    	oSymTab = SymbolTable::ReadText(config.crf_osymbols);
    	if (config.crf_isymbols==NULL ) {
    		cerr << "crf_osymbols required when crf_output_mlffile defined" << endl;
    		exit(-1);
    	}
    	iSymTab = SymbolTable::ReadText(config.crf_isymbols);
    	mlfstream.open(config.crf_output_mlffile);
    	mlfstream << "#!MLF!#" << endl;

        if (config.crf_align_mlffile != NULL) {
        	try {
        	mlfManager=new CRF_MLFManager(config.crf_align_mlffile, config.crf_olist, oSymTab);
        	}
    		catch (exception &e) {
    			cerr << "Exception: " << e.what() << endl;
    		}
        }
    }


	if (config.crf_lm_bin != NULL) {
		cout << "Reading in LM fst from file: " << config.crf_lm_bin << endl;
		VectorFst<LogArc>* log_fst=VectorFst<LogArc>::Read(config.crf_lm_bin);
		if (log_fst == NULL) {
			cerr << "ERROR: Failed opening file: " << config.crf_lm_bin << endl;
			exit(-1);
		}

		lm_fst=new VectorFst<StdArc>();
		log_msg("Mapping LM to Tropical Ring...");
		Map(*log_fst,lm_fst,LogToStdMapper());
		delete log_fst;

	    vector<int> disambig_list;
		if (config.crf_disambig == NULL) {
			cout << "crf_disambig is not set: no disambiguation symbols for the decoding graph." << endl;
		} else {
			ifstream disambig_file;
			disambig_file.open(config.crf_disambig);
			if (disambig_file.is_open()) {
				while (!disambig_file.eof()) {
					string s;
					getline(disambig_file, s);
					if (s != "") {
						int disambig_id = atoi(s.c_str());
						if (disambig_id == 0) {
							cerr << "ERROR: invalid disambiguation ID (" << s << ") from "
									<< config.crf_disambig << endl;
							exit(-1);
						}
						disambig_list.push_back(disambig_id);
					}
				}
			}
			disambig_file.close();
		}

		// changed by Ryan, originally commented out.
		// It's a big hack for WSJ to hard-code 55, 56, 57, 58.
		// Wouldn't work for other dataset. Need a true solution.
		//log_msg("Removing EOW tokens from final PhoneMap/LM/Dict lattice");
		log_msg("Removing disambiguation symbols from final PhoneMap/LM/Dict lattice");
		vector<pair<StdArc::Label, StdArc::Label> >* ipairs;
		vector<pair<StdArc::Label, StdArc::Label> >* opairs;
		ipairs = new vector<pair<StdArc::Label, StdArc::Label> >;
		opairs = new vector<pair<StdArc::Label, StdArc::Label> >;
		//pair<StdArc::Label,StdArc::Label> mypair1(55,0);
		//pair<StdArc::Label,StdArc::Label> mypair2(56,0);
		//pair<StdArc::Label,StdArc::Label> mypair3(57,0);
		//pair<StdArc::Label,StdArc::Label> mypair4(58,0);
		//(*ipairs).push_back(mypair1);
		//(*ipairs).push_back(mypair2);
		//(*ipairs).push_back(mypair3);
		//(*ipairs).push_back(mypair4);
		for (uint i = 0; i < disambig_list.size(); ++i) {
			pair<StdArc::Label, StdArc::Label> mypair(disambig_list[i], 0);
			(*ipairs).push_back(mypair);
			stringstream ss;
			ss << "Disambiguation symbol ID: " << disambig_list[i];
			log_msg(ss.str());
		}
		Relabel(lm_fst,*ipairs,*opairs);


	}

    ftrmaptype trn_ftrmap = STDSTATE;

	if (strcmp(config.crf_featuremap,"stdtrans")==0) { trn_ftrmap=STDTRANS;}
	if (strcmp(config.crf_featuremap,"stdsparse")==0) { trn_ftrmap=STDSPARSE;}
	if (strcmp(config.crf_featuremap,"stdsparsetrans")==0) { trn_ftrmap=STDSPARSETRANS;}
	if (strcmp(config.crf_featuremap,"file")==0) { trn_ftrmap=INFILE;}

	if ((trn_ftrmap == INFILE) && (config.crf_featuremap_file==NULL)) {
		cerr << "ERROR: crf_featuremap_file must be non-NULL when crf_featuremap set to 'file'" << endl;
		exit(-1);
	}

	CRF_FeatureStreamManager* str2=NULL;
	CRF_FeatureStreamManager* str3=NULL;

	// modified by Ryan, for context features
	CRF_FeatureStreamManager str1(1,"ftr1_file",config.ftr1_file,config.ftr1_format,config.hardtarget_file,config.hardtarget_window_offset,
							(size_t) config.ftr1_width, (size_t) config.ftr1_ftr_start, (size_t) config.ftr1_ftr_count,
							config.window_extent, config.ftr1_window_offset, config.ftr1_window_len,
							config.ftr1_left_context_len, config.ftr1_right_context_len, config.ftr1_extract_seg_ftr,
							config.ftr1_use_boundary_delta_ftr,
							config.ftr1_delta_order, config.ftr1_delta_win,
							config.crf_eval_range, 0,
							NULL,0,0,0,SEQUENTIAL);
	if (strcmp(config.ftr2_file,"") != 0) {
		str2=new CRF_FeatureStreamManager(1,"ftr2_file",config.ftr2_file,config.ftr2_format,config.hardtarget_file,config.hardtarget_window_offset,
							(size_t) config.ftr2_width, (size_t) config.ftr2_ftr_start, (size_t) config.ftr2_ftr_count,
							config.window_extent, config.ftr2_window_offset, config.ftr2_window_len,
							config.ftr2_left_context_len, config.ftr2_right_context_len, config.ftr2_extract_seg_ftr,
							config.ftr2_use_boundary_delta_ftr,
							config.ftr2_delta_order, config.ftr2_delta_win,
							config.crf_eval_range, 0,
							NULL,0,0,0,SEQUENTIAL);
		str1.join(str2);
	}
	if (strcmp(config.ftr3_file,"") != 0) {
		str3=new CRF_FeatureStreamManager(1,"ftr3_file",config.ftr3_file,config.ftr3_format,config.hardtarget_file,config.hardtarget_window_offset,
							(size_t) config.ftr3_width, (size_t) config.ftr3_ftr_start, (size_t) config.ftr3_ftr_count,
							config.window_extent, config.ftr3_window_offset, config.ftr3_window_len,
							config.ftr3_left_context_len, config.ftr3_right_context_len, config.ftr3_extract_seg_ftr,
							config.ftr3_use_boundary_delta_ftr,
							config.ftr3_delta_order, config.ftr3_delta_win,
							config.crf_eval_range, 0,
							NULL,0,0,0,SEQUENTIAL);
		str1.join(str3);

	}
    cout << "Feature File created" << endl;


	CRF_FeatureStream* crf_ftr_str = str1.trn_stream;
	CRF_Model my_crf(config.crf_label_size);
	cout << "LABELS: " << my_crf.getNLabs() << endl;

	// Added by Ryan, for segmental CRFs
	my_crf.setLabMaxDur(config.label_maximum_duration);
	my_crf.setNActualLabs(config.num_actual_labs);
	cout << "LABEL_MAXIMUM_DURATION: " << my_crf.getLabMaxDur() << endl;
	cout << "ACTUAL_LABELS: " << my_crf.getNActualLabs() << endl;

	// Added by Ryan, for segmental CRFs
	modeltype mtype;
	if (strcmp(config.crf_model_type,"stdframe")==0)
	{
			mtype = STDFRAME;
	}
	else if (strcmp(config.crf_model_type,"stdseg")==0)
	{
			mtype = STDSEG;
	}
	else if (strcmp(config.crf_model_type,"stdseg_no_dur")==0)
	{
			mtype = STDSEG_NO_DUR;
	}
	else if (strcmp(config.crf_model_type,"stdseg_no_dur_no_transftr")==0)
	{
			mtype = STDSEG_NO_DUR_NO_TRANSFTR;
	}
	else if (strcmp(config.crf_model_type,"stdseg_no_dur_no_segtransftr")==0)
	{
			mtype = STDSEG_NO_DUR_NO_SEGTRANSFTR;
	}
	else
	{
			mtype = STDFRAME;
	}
	my_crf.setModelType(mtype);
	cout << "MODEL_TYPE: " << config.crf_model_type << endl;

	// Added by Ryan, for segmental CRFs
	if (my_crf.getModelType() == STDFRAME && my_crf.getLabMaxDur() != 1)
	{
		string errstr="main() in CRFFstDecode caught exception: the maximum duration of labels must be 1 for \"stdframe\" CRF model.";
		throw runtime_error(errstr);
	}

	set_fmap_config(str1.getNumFtrs());
	my_crf.setFeatureMap(CRF_FeatureMap::createFeatureMap(&fmap_config));

	bool openchk=my_crf.readFromFile(config.weight_file);
	if (!openchk) {
		cerr << "ERROR: Failed opening file: " << config.weight_file << endl;
		exit(-1);
	}

	// just for debugging
	//cout << "Before processing segments..." << endl;

	// Added by Ryan. We are currently requiring the olist file since we need them to print out the 1-best or lattices.
	// TODO: Set it optional whether to output the 1-best or lattices.
	if (config.crf_olist == NULL) {
		cerr << "crf_olist required currently." << endl;
		exit(-1);
	}

	//CRF_ViterbiDecoder vd(crf_ftr_str,&my_crf);
	crf_ftr_str->rewind();
	QN_SegID segid = crf_ftr_str->nextseg();
	int count=0;
	QN_Range eval_sent_range(config.crf_eval_range);
	QN_Range::iterator eval_sent_range_iter = eval_sent_range.begin();
	while (segid != QN_SEGID_BAD) {
		if (count %100 == 0) {
			cout << "Processing segment " << count << endl;
		}
		if (config.crf_olist != NULL) {
			//cout << "Processing file: " << olist.at(count);
			//cout << " (" << count << ")" << endl;
			if (*eval_sent_range_iter >= olist.size()) {
				string errstr="main() in CRFDecode caught exception: "
						"eval sentence range goes out of the olist size.";
				throw runtime_error(errstr);
			}
			cout << "Processing file: " << olist.at(*eval_sent_range_iter);
			cout << " (" << *eval_sent_range_iter << " in the olist)";
			cout << " (" << count << " in the current test set)" << endl;
		}
		try {
			time_t rawtime;
			time(&rawtime);
			char* time1 = ctime(&rawtime);
			time1[strlen(time1)-1]='\0';
			log_msg(time1);

			// Changed by Ryan, for segmental CRFs
			////////////////////////////////////
			// the frame-level viterbi decoder
			////////////////////////////////////
////			CRF_ViterbiDecoder* vd = new CRF_ViterbiDecoder(crf_ftr_str,&my_crf);
			////////////////////////////////////
			// the segmental viterbi decoder
			////////////////////////////////////
			CRF_ViterbiDecoder_StdSeg_NoSegTransFtr<CRF_ViterbiNode>* vd;
			if (my_crf.getModelType() == STDFRAME ||
					my_crf.getModelType() == STDSEG_NO_DUR_NO_SEGTRANSFTR)
			{
				vd = new CRF_ViterbiDecoder_StdSeg_NoSegTransFtr<CRF_ViterbiNode>(crf_ftr_str,&my_crf);
			}
			else
			{
				string errstr="main() in CRFFstDecode caught exception: "
						"CRF_ViterbiDecoder for CRF models other than \"stdframe\" "
						"and \"stdseg_no_dur_no_segtransftr\" have not been implmented.";
				throw runtime_error(errstr);
			}
			/////////////////////////////////////////////////////////////
			// the segmental viterbi decoder with transition arc pruned
			/////////////////////////////////////////////////////////////
//			CRF_ViterbiDecoder_StdSeg_NoSegTransFtr<CRF_ViterbiNode_PruneTrans>* vd;
//			if (my_crf.getModelType() == STDFRAME ||
//					my_crf.getModelType() == STDSEG_NO_DUR_NO_SEGTRANSFTR)
//			{
//				vd = new CRF_ViterbiDecoder_StdSeg_NoSegTransFtr<CRF_ViterbiNode_PruneTrans>(crf_ftr_str,&my_crf);
//			}
//			else
//			{
//				string errstr="main() in CRFFstDecode caught exception: "
//						"CRF_ViterbiDecoder for CRF models other than \"stdframe\" "
//						"and \"stdseg_no_dur_no_segtransftr\" have not been implmented.";
//				throw runtime_error(errstr);
//			}

			vd->setIfOutputFullFst(config.crf_if_output_full_lat);

			VectorFst<StdArc>* best_lat=new VectorFst<StdArc>();

			// Added by Ryan
			// the fully composed lattice, including acoustic model, dictionary, and language model.
			VectorFst<StdArc>* out_full_lat = new VectorFst<StdArc>();

			// just for debugging
			//cout << "Before nStateDecode ..." << endl;

			// Changed by Ryan:
			// for CRF_ViterbiDecoder
//			int nodeCnt=vd->nStateDecode(best_lat,lm_fst,config.crf_decode_beam,config.crf_decode_min_hyp,
//					config.crf_decode_max_hyp,config.crf_decode_hyp_inc);
			// for CRF_ViterbiDecoder_StdSeg_NoSegTransFtr
			int nodeCnt=vd->nStateDecode(best_lat,lm_fst,out_full_lat,
					config.crf_decode_beam,config.crf_decode_min_hyp,
					config.crf_decode_max_hyp,config.crf_decode_hyp_inc);

			// just for debugging
			//cout << "After nStateDecode ..." << endl;

			if (config.crf_lat_outdir != NULL) {
				// Put code here to dump the lattice file out to the latdir in
				// OpenFst format
				//string fst_fname = string(config.crf_lat_outdir) + "/" + olist.at(count) + ".fst";
				if (*eval_sent_range_iter >= olist.size()) {
					string errstr="main() in CRFDecode caught exception: "
							"eval sentence range goes out of the olist size.";
					throw runtime_error(errstr);
				}
				string fst_fname = string(config.crf_lat_outdir) + "/" + olist.at(*eval_sent_range_iter) + ".fst";
				log_msg("*	*Writing FST lattice to "+fst_fname);

				// Changed by Ryan
				// output the fully-composed lattice
				if (config.crf_if_output_full_lat)
				{
					out_full_lat->Write(fst_fname);
				}
				// only output the best lattice
				else
				{
					best_lat->Write(fst_fname);
				}
			}

			// Added by Ryan
			// dump the lattice file out to the htk_lat_outdir in HTK SLF format
			if (config.htk_lat_outdir != NULL) {
				if (*eval_sent_range_iter >= olist.size()) {
					string errstr="main() in CRFDecode caught exception: "
							"eval sentence range goes out of the olist size.";
					throw runtime_error(errstr);
				}
				string utt_name = olist.at(*eval_sent_range_iter);
				//string slf_fname = string(config.htk_lat_outdir) + "/" + olist.at(count) + ".slf";
				string slf_fname = string(config.htk_lat_outdir) + "/" + utt_name + ".slf";
				log_msg("*	*Writing HTK lattice to " + slf_fname);

				// output the fully-composed lattice
				if (config.crf_if_output_full_lat)
				{
					FST2HTK_lat fst2htk_lat_converter;
					fst2htk_lat_converter.convert(*out_full_lat);
					fst2htk_lat_converter.Write(slf_fname, utt_name, oSymTab);
				}
				// only output the best lattice
				else
				{
					FST2HTK_lat fst2htk_lat_converter;
					fst2htk_lat_converter.convert(*best_lat);
					fst2htk_lat_converter.Write(slf_fname, utt_name, oSymTab);
				}
			}

			if (config.crf_output_mlffile != NULL) {
				VectorFst<StdArc>* shortest_fst = new VectorFst<StdArc>();
				log_msg("Finding Shortest Path");

				// Changed by Ryan
				// decode on the best lattice
				ShortestPath(*best_lat,shortest_fst,1);
				// decode on the fully-composed lattice
//				ShortestPath(*out_full_lat,shortest_fst,1);

				//log_msg("Projecting to Output Symbols");
				//Project(shortest_fst,PROJECT_OUTPUT);
				log_msg("Removing epsilons");
				RmEpsilon(shortest_fst);
				log_msg("Top sorting");
				TopSort(shortest_fst);

				//	Iterates over the FSTs states.
				if (oSymTab != NULL && mlfstream.is_open()) {
					 //mlfstream << "\"" << olist.at(count) << "\"" << endl;
					 //cout << "\"" << olist.at(count) << "\" : ";
					 if (*eval_sent_range_iter >= olist.size()) {
						 string errstr="main() in CRFDecode caught exception: "
								"eval sentence range goes out of the olist size.";
						 throw runtime_error(errstr);
					 }
					 mlfstream << "\"" << olist.at(*eval_sent_range_iter) << "\"" << endl;
					 cout << "\"" << olist.at(*eval_sent_range_iter) << "\" : ";
				}
				int phnStateStart=0;
				int wrdStart=0;
				int current=0;
				int lastId=-1;
				int lastPhnStateLabel=-1;
				int lastPhnActualLab=-1;
				int curWrd=-1;
				for (StateIterator<StdFst> siter(*shortest_fst); !siter.Done(); siter.Next())
				{
  					StateId state_id = siter.Value();
  					//Iterates over state state_id's arcs.
  					for (ArcIterator<StdFst> aiter(*shortest_fst, state_id); !aiter.Done(); aiter.Next())
					{
  						const StdArc &arc = aiter.Value();
  						if ((arc.olabel ==0) && (arc.ilabel != 0)) {
  							// In here we process an epsilon on the word path
  							if (config.crf_mlf_output_states) {
  								if (arc.ilabel != lastPhnStateLabel) {
  									if (lastPhnStateLabel != -1) {
  										// assuming the state/word label is at the end of the word path
  										if (config.crf_mlf_output_frames) {
  											//mlfstream << phnStateStart << "\t" << (current-1) << "\t";
  											mlfstream << phnStateStart << "\t" << current << "\t";  // start inclusive, end exclusive
  										}
  										mlfstream << iSymTab->Find(lastPhnStateLabel) << endl;
  									}
  									lastPhnStateLabel=arc.ilabel;

  									//
  									// Added by Ryan: currently only works for frame-level model, not segmental model yet!!!
  									//
  									lastPhnActualLab = (lastPhnStateLabel - 1) / config.crf_states;

  									phnStateStart=current;
  								}
  							}
  							current++;
  						}
  						else if (arc.olabel !=0 && arc.ilabel==0) {
  							// Here we process an epsilon on the phone path and a word
  							if (oSymTab != NULL && mlfstream.is_open()) {
  								// assuming the state/word label is at the end of the word path
  								if (config.crf_mlf_output_frames) {
  									//mlfstream << phnStateStart << "\t" << (current-1) << "\t";
									mlfstream << phnStateStart << "\t" << current << "\t";  // start inclusive, end exclusive
  								}
  								if (config.crf_mlf_output_states) {
  									mlfstream << iSymTab->Find(lastPhnStateLabel) << "\t";
  								}
  								mlfstream  << oSymTab->Find(arc.olabel) << endl;
  								cout << oSymTab->Find(arc.olabel) << " ";
  								phnStateStart=current;
  								lastPhnStateLabel=-1;

								//
								// Added by Ryan: currently only works for frame-level model, not segmental model yet!!!
								//
								lastPhnActualLab = -1;

  							}
  						}
  						else if (arc.olabel !=0 && arc.ilabel!=0) {
  							// And finally we process no epsilons on any line
							if (config.crf_mlf_output_states) {
  								if (arc.ilabel != lastPhnStateLabel) {
  									if (lastPhnStateLabel != -1) {
  										if (config.crf_mlf_output_frames) {
  											//mlfstream << phnStateStart << "\t" << (current-1) << "\t";
											mlfstream << phnStateStart << "\t" << current << "\t";  // start inclusive, end exclusive
  										}
  										mlfstream << iSymTab->Find(lastPhnStateLabel) << endl;
  									}
  									lastPhnStateLabel=arc.ilabel;

  									//
  									// Added by Ryan: currently only works for frame-level model, not segmental model yet!!!
  									//
  									lastPhnActualLab = (lastPhnStateLabel - 1) / config.crf_states;

  									phnStateStart=current;
  								}
  							}
							current++;
  							if (oSymTab != NULL && mlfstream.is_open()) {
  								if (config.crf_mlf_output_frames) {
  									// assuming the word label is at the end of the last phone
  									//mlfstream << phnStateStart << "\t" << (current-1) << "\t";
									mlfstream << phnStateStart << "\t" << current << "\t";  // start inclusive, end exclusive
  								}
  								if (config.crf_mlf_output_states) {
  									mlfstream << iSymTab->Find(lastPhnStateLabel) << "\t";
  								}
  								mlfstream  << oSymTab->Find(arc.olabel) << endl;
  								cout << oSymTab->Find(arc.olabel) << " ";
  								phnStateStart=current;
  								lastPhnStateLabel=-1;

								//
								// Added by Ryan: currently only works for frame-level model, not segmental model yet!!!
								//
								lastPhnActualLab = -1;

  							}
  							//current++;
  						}

  						/*
  						if (arc.olabel == 0) {
  							if (config.crf_mlf_output_states) {
  							if (arc.ilabel != lastId) {
  								if (lastId != -1) {
  									if (config.crf_mlf_output_frames) {
  										mlfstream << start << "\t" << (current-1) << "\t";
  									}
  									mlfstream << iSymTab->Find(lastId) << endl;
  								}
  								lastId=arc.ilabel;
  								start=current;
  							}
  							}
  							current++;
  						}
  						else {
  							if (oSymTab != NULL && mlfstream.is_open()) {
  								if (config.crf_mlf_output_frames) {
  									mlfstream << start << "\t" << (current-1) << "\t";
  								}
  								if (config.crf_mlf_output_states) {
  									mlfstream << iSymTab->Find(lastId) << "\t";
  								}
  								mlfstream  << oSymTab->Find(arc.olabel) << endl;
  								cout << oSymTab->Find(arc.olabel) << " ";
  								start=current;
  								lastId=-1;
  							}
  						}
  						*/
					}
				}
				if (oSymTab != NULL && mlfstream.is_open()) {
					if (lastPhnStateLabel != -1) {
							if (config.crf_mlf_output_states) {

									if (lastPhnStateLabel != -1) {
										if (config.crf_mlf_output_frames) {
											// assuming the state/word label is at the end of the word path
											//mlfstream << phnStateStart << "\t" << (current-1) << "\t";
											mlfstream << phnStateStart << "\t" << current << "\t";  // start inclusive, end exclusive
										}
										mlfstream << iSymTab->Find(lastPhnStateLabel) << endl;
									}

							}
					}
					mlfstream << "." << endl;
					cout << "." << endl;
				}
				delete shortest_fst;
			}

			time(&rawtime);
			char* time2 = ctime(&rawtime);
			time2[strlen(time2)-1]='\0';
			log_msg(time2);

			delete best_lat;
			delete vd;

			// Added by Ryan
			delete out_full_lat;
		}
		catch (exception &e) {
			cerr << "Exception: " << e.what() << endl;
			exit(0);
		}
		segid=crf_ftr_str->nextseg();
		count++;
		++eval_sent_range_iter;
	}
	delete lm_fst;

};
