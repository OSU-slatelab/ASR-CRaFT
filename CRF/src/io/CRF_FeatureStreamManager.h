#ifndef CRF_FEATURESTREAMMANAGER_H_
#define CRF_FEATURESTREAMMANAGER_H_
/*
 * CRF_FeatureStreamManager.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_FeatureStreamManager
 * Uses the feature stream model/classes defined for the ICSI Quicknet neural networks
 * package for compatibility with ICSI Quicknet.
 */

#include <QuickNet.h>
#include <iostream>
#include <stdio.h>

#include "../CRF.h"
#include "CRF_InFtrStream_RandPresent.h"
#include "CRF_InLabStream_RandPresent.h"
#include "CRF_FeatureStream.h"
//added by Ryan
#include "CRF_InFtrStream_SeqMultiWindow.h"
#include "CRF_InLabStream_SeqMultiWindow.h"

/*
 * class CRF_FeatureStreamManager
 *
 * Uses the underlying CRF_FeatureStream objects to define and control access to training and
 * cross-validation portions of the feature stream.  Acts as an interface to the various CRF_FeatureStream
 * classes.
 *
 * It is recommended to use this class rather than directly using one of the various CRF_FeatureStream
 * classes.
 */

class CRF_FeatureStreamManager
{
private:
	int debug;
	const char* dbgname;
	char* filename;
	const char* format;
	char* hardtarget_filename;
	size_t hardtarget_window_offset;
	size_t width;
	size_t first_ftr;
	size_t num_ftrs;
	size_t window_extent;
	size_t window_offset;
	size_t window_len;

	// added by Ryan, for context features
	size_t left_context_len;
	size_t right_context_len;
	bool extract_segment_features;
	bool use_boundary_delta_ftrs;

	int delta_order;
	int delta_win;
	char* train_sent_range;
	char* cv_sent_range;
	FILE* normfile;
	int norm_mode;
	double norm_am;
	double norm_av;
	size_t train_cache_frames;
	int train_cache_seed;
	seqtype train_seq_type;
	QNUInt32 rseed;
	size_t nthreads;
	CRF_FeatureStreamManager **children;
protected:
	int childnum;
public:
	// modified by Ryan, for context features
//	CRF_FeatureStreamManager(int debug, const char* debug_name,
//										char* ftr_fname, const char* ftr_file_fmt, char* ht_fname, size_t ht_offset,
//										size_t ftr_width, size_t first_ftr, size_t num_ftrs,
//										size_t win_ext, size_t win_off, size_t win_len,
//										int delta_o, int delta_w,
//										char* trn_rng, char* cv_rng,
//										FILE* nfile, int n_mode, double n_am, double n_av, seqtype ts,
//										QNUInt32 rseed=0, size_t n_threads=1);
	CRF_FeatureStreamManager(int debug, const char* debug_name,
									char* ftr_fname, const char* ftr_file_fmt, char* ht_fname, size_t ht_offset,
									size_t ftr_width, size_t first_ftr, size_t num_ftrs,
									size_t win_ext, size_t win_off, size_t win_len,
									size_t left_ctx_len, size_t right_ctx_len, bool extract_seg_ftr,
									bool use_bdy_delta_ftr,
									int delta_o, int delta_w,
									char* trn_rng, char* cv_rng,
									FILE* nfile, int n_mode, double n_am, double n_av, seqtype ts,
									QNUInt32 rseed=0, size_t n_threads=1);

	virtual ~CRF_FeatureStreamManager();
	void create();
	void display();
	void join(CRF_FeatureStreamManager*);
	void setUtt(QNUInt32);
	void setFiles(char*, const char*, char*, size_t, size_t, size_t);
	void setNorm(FILE*, int, double, double);
	void setRanges(char*, char*, size_t, int);
	void setWindow(size_t, size_t, size_t);
	void setDeltas(int,int);
	size_t getNumFtrs();
	int setDebug(int, const char*);
	inline CRF_FeatureStreamManager *getChild(size_t child) {
		return (children==NULL || child>=nthreads || child<0)?NULL:children[child];
	}
	inline size_t getNThreads() { return nthreads; }
	CRF_FeatureStream* trn_stream;
	CRF_FeatureStream* cv_stream;
	CRF_FeatureStream* old_trn_stream;


};

#endif /*QN_FEATURESTREAM_H_*/
