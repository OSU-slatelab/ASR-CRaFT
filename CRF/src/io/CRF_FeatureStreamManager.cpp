/*
 * CRF_FeatureStreamManager.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */

#include "CRF_FeatureStreamManager.h"

// ERROR & PROCESS LOGGING NEEDS TO BE ADDED TO THIS CLASS

/*
 * CRF_FeatureStreamManager constructor
 *
 * Input: dbg - debug flag
 *        dname - debug file anme
 *        fname - feature file name
 *        fmt - format of feature file
 *        ht_fname - label file name
 *        ht_offset - offset into label window for training
 *        width -
 *        first_ftr - offset into feature vector to start features
 *        num_ftrs - total number of features in feature vector to use, starting at first_ftr
 *        win_ext - size of window of feature vectors to use
 *        win_off - offset into window to consider "center" feature
 *        win_len - size of window of feature vectors to use
 *        delta_o - delta order to apply to feature file (not yet implemented)
 *        delta_w - window size to compute deltas (not yet implemented)
 *        trn_rng - range of feature vectors to use for training
 *        cv_rng - range of feature vectors to use for cross-validation
 *        nfile - file containing normalization terms (not yet implemented)
 *        n_mode - normalization mode (not yet implemented)
 *        n_am - (not yet implemented)
 *        n_av - (not yet implemented)
 *        train-seq_type - type of training to perform (sequential, random no replace, random replace)
 *        rseed - random seed
 *        nthreads - number of threads (used to split the training file over multiple views for multi-
 *                    threaded training aglorithms)
 *
 *  Normalization and delta parameters are included for future expansion for compatability with Quicknet
 *  feature streams, but are not implemented at this time.
 *
 */
// modified by Ryan, for context features
//CRF_FeatureStreamManager::CRF_FeatureStreamManager(int dbg, const char* dname,
//									char* fname, const char* fmt, char* ht_fname, size_t ht_offset,
//									size_t width, size_t first_ftr, size_t num_ftrs,
//									size_t win_ext, size_t win_off, size_t win_len,
//									int delta_o, int delta_w,
//									char* trn_rng, char* cv_rng,
//									FILE* nfile, int n_mode, double n_am, double n_av, seqtype ts,
//									QNUInt32 rseed, size_t n_threads)
CRF_FeatureStreamManager::CRF_FeatureStreamManager(int dbg, const char* dname,
									char* fname, const char* fmt, char* ht_fname, size_t ht_offset,
									size_t width, size_t first_ftr, size_t num_ftrs,
									size_t win_ext, size_t win_off, size_t win_len,
									size_t left_ctx_len, size_t right_ctx_len, bool extract_seg_ftr,
									bool use_bdy_delta_ftr,
									int delta_o, int delta_w,
									char* trn_rng, char* cv_rng,
									FILE* nfile, int n_mode, double n_am, double n_av, seqtype ts,
									QNUInt32 rseed, size_t n_threads)
	:debug(dbg),
	 dbgname(dname),
	 filename(fname),
	 format(fmt),
	 hardtarget_filename(ht_fname),
	 hardtarget_window_offset(ht_offset),
	 width(width),
	 first_ftr(first_ftr),
	 num_ftrs(num_ftrs),
	 window_extent(win_ext),
	 window_offset(win_off),
	 window_len(win_len),

	 // added by Ryan, for context features
	 left_context_len(left_ctx_len),
	 right_context_len(right_ctx_len),
	 extract_segment_features(extract_seg_ftr),
	 use_boundary_delta_ftrs(use_bdy_delta_ftr),

	 delta_order(delta_o),
	 delta_win(delta_w),
	 train_sent_range(trn_rng),
	 cv_sent_range(cv_rng),
	 normfile(nfile),
	 norm_mode(n_mode),
	 norm_am(n_am),
	 norm_av(n_av),
	 train_seq_type(ts),
	 rseed(rseed),
	 nthreads(n_threads)
{
	childnum=0;
	this->create();
}

/*
 * CRF_FeatureStreamManager destructor
 */
CRF_FeatureStreamManager::~CRF_FeatureStreamManager()
{
	if (this->trn_stream) {
		// the trn_stream feature file is owned by this manager
		delete this->trn_stream->ftr_stream;
		delete this->trn_stream->lab_stream;
		delete this->trn_stream;
	}
	if (this->cv_stream) {
		delete this->cv_stream->ftr_stream;
		delete this->cv_stream->lab_stream;
		delete this->cv_stream;
	}
	if (children) {
		for(int i=0;i<nthreads;i++) {
			delete children[i];
		}
		delete [] children;
	}
}

/*
 * CRF_FeatureStreamManager:create
 *
 * Takes the internal variables of the FeatureStreamManager and uses them to create CRF_FeatureStreams
 * under a variety of conditions.  Separates training and cross-validation streams, creates multiple
 * streams if the Manager was created with the multi-threaded option, determines if streams should be
 * accessed sequentially, randomly, or random with replacement, and associates the label file with
 * streams in the proper way.
 */
void CRF_FeatureStreamManager::create()
{
	QN_InFtrStream* ftr_str = NULL;     // Temporary stream holder.
    int index = 1;                      // training always requires indexed
    int buffer_frames = 500;
    QNUInt32 nseg;

    if (nthreads<=1) {
    	children=NULL;
    } else {
    	children=new CRF_FeatureStreamManager *[nthreads];
    	// do initiaization later
    }

    ftr_str = QN_build_ftrstream(this->debug, this->dbgname, this->filename, this->format,
                              this->width, index, this->normfile,
                              this->first_ftr, this->num_ftrs,
                              0, QN_ALL, // do utt selection ourselves
                              buffer_frames,
                              this->delta_order, this->delta_win,
                              this->norm_mode, this->norm_am, this->norm_av);

    // Create training and cross-validation streams.
    QN_InFtrStream* train_ftr_str = NULL;
    QN_InFtrStream* cv_ftr_str = NULL;

    // Using range strings
    if (cv_sent_range==0) {
    	train_ftr_str=ftr_str;
    	QN_InFtrStream_CutRange* fwd_ftr_str
    			= new QN_InFtrStream_CutRange(debug, dbgname, *ftr_str,
            	                              train_sent_range,
                	                          0);
		QN_InFtrStream_Cut* train_ftr_str_cut = (QN_InFtrStream_Cut*)fwd_ftr_str;
		train_ftr_str = (QN_InFtrStream*) train_ftr_str_cut;
		cv_ftr_str=NULL;
    }
    else {
    	QN_InFtrStream_CutRange* fwd_ftr_str
        	    = new QN_InFtrStream_CutRange(debug, dbgname, *ftr_str,
            	                              train_sent_range,
                	                          cv_sent_range);
    	QN_InFtrStream_Cut* train_ftr_str_cut = (QN_InFtrStream_Cut*)fwd_ftr_str;
    	cv_ftr_str = (QN_InFtrStream*) new QN_InFtrStream_Cut2(*train_ftr_str_cut);
    	train_ftr_str = (QN_InFtrStream*) train_ftr_str_cut;
    }

    nseg=train_ftr_str->num_segs();

    // Added by Ryan
    if (this->window_offset + this->window_len > this->window_extent)
    {
    	string errstr="CRF_FeatureStreamManager::create() caught exception: this->window_offset + this->window_len > this->window_extent.";
		throw runtime_error(errstr);
    }

    // Create training and CV windows.
    size_t bot_margin = this->window_extent - this->window_offset - this->window_len;

    //changed by Ryan
//    QN_InFtrStream_SeqWindow* train_winftr_str;
    CRF_InFtrStream_SeqMultiWindow* train_winftr_str;

    CRF_InFtrStream_RandPresent* train_randftr_str = NULL;

	switch ( this->train_seq_type )
	{
		case RANDOM_NO_REPLACE:
			train_randftr_str =
				new CRF_InFtrStream_RandPresent(this->debug, this->dbgname, this->dbgname,*train_ftr_str,RANDOM_NO_REPLACE,this->rseed);

			//changed by Ryan
//			train_winftr_str =
//        		new QN_InFtrStream_SeqWindow(this->debug, this->dbgname,
//                                     *train_randftr_str, this->window_len,
//                                      this->window_offset, bot_margin);
			// modified by Ryan, for context features
//			train_winftr_str =
//				new CRF_InFtrStream_SeqMultiWindow(this->debug, this->dbgname,
//                        *train_randftr_str, this->window_len,
//                         this->window_offset, bot_margin);
			train_winftr_str =
				new CRF_InFtrStream_SeqMultiWindow(this->debug, this->dbgname,
						*train_randftr_str, this->window_len,
						 this->window_offset, bot_margin,
						 this->left_context_len, this->right_context_len,
						 this->extract_segment_features, this->use_boundary_delta_ftrs);
			break;
		case RANDOM_REPLACE:
			train_randftr_str =
				new CRF_InFtrStream_RandPresent(this->debug, this->dbgname, this->dbgname,*train_ftr_str,RANDOM_REPLACE,this->rseed);

			//changed by Ryan
//			train_winftr_str =
//        		new QN_InFtrStream_SeqWindow(this->debug, this->dbgname,
//                                     *train_randftr_str, this->window_len,
//                                      this->window_offset, bot_margin);
			// modified by Ryan, for context features
//			train_winftr_str =
//				new CRF_InFtrStream_SeqMultiWindow(this->debug, this->dbgname,
//						*train_randftr_str, this->window_len,
//						 this->window_offset, bot_margin);
			train_winftr_str =
				new CRF_InFtrStream_SeqMultiWindow(this->debug, this->dbgname,
						*train_randftr_str, this->window_len,
						 this->window_offset, bot_margin,
						 this->left_context_len, this->right_context_len,
						 this->extract_segment_features, this->use_boundary_delta_ftrs);
			break;
		case SEQUENTIAL:
			//changed by Ryan
//   			train_winftr_str =
//   				new QN_InFtrStream_SeqWindow(this->debug, this->dbgname,
//   						*train_ftr_str, this->window_len,
//   						this->window_offset, bot_margin);
			// modified by Ryan, for context features
//			train_winftr_str =
//				new CRF_InFtrStream_SeqMultiWindow(this->debug, this->dbgname,
//						*train_ftr_str, this->window_len,
//						 this->window_offset, bot_margin);
			train_winftr_str =
				new CRF_InFtrStream_SeqMultiWindow(this->debug, this->dbgname,
						*train_ftr_str, this->window_len,
						 this->window_offset, bot_margin,
						 this->left_context_len, this->right_context_len,
						 this->extract_segment_features, this->use_boundary_delta_ftrs);
   			break;
   		default:
   			cerr << "Invalid training sequence type! ABORT!" << endl;
   			exit(-1);
	}

	//changed by Ryan
//    QN_InFtrStream_SeqWindow* cv_winftr_str;
	CRF_InFtrStream_SeqMultiWindow* cv_winftr_str;

	if (cv_ftr_str != NULL ) {
	//changed by Ryan
//    cv_winftr_str =
//        new QN_InFtrStream_SeqWindow(this->debug, this->dbgname,
//                                      *cv_ftr_str, this->window_len,
//	                                  this->window_offset, bot_margin
//                                      );
		// modified by Ryan, for context features
//		cv_winftr_str =
//				new CRF_InFtrStream_SeqMultiWindow(this->debug, this->dbgname,
//										*cv_ftr_str, this->window_len,
//										 this->window_offset, bot_margin);
		cv_winftr_str =
				new CRF_InFtrStream_SeqMultiWindow(this->debug, this->dbgname,
										*cv_ftr_str, this->window_len,
										 this->window_offset, bot_margin,
										 this->left_context_len, this->right_context_len,
										 this->extract_segment_features, this->use_boundary_delta_ftrs);
	}
	else {
		cv_winftr_str=NULL;
	}

//    this->trn_stream = train_winftr_str;
 //   this->cv_stream = cv_winftr_str;


    // Create Label Streams

	// Changed by Ryan
//    QN_InLabStream_SeqWindow* cv_winlab_str=NULL;
//    QN_InLabStream_SeqWindow* train_winlab_str=NULL;
	CRF_InLabStream_SeqMultiWindow* cv_winlab_str=NULL;
	CRF_InLabStream_SeqMultiWindow* train_winlab_str=NULL;

    if (this->hardtarget_filename != 0) {
    	FILE* hardtarget_fp=NULL;
		enum { LABFILE_BUF_SIZE = 0x8000 };
		hardtarget_fp = QN_open(this->hardtarget_filename, "r", LABFILE_BUF_SIZE,"hardtarget_file");


    	QN_InLabStream_ILab* lab_str = new QN_InLabStream_ILab(this->debug,this->dbgname,
    														hardtarget_fp, 1);



   		// Create training and cross-validation streams.
    	QN_InLabStream* train_lab_str = NULL;
    	QN_InLabStream* cv_lab_str = NULL;

    	if (cv_winftr_str != NULL) {
    		// Using range strings
    		QN_InLabStream_CutRange* fwd_lab_str
            	= new QN_InLabStream_CutRange(this->debug, this->dbgname, *lab_str,
                                          this->train_sent_range, this->cv_sent_range);
			QN_InLabStream_Cut* train_lab_str_cut = (QN_InLabStream_Cut*)fwd_lab_str;
			cv_lab_str = (QN_InLabStream*) new QN_InLabStream_Cut2(*train_lab_str_cut);
    		train_lab_str = (QN_InLabStream*)train_lab_str_cut;
    	}
    	else {
			QN_InLabStream_CutRange* fwd_lab_str
            	= new QN_InLabStream_CutRange(this->debug, this->dbgname, *lab_str,
                                          this->train_sent_range, 0);
			QN_InLabStream_Cut* train_lab_str_cut = (QN_InLabStream_Cut*)fwd_lab_str;
    		train_lab_str = (QN_InLabStream*)train_lab_str_cut;
       		cv_lab_str=NULL;
    	}

	    // Create training and CV windows.
	    CRF_InLabStream_RandPresent* train_randlab_str =NULL;

	    // Changed by Ryan
//	    const size_t window_len = 1;
//    	bot_margin = this->window_extent - this->hardtarget_window_offset - window_len;
		bot_margin = this->window_extent - this->hardtarget_window_offset - this->window_len;

	    // Added by Ryan
//	    if (this->hardtarget_window_offset + window_len > this->window_extent)
//		{
//			string errstr="CRF_FeatureStreamManager::create() caught exception: this->hardtarget_window_offset + window_len > this->window_extent.";
//			throw runtime_error(errstr);
//		}
	    if (this->hardtarget_window_offset + this->window_len > this->window_extent)
		{
			string errstr="CRF_FeatureStreamManager::create() caught exception: this->hardtarget_window_offset + this->window_len > this->window_extent.";
			throw runtime_error(errstr);
		}


		switch ( this->train_seq_type )
		{
			case RANDOM_NO_REPLACE:
				train_randlab_str =
					new CRF_InLabStream_RandPresent(this->debug, this->dbgname, this->dbgname,*train_lab_str,RANDOM_NO_REPLACE, this->rseed );
				//train_winlab_str =
        		//	new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
                //       	              *train_randlab_str, window_len,
                //            	          this->window_offset, bot_margin);

				// Changed by Ryan
//				train_winlab_str =
//        			new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
//                        	              *train_randlab_str, window_len,
//                            	          this->hardtarget_window_offset, bot_margin);
				train_winlab_str =
					new CRF_InLabStream_SeqMultiWindow(this->debug, this->dbgname,
										  *train_randlab_str, this->window_len,
										  this->hardtarget_window_offset, bot_margin);
				break;
			case RANDOM_REPLACE:
				train_randlab_str =
					new CRF_InLabStream_RandPresent(this->debug, this->dbgname, this->dbgname,*train_lab_str,RANDOM_REPLACE, this->rseed );
				//train_winlab_str =
        		//	new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
                //        	              *train_randlab_str, window_len,
                //            	          this->window_offset, bot_margin);

				// Changed by Ryan
//				train_winlab_str =
//        			new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
//                        	              *train_randlab_str, window_len,
//                            	          this->hardtarget_window_offset, bot_margin);
				train_winlab_str =
					new CRF_InLabStream_SeqMultiWindow(this->debug, this->dbgname,
										  *train_randlab_str, this->window_len,
										  this->hardtarget_window_offset, bot_margin);
				break;
			case SEQUENTIAL:
				//train_winlab_str =
				//	new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
				//							*train_lab_str, window_len,
				//							this->window_offset, bot_margin);

				// Changed by Ryan
//				train_winlab_str =
//					new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
//											*train_lab_str, window_len,
//											this->hardtarget_window_offset, bot_margin);
				train_winlab_str =
					new CRF_InLabStream_SeqMultiWindow(this->debug, this->dbgname,
										  *train_lab_str, this->window_len,
										  this->hardtarget_window_offset, bot_margin);
				break;
   			default:
   				cerr << "Invalid training sequence type! ABORT!" << endl;
   				exit(-1);
		}

		if (cv_lab_str != NULL ) {
//    		cv_winlab_str =
//        		new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
//           	                          *cv_lab_str, window_len,
//                	                      this->window_offset, bot_margin
//                    	                  );

			// Changed by Ryan
//    		cv_winlab_str =
//        		new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
//           	                          *cv_lab_str, window_len,
//                	                      this->hardtarget_window_offset, bot_margin
//                    	                  );
    		cv_winlab_str =
				new CRF_InLabStream_SeqMultiWindow(this->debug, this->dbgname,
									  *cv_lab_str, this->window_len,
									  this->hardtarget_window_offset, bot_margin);

		}
		else {
			cv_winlab_str=NULL;
		}
    }
	else {
		cv_winlab_str=NULL;
		train_winlab_str=NULL;
	}

	this->trn_stream = new CRF_FeatureStream(train_winftr_str, train_winlab_str,this->debug);
	if (cv_winftr_str != NULL) {
		this->cv_stream = new CRF_FeatureStream(cv_winftr_str, cv_winlab_str,this->debug);
	}
	else {
		this->cv_stream = NULL;
	}

	// now initialize children, if any
	if (nthreads>1) {
		QNUInt32 nseg_per_child=nseg/nthreads;
		//cout << "nseg=" << nseg << " nseg_per_child=" << nseg_per_child << endl;
		for(size_t i=0;i<nthreads;i++) {
			if (this->normfile) {
				if(fseek(this->normfile,0L,SEEK_SET)<0) {
					cerr << "Can't rewind normfile in setting up threads" << endl;
					exit(1);
				}
			}
			// modified by Ryan, for context features
//			children[i]=new CRF_FeatureStreamManager(this->debug, this->dbgname, this->filename,
//													 this->format, this->hardtarget_filename,
//													 this->hardtarget_window_offset,
//													 this->width, this->first_ftr, this->num_ftrs,
//													 this->window_extent,this->window_offset, this->window_len,
//													 this->delta_order, this->delta_win, this->train_sent_range,
//													 this->cv_sent_range,
//													 this->normfile, // WARNING needs rewinding
//													 this->norm_mode,this->norm_am, this->norm_av,
//													 this->train_seq_type, this->rseed, 1);
			children[i]=new CRF_FeatureStreamManager(this->debug, this->dbgname, this->filename,
													 this->format, this->hardtarget_filename,
													 this->hardtarget_window_offset,
													 this->width, this->first_ftr, this->num_ftrs,
													 this->window_extent,this->window_offset, this->window_len,
													 this->left_context_len,this->right_context_len,this->extract_segment_features,
													 this->use_boundary_delta_ftrs,
													 this->delta_order, this->delta_win, this->train_sent_range,
													 this->cv_sent_range,
													 this->normfile, // WARNING needs rewinding
													 this->norm_mode,this->norm_am, this->norm_av,
													 this->train_seq_type, this->rseed, 1);
			children[i]->childnum=i;

			// the stuff below doesn't work because the streams aren't indexed
			// need to modify crf_featurestream instead


			// make each child cover a portion of the stream
			//cout << "stream " << i << " view " << (i*nseg_per_child) << " " << ((i==nthreads-1)?QN_ALL:nseg_per_child) << endl;
			children[i]->trn_stream->view(i*nseg_per_child,(i==nthreads-1)?QN_ALL:nseg_per_child);

		}

	}

    QN_OUTPUT("End of Stream Creation");
}

/*
 * CRF_FeatureStreamManager::join
 *
 * Input: instr - input feature stream
 *
 * Alters the current streams managed to be concatenated streams with the initial stream concatenated with
 * the new stream.  Manages the concatenation of multi-threaded streams as well as training and cross-validation
 * streams.
 */
void CRF_FeatureStreamManager::join(CRF_FeatureStreamManager* instr) {
	this->old_trn_stream=this->trn_stream;
	this->trn_stream = this->trn_stream->join(instr->trn_stream);
	//cout << "Joined numsegs: " << this->trn_stream->num_segs();
	cout << "Old trn_stream id: " << this->old_trn_stream << endl;
	cout << "Joined trn_stream id: " << this->trn_stream << endl;
	if (this->cv_stream != NULL) {
		this->cv_stream = this->cv_stream->join(instr->cv_stream);
	}
	if (this->nthreads>1 || instr->nthreads>1) {
		if (this->nthreads != instr->nthreads) {
			cerr << "Tried to join streams with different number of threads! " << this->nthreads << "!=" << instr->nthreads << endl;
			exit(1);
		}
		for(size_t i=0;i<this->nthreads;i++) {
			this->children[i]->trn_stream = this->children[i]->trn_stream->join(instr->children[i]->trn_stream);
			if (this->children[i]->cv_stream != NULL) {
				this->children[i]->cv_stream = this->children[i]->cv_stream->join(instr->children[i]->cv_stream);
			}
		}
	}
	//this->display();
}

/*
 * CRF_FeatureStreamManager::join
 *
 * Displays the streams being managed by the Manager to standard output.
 * Debugging function.
 */
void CRF_FeatureStreamManager::display()
{
//	this->rewind(); // Reset it to the beginning
	cout << "trn_stream id: " << this->trn_stream << endl;
	cout << "Old Partition" << endl;
	//this->old_trn_stream->display();
	cout << "Training Partition" << endl;
	this->trn_stream->display();
	cout << "CV Partition" << endl;
	if (this->cv_stream != NULL) {
		this->cv_stream->display();
	}
//	this->rewind(); // And do it again once we've displayed the files
}


/*
 * CRF_FeatureStreamManager::setFiles
 *
 * Input: see the constructor above
 *
 * Mutator function to set file information after the object has been instantiated.  Should be followed
 * by a call to create().
 */
void CRF_FeatureStreamManager::setFiles(char* filename, const char* format, char* hardtarget_file,
								size_t width, size_t first_ftr, size_t num_ftrs)
{
	this->filename=filename;
	this->format=format;
	this->hardtarget_filename=hardtarget_file;
	this->width=width;
	this->first_ftr=first_ftr;
	this->num_ftrs=num_ftrs;
}

/*
 * CRF_FeatureStreamManager::setNorm
 *
 * Input: see the constructor above
 *
 * Mutator function to set normalization information after the object has been instantiated.  Should be
 * followed by a call to create().
 */

void CRF_FeatureStreamManager::setNorm(FILE* normfile, int norm_mode, double norm_am, double norm_av)
{
	this->normfile=normfile;
	this->norm_mode=norm_mode;
	this->norm_am=norm_am;
	this->norm_av=norm_av;
}

/*
 * CRF_FeatureStreamManager::setRanges
 *
 * Input: see the constructor above
 *
 * Mutator function to set range information after the object has been instantiated.  Should be
 * followed by a call to create().
 */

void CRF_FeatureStreamManager::setRanges(char* train_range, char* cv_range,
 							     size_t train_cache_frames, int train_cache_seed)
{
	this->train_sent_range=train_range;
	this->cv_sent_range=cv_range;
	this->train_cache_frames=train_cache_frames;
	this->train_cache_seed=train_cache_seed;
}

/*
 * CRF_FeatureStreamManager::setWindow
 *
 * Input: see the constructor above
 *
 * Mutator function to set windowing information after the object has been instantiated.  Should be
 * followed by a call to create().
 */

void CRF_FeatureStreamManager::setWindow(size_t window_extent, size_t window_offset,
									size_t window_len)
{
	this->window_extent=window_extent;
	this->window_offset=window_offset;
	this->window_len=window_len;
}

/*
 * CRF_FeatureStreamManager::setDeltas
 *
 * Input: see the constructor above
 *
 * Mutator function to set delta information after the object has been instantiated.  Should be
 * followed by a call to create().
 */
void CRF_FeatureStreamManager::setDeltas(int delta_order, int delta_win)
{
	this->delta_order=delta_order;
	this->delta_win=delta_win;
}

/*
 * CRF_FeatureStreamManager::setDebug
 *
 * Input: see the constructor above
 *
 * Mutator function to set debug information after the object has been instantiated.
 */
int CRF_FeatureStreamManager::setDebug(int dbgLvl, const char* dbg) {
	this->debug=dbgLvl;
	this->dbgname=dbg;
	return this->debug;
}

/*
 * CRF_FeatureStreamManager::getNumFtrs
 *
 *
 * Accessor function to get the number of features.
 */
size_t CRF_FeatureStreamManager::getNumFtrs() {
	return this->trn_stream->num_ftrs();
}

/*
 * CRF_FeatureStreamManager::setUtt
 *
 * Input: utt - segment to set the stream position to
 *
 * Mutator function to set the training stream to a particular segment.
 * used in debugging.
 */
void CRF_FeatureStreamManager::setUtt(QNUInt32 utt) {
	this->trn_stream->set_pos(utt,0);
}

