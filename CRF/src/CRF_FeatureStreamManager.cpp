#include "CRF_FeatureStreamManager.h"

// ERROR & PROCESS LOGGING NEEDS TO BE ADDED TO THIS CLASS

CRF_FeatureStreamManager::CRF_FeatureStreamManager(int dbg, const char* dname,
									char* fname, const char* fmt, char* ht_fname, size_t ht_offset,
									size_t width, size_t first_ftr, size_t num_ftrs,
									size_t win_ext, size_t win_off, size_t win_len,
									int delta_o, int delta_w,
									char* trn_rng, char* cv_rng,
									FILE* nfile, int n_mode, double n_am, double n_av, seqtype ts, QNUInt32 rseed)
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
	 delta_order(delta_o),
	 delta_win(delta_w),
	 train_sent_range(trn_rng),
	 cv_sent_range(cv_rng),
	 normfile(nfile),
	 norm_mode(n_mode),
	 norm_am(n_am),
	 norm_av(n_av),
	 train_seq_type(ts),
	 rseed(rseed)
{
	this->create();
}


CRF_FeatureStreamManager::~CRF_FeatureStreamManager()
{
}

// Modify this so that we can:
//  * Create feature streams with NO cv stream
//  * Create feature streams with NO corresponding label stream
//  * Create feature streams with one of three types of selection (Sequential, Random, Random no Replace)

void CRF_FeatureStreamManager::create()
{
	QN_InFtrStream* ftr_str = NULL;     // Temporary stream holder.
    int index = 1;                      // training always requires indexed
    int buffer_frames = 500;

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
    
    // Create training and CV windows.
    size_t bot_margin = this->window_extent - this->window_offset - this->window_len;
	QN_InFtrStream_SeqWindow* train_winftr_str;
	CRF_InFtrStream_RandPresent* train_randftr_str = NULL;
	
	switch ( this->train_seq_type )
	{
		case RANDOM_NO_REPLACE:
			train_randftr_str =
				new CRF_InFtrStream_RandPresent(this->debug, this->dbgname, this->dbgname,*train_ftr_str,RANDOM_NO_REPLACE,this->rseed);
			train_winftr_str =
        		new QN_InFtrStream_SeqWindow(this->debug, this->dbgname,
                                     *train_randftr_str, this->window_len,
                                      this->window_offset, bot_margin);
			break;
		case RANDOM_REPLACE:
			train_randftr_str =
				new CRF_InFtrStream_RandPresent(this->debug, this->dbgname, this->dbgname,*train_ftr_str,RANDOM_REPLACE,this->rseed);
			train_winftr_str =
        		new QN_InFtrStream_SeqWindow(this->debug, this->dbgname,
                                     *train_randftr_str, this->window_len,
                                      this->window_offset, bot_margin);
			break;
		case SEQUENTIAL:
   			train_winftr_str =
   				new QN_InFtrStream_SeqWindow(this->debug, this->dbgname,
   						*train_ftr_str, this->window_len,
   						this->window_offset, bot_margin);
   			break;
   		default:
   			cerr << "Invalid training sequence type! ABORT!" << endl;
   			exit(-1);
	}            

   QN_InFtrStream_SeqWindow* cv_winftr_str;
                                      
	if (cv_ftr_str != NULL ) {
    cv_winftr_str =
        new QN_InFtrStream_SeqWindow(this->debug, this->dbgname,
                                      *cv_ftr_str, this->window_len,
	                                  this->window_offset, bot_margin
                                      );
	}
	else {
		cv_winftr_str=NULL;
	}

//    this->trn_stream = train_winftr_str;
 //   this->cv_stream = cv_winftr_str;

    
    // Create Label Streams
    QN_InLabStream_SeqWindow* cv_winlab_str=NULL;
    QN_InLabStream_SeqWindow* train_winlab_str=NULL;
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
	    const size_t window_len = 1;
    	bot_margin = window_extent - this->hardtarget_window_offset - window_len;

		switch ( this->train_seq_type )
		{
			case RANDOM_NO_REPLACE:
				train_randlab_str =
					new CRF_InLabStream_RandPresent(this->debug, this->dbgname, this->dbgname,*train_lab_str,RANDOM_NO_REPLACE, this->rseed );
				//train_winlab_str =
        		//	new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
                //       	              *train_randlab_str, window_len,
                //            	          this->window_offset, bot_margin);
				train_winlab_str =
        			new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
                        	              *train_randlab_str, window_len,
                            	          this->hardtarget_window_offset, bot_margin);
				break;
			case RANDOM_REPLACE:
				train_randlab_str =
					new CRF_InLabStream_RandPresent(this->debug, this->dbgname, this->dbgname,*train_lab_str,RANDOM_REPLACE, this->rseed );
				//train_winlab_str =
        		//	new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
                //        	              *train_randlab_str, window_len,
                //            	          this->window_offset, bot_margin);
				train_winlab_str =
        			new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
                        	              *train_randlab_str, window_len,
                            	          this->hardtarget_window_offset, bot_margin);
				break;
			case SEQUENTIAL:
				//train_winlab_str =
				//	new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
				//							*train_lab_str, window_len,
				//							this->window_offset, bot_margin);
				train_winlab_str =
					new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
											*train_lab_str, window_len,
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
    		cv_winlab_str =
        		new QN_InLabStream_SeqWindow(this->debug, this->dbgname,
           	                          *cv_lab_str, window_len,
                	                      this->hardtarget_window_offset, bot_margin
                    	                  );

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
    QN_OUTPUT("End of Stream Creation");
}

void CRF_FeatureStreamManager::join(CRF_FeatureStreamManager* instr) {
	this->trn_stream = this->trn_stream->join(instr->trn_stream);
	if (this->cv_stream != NULL) {
		this->cv_stream = this->cv_stream->join(instr->cv_stream);
	}
}


void CRF_FeatureStreamManager::display()
{
//	this->rewind(); // Reset it to the beginning
	cout << "Training Partition" << endl;
	this->trn_stream->display();
	cout << "CV Partition" << endl;
	if (this->cv_stream != NULL) {
		this->cv_stream->display();
	}
//	this->rewind(); // And do it again once we've displayed the files
}



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

void CRF_FeatureStreamManager::setNorm(FILE* normfile, int norm_mode, double norm_am, double norm_av)
{
	this->normfile=normfile;
	this->norm_mode=norm_mode;
	this->norm_am=norm_am;
	this->norm_av=norm_av;
}

void CRF_FeatureStreamManager::setRanges(char* train_range, char* cv_range,
 							     size_t train_cache_frames, int train_cache_seed)
{
	this->train_sent_range=train_range;
	this->cv_sent_range=cv_range;
	this->train_cache_frames=train_cache_frames;
	this->train_cache_seed=train_cache_seed;
}

void CRF_FeatureStreamManager::setWindow(size_t window_extent, size_t window_offset,
									size_t window_len)
{
	this->window_extent=window_extent;
	this->window_offset=window_offset;
	this->window_len=window_len;
}

void CRF_FeatureStreamManager::setDeltas(int delta_order, int delta_win)
{
	this->delta_order=delta_order;
	this->delta_win=delta_win;	
}


int CRF_FeatureStreamManager::setDebug(int dbgLvl, const char* dbg) {
	this->debug=dbgLvl;
	this->dbgname=dbg;
	return this->debug;
}

size_t CRF_FeatureStreamManager::getNumFtrs() {
	return this->trn_stream->num_ftrs();
}

void CRF_FeatureStreamManager::setUtt(QNUInt32 utt) {
	this->trn_stream->set_pos(utt,0);
}

