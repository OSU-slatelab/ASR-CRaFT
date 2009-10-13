#include "CRF_LocalPosteriorBuilder.h"

CRF_LocalPosteriorBuilder::CRF_LocalPosteriorBuilder(CRF_Model* crf_in)
	: CRF_StdGradBuilderLog(crf_in)
{
}

CRF_LocalPosteriorBuilder::~CRF_LocalPosteriorBuilder()
{
}

CRF_Seq* CRF_LocalPosteriorBuilder::buildFtrSeq(CRF_FeatureStream *ftr_strm)
{
	CRF_Seq* cur_seq = NULL;
	CRF_Seq* seq_head = NULL;
	CRF_Seq* next_seq = NULL;
	CRF_Seq* last_seq = NULL;
	double* alpha;
	QNUInt32 uc=0;
	QNUInt32 lambda_len = this->crf->getLambdaLen();	
	double logLi = 0.0;
	size_t bunch_size = 3;
	size_t num_ftrs=ftr_strm->num_ftrs();
	
	if (this->ftr_buf==NULL) { // First pass through, initialize the buffers
		this->ftr_buf = new float[num_ftrs*bunch_size];
		this->lab_buf = new QNUInt32[bunch_size];
	}
	
			
	size_t ftr_count;
	
	for (QNUInt32 i=0; i<this->lambda_len; i++) {
		this->ExpF[i]=0.0;
	}
			
	int seq_cnt=0;	
	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		ftr_count=ftr_strm->read(bunch_size,ftr_buf,lab_buf);
			
		for (QNUInt32 i=0; i<ftr_count; i++) {
			// Now, separate the bunch into individual frames
			float* new_buf = new float[num_ftrs];
			for (QNUInt32 j=0; j<num_ftrs; j++) {
				int idx=i*num_ftrs+j;
				new_buf[j]=ftr_buf[idx];
				//cout << " " << new_buf[j];
			}
			//cout << endl;
			// Store the current frame/label information in a sequence node
			//	* sequence nodes create a doubly-linked list, with the previous node known at creation time
			//  * new_buf will be deleted when this sequence object gets deleted
			next_seq = new CRF_Seq(new_buf,num_ftrs,lab_buf[i],this->num_labs,cur_seq);
			if (seq_head == NULL)
			{
				seq_head=next_seq;
			}
			if (cur_seq != NULL) {
				cur_seq->setNext(next_seq);
			}
			cur_seq=next_seq;
			
			//double value=this->computeTransMatrix(cur_seq,crf);
			double value=this->computeTransMatrixLog(cur_seq); // logspace computation
			double alpha_scale;
			try {
				alpha_scale=this->computeAlphaLog(cur_seq);
			}
			catch (exception &e) {
				string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq caught exception: "+string(e.what())+" while computing alpha";
				throw runtime_error(errstr);
				return(NULL);
			}
			
			seq_cnt++;				
		}
		
	} while (ftr_count >= bunch_size);
	
	alpha=cur_seq->getAlpha();
	double Zx;
	try {
		Zx=logAdd(alpha,this->num_labs);
	}
	catch (exception &e) {
		string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq caught exception: "+string(e.what())+" while accumulating Zx";
		throw runtime_error(errstr);
		return(NULL);
	}
		
	int s_counter=0;
	while (cur_seq != NULL) {
		last_seq=cur_seq;

		double beta_scale;
		try {
			beta_scale=this->computeBetaLog(cur_seq);
		}
		catch (exception &e) {
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq caught exception: "+string(e.what())+" while computing beta";
			throw runtime_error(errstr);
			return(NULL);
		}
		
		
		double* alpha_beta=cur_seq->getAlphaBeta();
		double* beta=cur_seq->getBeta();
		double* alpha=cur_seq->getAlpha();
		//double scale=cur_seq->getScale();
		
		double alpha_beta_tot=0.0;
		
		for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
			try {
				alpha_beta[clab]=alpha[clab]+beta[clab]-Zx; // logspace computation
				alpha_beta_tot += expE(alpha_beta[clab]);
			}
			catch (overflow_error &e) {
				string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq caught exception: "+string(e.what())+" while computing alpha_beta";
				throw runtime_error(errstr);
				return(NULL);
			}
		}
		if ((alpha_beta_tot >1.1))  {
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq: Probability sums greater than 1.0";
			throw runtime_error(errstr);
			return(NULL);
		}
		if (alpha_beta_tot < 0.9) {
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq: Probability sums less than 1.0";
			throw runtime_error(errstr);
			return(NULL);
		} 
		cur_seq=cur_seq->getPrev();
		s_counter++;
	}
	
	return seq_head;	
}
