#include "CRF_StdGradBuilder.h"

CRF_StdGradBuilder::CRF_StdGradBuilder(CRF_Model* crf_in) 
	: CRF_GradBuilder(crf_in)
{
	for (QNUInt32 i=0; i<this->num_labs; i++) {
		this->alpha_base[i]=1.0/this->num_labs;
	}
	this->lambda_len = crf_in->getLambdaLen();
	this->ExpF = new double[this->lambda_len];
}

double CRF_StdGradBuilder::buildGradient(CRF_FeatureStream* ftr_strm, double* grad)
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
	if (this->ftr_buf==NULL) {  // First pass through initialize the buffers
		this->ftr_buf = new float[num_ftrs*bunch_size];
		this->lab_buf = new QNUInt32[bunch_size];
	}
	
	size_t ftr_count;
	
	//cout << "Utt " << uc << " begin processing" << endl;
	for (QNUInt32 i=0; i<lambda_len; i++) {
		this->ExpF[i]=0.0;
	}
	
	try {
		logLi-=logE(this->num_labs);  // Correct the log liklihood scaling factor in alpha_base
	}
	catch (exception &e) {
		string errstr="CRF_StdGradBuilder::buildGradient caught exception: "+string(e.what())+" while taking log of initial scale factor";
		throw runtime_error(errstr);
		return(-1);
	}
	
	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		ftr_count=ftr_strm->read(bunch_size,this->ftr_buf,this->lab_buf);
			
		for (QNUInt32 i=0; i<ftr_count; i++) {
			//cout << "\tLabel: " << lab_buf[i] << "\tFeas:";
			// Now, separate the bunch into individual frames
			float* new_buf = new float[num_ftrs];
			for (QNUInt32 j=0; j<num_ftrs; j++) {
				int idx=i*num_ftrs+j;
				new_buf[j]=this->ftr_buf[idx];
				//cout << " " << new_buf[j];
			}
			//cout << endl;
			// Store the current frame/label information in a sequence node
			//	* sequence nodes create a doubly-linked list, with the previous node known at creation time
			//  * new_buf will be deleted when this sequence object gets deleted
			next_seq = new CRF_Seq(new_buf,num_ftrs,this->lab_buf[i],this->crf->getNLabs(),cur_seq);
			if (seq_head == NULL)
			{
				seq_head=next_seq;
			}
			if (cur_seq != NULL) {
				cur_seq->setNext(next_seq);
			}
			cur_seq=next_seq;
			double value=this->computeTransMatrix(cur_seq);
			//double scale=this->computeAlpha(cur_seq);
			double scale=this->computeAlpha(cur_seq,cur_seq->getPrev());
			
			try {	
				logLi-=logE(scale);
			}
			catch (exception &e) {
				string errstr="CRF_StdGradBuilder::buildGradient caught exception: "+string(e.what())+" while taking log of scale factor";
				throw runtime_error(errstr);
				return(-1);
			}
			
			// End of Loop:
			//	alpha[i] = alpha[i-1]*M[i]				
		}
	} while (ftr_count >= bunch_size);
		
	alpha=cur_seq->getAlpha();
	double Zx = 0.0;
	for (QNUInt32 i=0; i<this->num_labs; i++) {
		Zx+=alpha[i];
	}
	//cout << "Zx set to " << Zx << endl;

	int s_counter=0;	
	while (cur_seq != NULL) {
		last_seq=cur_seq;

		//double scale=this->computeBeta(cur_seq);
		double scale=this->computeBeta(cur_seq,cur_seq->getNext());
		logLi += this->computeExpF(cur_seq,grad,Zx);

		cur_seq=cur_seq->getPrev();
		s_counter++;
	}
		
	//cout << "Utt " << uc << " grad: ";
	for (QNUInt32 i=0; i<lambda_len; i++) {
		//ExpF[i]=ExpF[i]/Zx;
		grad[i]-=this->ExpF[i];
		//cout << grad[i] << " ";		
	}
	//cout << endl;
	//logLi-=log(Zx);
	try {
		logLi-=logE(Zx);
	}
	catch (exception &e) {
		string errstr="CRF_StdGradBuilder::buildGradient caught exception: "+string(e.what())+" while taking log of scale factor";
		throw runtime_error(errstr);
		return(-1);
	}
	
		
	//cout << "Utt " << uc << " LogLi: " << logLi << endl;

	cur_seq=seq_head;
	// A short loop to delete the nodes we've allocated
	while (cur_seq != NULL) {
		last_seq=cur_seq;
		cur_seq=cur_seq->getNext();
		delete last_seq;
	}
	uc++;
	
	return logLi;	
}

double CRF_StdGradBuilder::computeTransMatrix(CRF_Seq* seq_ptr) 
{
	double ret_val = this->computeTransMatrixLog(seq_ptr);
	QNUInt32 nLabs=this->crf->getNLabs();
	double* Mstate_ptr;
	double* Mtrans_ptr;
	for (QNUInt32 i=0; i<nLabs; i++) {
		Mstate_ptr=seq_ptr->getMstate();
		Mtrans_ptr=seq_ptr->getMtrans();
		try {
			Mstate_ptr[i]=expE(Mstate_ptr[i]);
		}
		catch (overflow_error &e) {
			cerr << "ERROR computing M_state: Overflow in cell " << i << " " << e.what() << endl;
			cout << "ERROR computing M_state: Overflow in cell " << i << " " << e.what() << endl;
			exit(-1); 
		}
		for (QNUInt32 j=0; j<nLabs; j++) {
			try {
				Mtrans_ptr[j*nLabs+i]=expE(Mtrans_ptr[j*nLabs+i]);
			}
			catch (overflow_error &e) {
				cerr << "ERROR computing M_trans: Overflow in cell " << i << "," << j << " " << e.what() << endl;
				cout << "ERROR computing M_trans: Overflow in cell " << i << "," << j << " " << e.what() << endl;
				exit(-1); 
			}
		}
	}
	return ret_val;
}

double CRF_StdGradBuilder::computeTransMatrixLog(CRF_Seq* seq_ptr)
{
	double tmp=0.0;
	
	float* ftr_buff = seq_ptr->getFtr();
	QNUInt32 nLabs = this->crf->getNLabs();
	double* Mtrans = seq_ptr->getMtrans();
	double* Mstate = seq_ptr->getMstate();
	for (QNUInt32 i=0; i<nLabs; i++) {
		Mstate[i]=LOG0;
		for (QNUInt32 j=0; j<nLabs; j++) {
			Mtrans[j*nLabs+i]=LOG0;
		}
	}
	QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf->getLambda();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		double value=0.0;
		Mstate[clab]=this->crf->getFeatureMap()->computeRi(ftr_buff,lambda,lc,clab);
		for (QNUInt32 plab=0; plab<nLabs; plab++) {
			QNUInt32 idx=plab*nLabs+clab;
			Mtrans[idx]=this->crf->getFeatureMap()->computeMij(ftr_buff,lambda,lc,plab,clab);
		}
	}
	return tmp;
}

double CRF_StdGradBuilder::computeAlpha(CRF_Seq* cur_seq, CRF_Seq* prev_seq)
{
	double* Mtrans=cur_seq->getMtrans();
	double* Mstate=cur_seq->getMstate();
	double* alpha;
	double scale=0.0;
		
	if (prev_seq == NULL) {
		alpha = this->alpha_base; // alpha_tmp = 1
	}
	else {
		alpha=prev_seq->getAlpha(); // alpha_tmp = alpha[i-1]
	}
	double* new_alpha=cur_seq->getAlpha();  // alpha[i] = 0
	
				
	// cblas_sgemm(Row|Column, TransA?, TransB?, M, N, K, scale_A, A, lda, B, ldb, scale_B, C, ldc)
	// Where M is the rows of A (if TransA?==no)
	//       N is the columns of B (if TransB?==no)
	//       K is the columns of A and the rows of B (if TransA?==no && TransB? == no)
	//       lda is the linear depth of A (how many elements you must stride to get to the end of a row
	//       ldb is the linear depth of B
	//       ldc is the linear depth of C
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, this->num_labs, this->num_labs, 1.0f, alpha, this->num_labs, Mtrans, this->num_labs, 1.0f, new_alpha,this->num_labs);
				
	for (QNUInt32 i=0; i<this->num_labs; i++) {
		new_alpha[i]=new_alpha[i]*Mstate[i];
		scale=scale+new_alpha[i];
	}
				
	if ( (scale<1.0) && (scale>-1.0)) scale=1.0; //Don't inflate scores
		cur_seq->setScale(scale);
	for (QNUInt32 i=0; i<this->num_labs; i++) {
		new_alpha[i]=new_alpha[i]/scale;
	}
	
	return scale;
}

double CRF_StdGradBuilder::computeBeta(CRF_Seq* cur_seq, CRF_Seq* next_seq)
{
	double scale=cur_seq->getScale();
	double* Mtrans;
	double* Mstate;
	double* new_beta;
		
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]
		
	if (next_seq == NULL) {
		new_beta = cur_seq->getBeta();
		for (QNUInt32 i=0; i<this->num_labs; i++) {
			new_beta[i]=1.0/scale;
			//beta[i]=1.0;
		}
	}
	else {
		Mtrans=cur_seq->getNext()->getMtrans();
		Mstate=cur_seq->getNext()->getMstate();
		
		double* beta=next_seq->getBeta();
		new_beta=cur_seq->getBeta();
			
		for (QNUInt32 i=0; i<this->num_labs; i++) {
			this->tmp_beta[i]=(beta[i]*Mstate[i])/scale;
		}
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, this->num_labs, 1, this->num_labs, 1.0f, Mtrans, this->num_labs, this->tmp_beta, this->num_labs, 1.0f, new_beta,1);
	}
	return scale;
}

double CRF_StdGradBuilder::computeExpF(CRF_Seq* cur_seq, double* grad, double Zx)
{
	double logLi=0.0;
	double alpha_beta;
	QNUInt32 cur_lab=cur_seq->getLab();
	double* cur_alpha=cur_seq->getAlpha();
	double* prev_alpha;
	double prev_scale;
	QNUInt32 prev_lab;
	if (cur_seq->getPrev() != NULL) {
		prev_lab=cur_seq->getPrev()->getLab();
		prev_alpha=cur_seq->getPrev()->getAlpha();
		prev_scale=cur_seq->getPrev()->getScale();
	}
	else {
		prev_lab=this->num_labs+1;
		prev_alpha=this->alpha_base;
		prev_scale=this->num_labs;
	}
	double* Mtrans=cur_seq->getMtrans();
	double* Mstate=cur_seq->getMstate();
	double* new_beta=cur_seq->getBeta();
	float* ftr_buff = cur_seq->getFtr();
	double scale=cur_seq->getScale();
	
	QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;
		
	for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
		alpha_beta=cur_alpha[clab]*new_beta[clab]*scale/Zx;
		alpha_beta_tot += alpha_beta;
		bool match=(clab==cur_lab);
		logLi+=this->crf->getFeatureMap()->computeExpFState(ftr_buff,lambda,lc,this->ExpF,grad,alpha_beta,match,clab);
		for (QNUInt32 plab=0; plab<this->num_labs; plab++) {
			QNUInt32 idx = plab*this->num_labs+clab;
			alpha_beta=prev_alpha[plab]*Mtrans[idx]*Mstate[clab]*new_beta[clab]/Zx;
			alpha_beta_trans_tot+=alpha_beta;
			match=((clab==cur_lab)&&(plab==prev_lab));
			logLi+=this->crf->getFeatureMap()->computeExpFTrans(ftr_buff,lambda,lc,this->ExpF,grad,alpha_beta,match,clab,plab);
		}
	}
	
	// NOTE:  The error blocks below should be reworked to throw exceptions
			
	if ((alpha_beta_tot >1.1))  {
		cout << "ERROR: Probability sums greater than 1.0 " << alpha_beta_tot << endl;
		exit(-1);
	}
	if (alpha_beta_tot < 0.9) {
		cout << "ERROR: Probability sums less than 1.0 " << alpha_beta_tot << endl;
		exit(-1);			
	} 
	if (alpha_beta_trans_tot > 1.1) {
		cout << "ERROR: Trans Probability sums greater than 1.0 " << alpha_beta_trans_tot << endl;
		exit(-1);
	}
	if (alpha_beta_trans_tot < 0.9) {
		cout << "ERROR: Trans Probability sums less than 1.0 " << alpha_beta_trans_tot << endl;
		exit(-1);			
	}
	return logLi; 
}
