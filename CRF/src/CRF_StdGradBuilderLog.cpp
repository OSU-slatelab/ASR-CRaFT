#include "CRF_StdGradBuilderLog.h"

CRF_StdGradBuilderLog::CRF_StdGradBuilderLog(CRF_Model* crf_in)
	: CRF_StdGradBuilder(crf_in)
{
	for (QNUInt32 i=0; i<this->num_labs; i++) {
		this->alpha_base[i]=0.0;
	}
	this->logAddAcc = new double[this->num_labs];
}

//CRF_StdGradBuilderLog::~CRF_StdGradBuilderLog()
//{
//}



double CRF_StdGradBuilderLog::buildGradient(CRF_FeatureStream *ftr_strm, double* grad)
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
				string errstr="CRF_StdGradBuilderLog::buildGradient caught exception: "+string(e.what());
				throw runtime_error(errstr);
				return(-1);
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
		string errstr="CRF_StdGradBuilderLog::buildGradient caught exception: "+string(e.what())+" while accumulating Zx";
		throw runtime_error(errstr);
		return(-1);
	}
	
	int s_counter=0;
	while (cur_seq != NULL) {
		last_seq=cur_seq;


		double beta_scale;
		try {
			beta_scale=this->computeBetaLog(cur_seq);
		}
		catch (exception &e) {
			string errstr="CRF_StdGradBuilderLog::buildGradient caught exception: "+string(e.what());
			throw runtime_error(errstr);
			return(-1);
		}
		
		logLi += this->computeExpFLog(cur_seq,grad,Zx);
		
		cur_seq=cur_seq->getPrev();
		s_counter++;
	}
	for (QNUInt32 i=0; i<lambda_len; i++) {
		//ExpF[i]=ExpF[i]/Zx; // Removed - Zx included in alpha_beta above now
		grad[i]-=this->ExpF[i];
	}
	logLi -=Zx;

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

double CRF_StdGradBuilderLog::computeAlphaLog(CRF_Seq* cur_seq) 
{	
	double alpha_tot=0.0;
	double* alpha;
	double* Mtrans=cur_seq->getMtrans();
	double* Mstate=cur_seq->getMstate();
	
	if (cur_seq->getPrev() == NULL) {
		alpha = this->alpha_base; // alpha_tmp = 1
	}
	else {
		alpha=cur_seq->getPrev()->getAlpha(); // alpha_tmp = alpha[i-1]
	}
	double* new_alpha=cur_seq->getAlpha();  // alpha[i] = 0
				
	// cblas_sgemm(Row|Column, TransA?, TransB?, M, N, K, scale_A, A, lda, B, ldb, scale_B, C, ldc)
	// Where M is the rows of A (if TransA?==no)
	//       N is the columns of B (if TransB?==no)
	//       K is the columns of A and the rows of B (if TransA?==no && TransB? == no)
	//       lda is the linear depth of A (how many elements you must stride to get to the end of a row
	//       ldb is the linear depth of B
	//       ldc is the linear depth of C
			
	//cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, nLabs, nLabs, 1.0f, alpha, nLabs, M, nLabs, 1.0f, new_alpha,nLabs);
			
	// Replace the above with a computation in logspace
						
	/*double tmp=0;
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		new_alpha[clab]=alpha[0]+M[0+clab];  // Ensure that we start with a non-empty value for new_alpha
		for (QNUInt32 plab=1; plab<nLabs; plab++) {
			//tmp=alpha[plab]*M[plab*nLabs+clab];
			//new_alpha[clab]=new_alpha[clab]+tmp;
			tmp=alpha[plab]+M[plab*nLabs+clab]; // logspace compuation
			new_alpha[clab]=logAdd(new_alpha[clab],tmp); // logspace computation
		}
		//new_alpha[clab]*=R[clab];
		new_alpha[clab]+=R[clab]; //logspace computation
		//cout << " Setting new_alpha[i] to " << new_alpha[i] << endl;
	}
	//cout << "Ending alpha computation" << endl;
	 */
			 
	 // This loop is a transformation of the above that makes use of the logAdd function applied to
	 // an array.
			 
	 for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
	 	this->logAddAcc[0]=alpha[0]+Mtrans[0+clab];
	 	double maxv=this->logAddAcc[0];
	 	for (QNUInt32 plab=1; plab<this->num_labs; plab++) {
	 		this->logAddAcc[plab]=alpha[plab]+Mtrans[plab*this->num_labs+clab];
	 		if (this->logAddAcc[plab]>maxv) {
	 			maxv=logAddAcc[plab];
	 		}
	 	}
	 	try {
	 		new_alpha[clab]=logAdd(this->logAddAcc,maxv,this->num_labs);
	 	}
	 	catch (exception &e) {
	 		string errstr="CRF_StdGradBuilderLog::computeAlphaLog caught exception: "+string(e.what())+" while computing alpha";
			throw runtime_error(errstr);
			return(-1);
	 	}
	 	new_alpha[clab]+=Mstate[clab];
	 	//alpha_tot+=new_alpha[clab];
	 }
	
	// Scaling removed in logspace computation for the moment
	//double scale=0.0;
	/*double scale=new_alpha[0];
	for (QNUInt32 i=1; i<nLabs; i++) {
		//scale=scale+new_alpha[i];
		//cout << " ++Scale: " << scale << endl;
		scale=logAdd(scale,new_alpha[i]); //logspace computation
	}*/
	//cout << "Scale: " << scale << endl;
	//if ( (scale<1.0) && (scale>-1.0)) scale=1.0; //Don't inflate scores
	/*if (scale <0.0) scale=0.0; // logspace computation - don't inflate scores
	scale=0.0; // take scaling out of the problem entirely.
	cur_seq->setScale(scale);*/
	/*for (QNUInt32 i=0; i<nLabs; i++) {
		//new_alpha[i]=new_alpha[i]/scale;
		new_alpha[i]=new_alpha[i]-scale;
		//cout << "new_alpha[" << i << "]: " << new_alpha[i] << endl;
	}*/
			
	/*double alpha_tot=new_alpha[0];
	for (QNUInt32 i=1; i<nLabs; i++) {
		//alpha_tot+=new_alpha[i];
		alpha_tot=logAdd(alpha_tot,new_alpha[i]);
	}
			
	alpha_tot=expE(alpha_tot); //logspace computation
	if ((alpha_tot > 1.000005) || (alpha_tot < 0.999)) {
		cout << "ERROR: alpha does not total to 1.0 after scaling " << alpha_tot << endl;
		exit(-1);
	}*/
	
	return alpha_tot;
}

double CRF_StdGradBuilderLog::computeBetaLog(CRF_Seq* cur_seq)
{
	double beta_tot=0.0;
	double* Mtrans;
	double* Mstate;
	double* new_beta;
	//double scale=cur_seq->getScale();
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]
	//	* Compute ExpF and logLi using alpha, beta, M_i and current features
	//	* Destroy the sequence node once it has been used
		
	if (cur_seq->getNext() == NULL) {  // Set beta for the last frame to all 1s
		new_beta = cur_seq->getBeta();
		for (QNUInt32 i=0; i<this->num_labs; i++) {
			//beta[i]=1.0/scale;
			//beta[i]=-1*scale;  //logspace computation
			new_beta[i]=0.0; // Remove scaling
		}
	}
	else {
		Mtrans=cur_seq->getNext()->getMtrans();
		Mstate=cur_seq->getNext()->getMstate();
		
		double* beta=cur_seq->getNext()->getBeta();
		new_beta=cur_seq->getBeta();
			
		for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
			//tmp_beta[clab]=(beta[clab]*R[clab])/scale;
			//tmp_beta[clab]=beta[clab]+R[clab]-scale; //logspace computation
			this->tmp_beta[clab]=beta[clab]+Mstate[clab]; //logspace computation
		}
		//cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLabs, 1, nLabs, 1.0f, M, nLabs, tmp_beta, nLabs, 1.0f, new_beta,1);
		// Transform above matrix multiply into logspace
		/*double tmp=0;			
		for (QNUInt32 plab=0; plab<nLabs; plab++) {
			new_beta[plab]=M[plab*nLabs+0]+tmp_beta[0]; //logspace computation
			for (QNUInt32 clab=1; clab<nLabs; clab++) {
				//tmp=M[plab*nLabs+clab]*tmp_beta[clab];
				//new_beta[plab]+=tmp;
				tmp=M[plab*nLabs+clab]+tmp_beta[clab];
				new_beta[plab]=logAdd(new_beta[plab],tmp);
			}
		}*/
			
		for (QNUInt32 plab=0; plab<this->num_labs; plab++) {
			this->logAddAcc[0]=Mtrans[plab*this->num_labs+0]+this->tmp_beta[0];
		 	double maxv=this->logAddAcc[0];
			for (QNUInt32 clab=1; clab<this->num_labs; clab++) {
				this->logAddAcc[clab]=Mtrans[plab*this->num_labs+clab]+this->tmp_beta[clab];
				if (this->logAddAcc[clab]>maxv) {
					maxv=this->logAddAcc[clab];
				}
			}
			try {
				new_beta[plab]=logAdd(this->logAddAcc,maxv,this->num_labs);
			}
			catch (exception &e) {
	 			string errstr="CRF_StdGradBuilderLog::computeBetaLog caught exception: "+string(e.what())+" while computing beta";
				throw runtime_error(errstr);
				return(-1);				
			}
		}
			
	}
	return beta_tot;
}

double CRF_StdGradBuilderLog::computeExpFLog(CRF_Seq* cur_seq, double* grad, double Zx)
{
	//vector<QNUInt32>::iterator featureIter;
	double logLi=0.0;
	double alpha_beta;
	QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;

	double* lambda = this->crf->getLambda();	
	double* Mtrans=cur_seq->getMtrans();
	double* Mstate=cur_seq->getMstate();
	float* ftr_buff = cur_seq->getFtr();
	QNUInt32 cur_lab=cur_seq->getLab();
	double* cur_alpha=cur_seq->getAlpha();
	double* new_beta=cur_seq->getBeta();
	QNUInt32 prev_lab;
	double* prev_alpha;
	if (cur_seq->getPrev() != NULL) {
		prev_lab=cur_seq->getPrev()->getLab();
		prev_alpha=cur_seq->getPrev()->getAlpha();
		//prev_scale=cur_seq->getPrev()->getScale();
	}
	else {
		//prev_lab=-1;
		prev_lab=this->num_labs+1; // Unsigned means we can't set it to -1, set it to max label plus 1 instead
		prev_alpha=this->alpha_base;
		//prev_scale=nLabs;
		//prev_scale=log(nLabs);  //logspace computation
		//prev_scale=0.0; // Remove scaling from consideration
	}
			
	for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
		try {
			//alpha_beta=cur_alpha[clab]*new_beta[clab]*scale;  // Move the Zx correction out
			//alpha_beta=expE(cur_alpha[clab]+new_beta[clab]+scale-Zx); // logspace computation
			alpha_beta=expE(cur_alpha[clab]+new_beta[clab]-Zx); // logspace computation
			alpha_beta_tot += alpha_beta;
			//cout << "Ending alpha_beta computation " << alpha_beta << " " << alpha_beta_tot << endl;
		}
		catch (overflow_error &e) {
			cerr << "Error when computing alpha_beta for current label " << clab;
			cerr << " Value: " << cur_alpha[clab]+new_beta[clab]-Zx << endl;
			cerr << " Alpha: " << cur_alpha[clab] << " Beta: " << new_beta[clab] << endl;
			cerr << " Zx: " << Zx << endl;
			cerr << e.what() << endl;
			exit(-1);
		}
		//CRF_RangeNode* rngHead = this->crf->getRangeMachine()->getRange(clab);
		//CRF_RangeNode* rngObj = rngHead;
		//while (rngObj != NULL) {
			//for (QNUInt32 k=rngObj->getStart(); k<=rngObj->getEnd(); k=k+rngObj->getIncrement()) {
		//vector<QNUInt32> rngVec = this->crf->getRangeMachine()->getRange(clab);
		//for (featureIter = rngVec.begin(); featureIter != rngVec.end(); featureIter++)
		//for (QNUInt32 k=0; k<rngVec.size(); k++) // This is the wrong way to do it, but it may be faster
		//{		
		//	this->ExpF[lc]+=alpha_beta*ftr_buff[rngVec[k]];
		//	if (clab==cur_lab) {
		//		grad[lc]+=ftr_buff[rngVec[k]]; 
		//		logLi += lambda[lc]*ftr_buff[rngVec[k]];
		//		//cout << "logLi is now " << logLi << endl;
		//	}
		//	lc++;
			//}
			//rngObj = rngObj->getNext();
		//} 
		//if (this->crf->getRangeMachine()->getBias(clab)) {
		//	this->ExpF[lc]+=alpha_beta;
		//	if (clab==cur_lab) {
		//		grad[lc]+=1;
		//		logLi += lambda[lc];
				//cout << "logLi is now " << logLi << endl;
		//	}
		//	lc++;
		//}
		bool match=(clab==cur_lab);
		logLi+=this->crf->getFeatureMap()->computeExpFState(ftr_buff,lambda,lc,this->ExpF,grad,alpha_beta,match, clab);		
		for (QNUInt32 plab=0; plab<this->num_labs; plab++) {
			//rngObj = this->crf->getRangeMachine()->getRange(plab,clab);
			//rngVec = this->crf->getRangeMachine()->getRange(plab,clab);
			
			QNUInt32 idx = plab*this->num_labs+clab;
			
			// For all computations like this, i is the current label and j is the previous label
			//	(row is previous, column is current in the M matrix)
			//alpha_beta=alpha[plab]*M[idx]*R[clab]*new_beta[clab]/Zx;
				
			try {
				//alpha_beta=alpha[plab]*M[idx]*R[clab]*new_beta[clab]; // Move Zx out
				alpha_beta=expE(prev_alpha[plab]+Mtrans[idx]+Mstate[clab]+new_beta[clab]-Zx);
				alpha_beta_trans_tot+=alpha_beta;
				//cout << "Ending alpha_beta trans computation " << alpha_beta << " " << alpha_beta_trans_tot << endl;
			}
			catch (overflow_error &e) {
				cerr << "Error when computing alpha_beta for current label " << clab << " prior label " << plab;
				cerr << " Value: " << prev_alpha[plab]+Mtrans[idx]+Mstate[clab]+new_beta[clab]-Zx << endl;
				cerr << " Alpha: " << cur_alpha[clab] << " Beta: " << new_beta[clab] << endl;
				cerr << " Zx: " << Zx << endl;
				cerr << " Total Alpha Beta: " << alpha_beta_trans_tot << endl;			
				cerr << e.what() << endl;
				exit(-1);
			}
			//cout << "Alpha Beta plab " << plab << " " << idx << " " << clab << " " << alpha_beta << endl;
			//while (rngObj != NULL) {
			//	for (QNUInt32 k=rngObj->getStart(); k<=rngObj->getEnd(); k=k+rngObj->getIncrement()) {
			//for (featureIter = rngVec.begin(); featureIter != rngVec.end(); featureIter++)
			//for (QNUInt32 k=0; k<rngVec.size(); k++) // This is the wrong way to do it, but it may be faster
			//{
			//		this->ExpF[lc]+=alpha_beta*ftr_buff[rngVec[k]];
			//		if ((clab==cur_lab) && (plab==prev_lab)) {
			//			grad[lc]+=ftr_buff[rngVec[k]];
			//			logLi += lambda[lc]*ftr_buff[rngVec[k]];
						//cout << "logLi is now " << logLi << endl;
			//		}
			//		lc++;
				//}
				//rngObj=rngObj->getNext();
			//}
			//if (this->crf->getRangeMachine()->getBias(plab,clab)) {
			//	this->ExpF[lc] += alpha_beta;
			//	if ((clab==cur_lab) && (plab==prev_lab)) {
			//		grad[lc]+=1;
			//		logLi+=lambda[lc];
					//cout << "logLi is now " << logLi << en	dl;
			//	}
			//	lc++;
			//}
			match=((clab==cur_lab)&&(plab==prev_lab));
			logLi+=this->crf->getFeatureMap()->computeExpFTrans(ftr_buff,lambda,lc,this->ExpF,grad,alpha_beta,match, plab, clab);
		}
	}
		
	//cout << " Alpha Beta Total: " << alpha_beta_tot << " Trans Total: " << alpha_beta_trans_tot << endl;
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
		cout << "ERROR: Trans Probability sums less than 1.0" << endl;
		exit(-1);			
	} 
	
	return logLi;
}
