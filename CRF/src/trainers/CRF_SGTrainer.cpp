#include "CRF_SGTrainer.h"

CRF_SGTrainer::CRF_SGTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str_mgr, char* wt_fname)
	: CRF_Trainer(crf_in, ftr_str_mgr, wt_fname)
{
}


void CRF_SGTrainer::train()
{
	double logLi = 0.0;
	double tmpLogLi;
	double tmp_Zx;
	double totLogLi = 0.0;
	int iCounter=0;
	int uCounter=0;
	CRF_FeatureStream *ftr_str=this->ftr_strm_mgr->trn_stream;
	//int uCounter = this->crf_ptr->getPresentations();
	ofstream ofile;


	/* This is now meaningless - we're always using logspace now, and the StateVector
	 * figures out on its own what kind of state topology to use
	 * QNUInt32 nStates = this->crf_ptr->getFeatureMap()->getNumStates();
	 * if (this->useLogspace) {
		gbuild=new CRF_NewGradBuilderLog(this->crf_ptr);
		//gbuild=new CRF_FerrGradBuilder(this->crf_ptr);
		cout << "Using Logspace training..." << endl;
	}
	else {
		gbuild=new CRF_NewGradBuilder(this->crf_ptr);
		cout << "Using Realspace training..." << endl;
	}*/
/*	switch (this->objective) {
	case EXPF :
		gbuild = new CRF_NewGradBuilder(this->crf_ptr);
		break;
	case EXPFSOFT :
		gbuild = new CRF_NewGradBuilderSoft(this->crf_ptr);
		break;
	case FERR :
		gbuild = new CRF_FerrGradBuilder(this->crf_ptr);
		break;
	default :
		gbuild = new CRF_NewGradBuilder(this->crf_ptr);
		break;
	}
	*/
	//gbuild = new CRF_NewGradBuilder(this->crf_ptr);
	//gbuild = new CRF_NewGradBuilderSoft(this->crf_ptr);
	CRF_GradBuilder* gbuild = CRF_GradBuilder::create(this->crf_ptr,this->objective);
	gbuild->setNodeList(new CRF_StateVector());
	cout << "Using Logspace training..." << endl;
	/*if (this->useLogspace) {
		if (nStates==1) {
			//gbuild=new CRF_StdGradBuilderLog(this->crf_ptr);
			gbuild=new CRF_NewGradBuilderLog(this->crf_ptr);
			gbuild->setNodeList(new CRF_StdStateVectorLog());
			cout << "Using Logspace training..." << endl;
		}
		else {
			gbuild=new CRF_NewGradBuilderLog(this->crf_ptr);
			gbuild->setNodeList(new CRF_StdNStateVectorLog());
			//cerr << "Logspace n-state training not yet implemented" << endl;
			//exit(-1);
		}
	}
	else{
		if (nStates==1) {
			//gbuild=new CRF_StdGradBuilder(this->crf_ptr);
			gbuild=new CRF_NewGradBuilder(this->crf_ptr);
			gbuild->setNodeList(new CRF_StdStateVector());
			cout << "Using normal space training..." << endl;
		}
		else {
			gbuild=new CRF_NewGradBuilder(this->crf_ptr);
			gbuild->setNodeList(new CRF_StdNStateVector());
			//cerr << "N-state training not yet implemented" << endl;
			//exit(-1);
			//gbuild=new CRF_NstateGradBuilder(this->crf_ptr);
		}
	}*/
	double* lambda=this->crf_ptr->getLambda();
	QNUInt32 lambdaLen = this->crf_ptr->getLambdaLen();
	//double* lambdaAcc=new double[lambdaLen];
	double* lambdaAcc=this->crf_ptr->getLambdaAcc();
	double* lambdaSqrAcc=new double[lambdaLen];
	double* lambdaAvg=new double[lambdaLen];
	double* lambdaVar=new double[lambdaLen];
	//int accCnt=0;
	int accCnt = this->crf_ptr->getPresentations();
	double* grad=new double[lambdaLen];
	for (QNUInt32 i=0; i<lambdaLen; i++) {
		grad[i]=0.0;
		lambdaSqrAcc[i]=0.0;
		//lambdaAcc[i]=0.0;
	}
	bool start=true;
	float invSquareVar=0.0;
	if (this->useGvar) {
		invSquareVar=1/this->gvar;
	}
	//invSquareVar=1/this->gvar;
	//if (this->useGvar) {
	//	invSquareVar=1/(this->gvar*this->gvar);
	//}
	//cout << "Iteration: " << iCounter << " starting" << endl;
	ftr_str->rewind();
	QN_SegID segid = ftr_str->nextseg();

	//QNUInt32 nSegs = this->ftr_strm->num_segs();

    time_t rawtime;
    //QNUInt32 seg,frm;

	while (iCounter<this->maxIters) {
		if (start) {
			cout << "Iteration: " << iCounter << " starting LR: " << this->lr << endl;
			start=false;
		}
		if (uCounter % this->uttRpt == 0) {
			//int pos = this->ftr_strm->get_pos(&seg,&frm);
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time << " Beginning Utt: " << uCounter << " SegID: " << segid <<  endl;
			//delete time;
		}
		tmp_Zx=uCounter%this->uttRpt; // Added for Ferr testing
		tmpLogLi=gbuild->buildGradient(ftr_str,grad,&tmp_Zx);
		//tmpLogLi=logLi;
		logLi=tmpLogLi - tmp_Zx;
		totLogLi += logLi;
		if (uCounter % this->uttRpt == 0) {
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time << " Finished Utt: " << uCounter << " logLi: " << logLi;
			cout << " Avg LogLi: " << totLogLi/(uCounter+1) << " Zx: " << tmp_Zx;
			cout << " Tmp LogLi: " << tmpLogLi << endl;
		}
		uCounter++;
		for (QNUInt32 i=0; i<lambdaLen; i++) {
			if (this->useGvar) {
				//grad[i]-=lambda[i]*invSquareVar;
				grad[i]-=grad[i]*invSquareVar;
			}
			lambda[i]+=this->lr*grad[i];
			lambdaAcc[i]+=lambda[i];
			lambdaSqrAcc[i]+=lambda[i]*lambda[i];
			grad[i]=0.0;
		}
		accCnt++;
		//cout << "Utterance " << uCounter << ": Lambda[0]: " << lambda[0] << endl;
		//cout << "Utterance " << uCounter << ": Acc Lambda[0]: " << (lambdaAcc[0]/accCnt) << endl;
		segid = ftr_str->nextseg();
		if (segid == QN_SEGID_BAD) {
			cout << "Iteration: " << iCounter << " ending" << endl;
			cout << "Iteration: " << iCounter << " applying lambda updates" << endl;

			string fname;
			stringstream ss;
			ss << this->weight_fname << ".i" << iCounter << ".out";
			ss >> fname;
			cout << "Writing Iteration " << iCounter << " weights to file " << fname << endl;
			bool chkwrite=this->crf_ptr->writeToFile(fname.c_str());
			if (!chkwrite) {
				cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
				exit(-1);
			}
			stringstream savg;
			savg << this->weight_fname << ".i" << iCounter << ".avg.out";
			savg >> fname;
			cout << "Writing Iteration " << iCounter << " average weights to file " << fname << endl;
			cout << "Writing after " << (float) accCnt << " samples" << endl;
			for (QNUInt32 i=0; i<lambdaLen; i++) {
				lambdaAvg[i]=(lambdaAcc[i]/(float)accCnt);
				lambdaVar[i]=(lambdaSqrAcc[i]/(float)accCnt)-(lambdaAvg[i]*lambdaAvg[i]);
			}
			chkwrite=this->crf_ptr->writeToFile(fname.c_str(),lambdaAvg,lambdaLen);
			if (!chkwrite) {
				cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
				exit(-1);
			}
			stringstream svar;
			svar << this->weight_fname << ".i" << iCounter << ".var.out";
			svar >> fname;
			cout << "Writing Iteration " << iCounter << " weight variance to file " << fname << endl;
			cout << "Writing after " << (float) accCnt << " samples" << endl;
			chkwrite=this->crf_ptr->writeToFile(fname.c_str(),lambdaVar,lambdaLen);
			if (!chkwrite) {
				cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
				exit(-1);
			}
			ftr_str->rewind();
			ftr_str->nextseg();
			iCounter++;
			uCounter=0;
			totLogLi=0.0;
			start=true;
		}
	}
	cout << "Writing Final Iteration weights to file " << this->weight_fname << endl;
	this->crf_ptr->writeToFile(this->weight_fname);
	string fname;
	stringstream savg;
	savg << this->weight_fname << ".avg.out";
	savg >> fname;
	cout << "Writing Final Iteration average weights to file " << fname << endl;
	this->crf_ptr->writeToFile(fname.c_str(),lambdaAvg,lambdaLen);
	stringstream svar;
	svar << this->weight_fname << ".var.out";
	svar >> fname;
	cout << "Writing Final Iteration weight variance to file " << fname << endl;
	delete[] lambdaAvg;
	delete[] lambdaVar;
	delete[] lambdaSqrAcc;
	delete[] grad;
}
