#include "CRF_LBFGSTrainer.h"

CRF_LBFGSTrainer::CRF_LBFGSTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str_mgr, char* wt_fname)
	: CRF_Trainer(crf_in, ftr_str_mgr, wt_fname)
{
}

void CRF_LBFGSTrainer::train()
{
	ofstream ofile;
	QNUInt32 nStates = this->crf_ptr->getFeatureMap()->getNumStates();

	if (this->ftr_strm_mgr->getNThreads()==1)
		gaccum = new CRF_GradAccumulator(this->crf_ptr,this->useLogspace,nStates);
	else
		gaccum = new CRF_Pthread_GradAccumulator(this->crf_ptr,this->useLogspace,nStates);
	this->gaccum->setUttReport(this->uttRpt);

	//	gbuild = CRF_GradBuilder::create(this->crf_ptr,this->useLogspace,nStates);


	double* lambda=this->crf_ptr->getLambda();
	QNUInt32 lambdaLen = this->crf_ptr->getLambdaLen();
	this->grad=new double[lambdaLen];

	// lbfgs likes to have their own malloc (sigh)
	lbfgsfloatval_t *lambdaCoef= lbfgs_malloc(lambdaLen);
	lbfgsfloatval_t logLikelihood = 0.0;

	for (QNUInt32 i=0; i<lambdaLen; i++) {
		lambdaCoef[i]=lambda[i];
		grad[i]=0.0;
	}

	this->start=true;
	this->invSquareVar=0.0;
	if (this->useGvar) {
		invSquareVar=1/(this->gvar*this->gvar);
	}
	lbfgs_parameter_t lbfgs_params;
	lbfgs_parameter_init(&lbfgs_params);
	lbfgs_params.epsilon=0.001;
	lbfgs_params.m=7;
	lbfgs_params.xtol=1.0e-16;

	//lbfgs_params.linesearch=LBFGS_LINESEARCH_BACKTRACKING;

	// Start a lbfgs loop.  Callback to gradient evaluator
	int ret=lbfgs(lambdaLen,lambdaCoef,&logLikelihood,_evaluateGradient,_progress,this,NULL);

	if (ret==0) {
		cout << "Writing Final Iteration weights to file " << this->weight_fname << endl;
		this->crf_ptr->writeToFile(this->weight_fname);
	} else {
		cerr << "LBFGS returned error: " << ret << endl;
		//exit(-1);
	}

	lbfgs_free(lambdaCoef);
	//delete gbuild;
	delete gaccum;
	delete[] grad;
}


// Compute the gradient and return the log likelihood
lbfgsfloatval_t CRF_LBFGSTrainer::evaluateGradient(const lbfgsfloatval_t *lambda,
        lbfgsfloatval_t *gradient,
        const int lambdaLen,
        const lbfgsfloatval_t step
        ) {

	//QN_SegID segid = QN_SEGID_BAD;
	//CRF_FeatureStream *ftr_str=this->ftr_strm_mgr->trn_stream;

	/*
	if (this->gbuild) {
		// Rewind the stream

		ftr_str->rewind();
		segid = ftr_str->nextseg();

		// SHOULD CODE BETTER TO THROW EXCEPTION
		if (segid == QN_SEGID_BAD) {
			cerr << "Feature stream contains no utterances!" << endl;
			exit(1);
		}
	}
	*/

    time_t rawtime;
    double totLogLi=0.0;
    unsigned int uCounter=0;
    //double tmp_Zx;

    //double *crflambda=this->crf_ptr->getLambda();

    for (QNUInt32 i=0; i<lambdaLen; i++) {
    	(this->crf_ptr->getLambda())[i]=(double)(lambda[i]);
    	this->grad[i]=0.0;
    }


	// Only write previous iteration's weights if this is not first iteration
	if (this->start) {
		this->start=false;
	} else {
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
	}

	iCounter++;
	cout << "Iteration: " << this->iCounter  << endl;

	/*
	if (this->gbuild) {
		// now loop over all utterances
		do {
			double tmpLogLi=this->gbuild->buildGradient(ftr_str,this->grad,&tmp_Zx);
			double logLi=tmpLogLi - tmp_Zx;
			totLogLi += logLi;

			if (uCounter % this->uttRpt == 0) {
				time(&rawtime);
				char* time = ctime(&rawtime);
				time[strlen(time)-1]='\0';
				cout << time << " Finished Utt: " << uCounter << " logLi: " << logLi;
				cout << " Avg LogLi: " << totLogLi/(uCounter+1) << " Zx: " << tmp_Zx;
				cout << " Numerator: " << tmpLogLi << endl;
			}
			uCounter++;


			segid = ftr_str->nextseg();
		} while (segid != QN_SEGID_BAD);
		for (QNUInt32 i=0; i<lambdaLen; i++) {
			// note: lbfgs wants to minimize, so we invert sign to maximize
			gradient[i]=-(lbfgsfloatval_t)(this->grad[i]);
		}

		for (QNUInt32 i=0; i<lambdaLen; i++) {
			if (this->useGvar) {
				// add penalty term (because of inversion)
				gradient[i]+=lambda[i]*invSquareVar;
				totLogLi-=((lambda[i]*lambda[i])*invSquareVar)/2;
			}
			cout << this->crf_ptr->getFeatureMap()->getMapDescriptor(i) << ":" << gradient[i] << ":" << this->grad[i] << ":" << i << endl;
		}

	} else {
	*/
	totLogLi=this->gaccum->accumulateGradient(this->ftr_strm_mgr,this->ftr_strm_mgr->getNThreads(),this->grad,&uCounter);
	//totLogLi=this->gaccum->accumulateGradient(this->ftr_strm_mgr,1,this->grad,&uCounter);

	// invert gradient for use in lbfgs
	if (this->useGvar) {
		for (QNUInt32 i=0; i<lambdaLen; i++) {
			// add penalty term (because of inversion)
			gradient[i]=(lbfgsfloatval_t)(-this->grad[i]+lambda[i]*invSquareVar);
			totLogLi-=((lambda[i]*lambda[i])*invSquareVar)/2;
		}
		//cout << this->crf_ptr->getFeatureMap()->getMapDescriptor(i) << ":" << gradient[i] << ":" << this->grad[i] << ":" << i << endl;
	} else {
		for (QNUInt32 i=0; i<lambdaLen; i++) {
			gradient[i]=(lbfgsfloatval_t)(-this->grad[i]);
		}
	}
/*}*/

    time(&rawtime);
    char* time = ctime(&rawtime);
    time[strlen(time)-1]='\0';
    cout << time << " End iteration: " << iCounter << " totLogLi: " << totLogLi;
    cout << " Avg LogLi: " << totLogLi/(uCounter) << " ucounter: " << uCounter << endl;

    /*if (this->gbuild) {
    	// clean up stream
    	ftr_str->rewind();
    	ftr_str->nextseg();
    }*/

	// return the negative log likelihood, which wants to be minimized by lbfgs
    return -totLogLi;
}

 lbfgsfloatval_t CRF_LBFGSTrainer::_evaluateGradient(
       void *instance,
       const lbfgsfloatval_t *x,
       lbfgsfloatval_t *g,
       const int n,
       const lbfgsfloatval_t step
       )
   {
       return reinterpret_cast<CRF_LBFGSTrainer *>(instance)->evaluateGradient(x, g, n, step);
   }

// note: progress is not used because lbfgs does not call progress
// on every iteration while in a linesearch
int CRF_LBFGSTrainer::progress(
        const lbfgsfloatval_t *lambda,
        const lbfgsfloatval_t *gradient,
        const lbfgsfloatval_t logLikelihood,
        const lbfgsfloatval_t lambdaNorm,
        const lbfgsfloatval_t gradientNorm,
        const lbfgsfloatval_t step,
        int lambdaLen,
        int lbfgsIter,
        int numberOfEvals
        )
    {

	cout << "PROGRESS called" << endl;
	cout << "Iteration: " << this->iCounter << " ending" << endl;
	cout << "Iteration: " << this->iCounter << " applying lambda updates" << endl;


	return (iCounter>=this->maxIters)?1:0;
}

int CRF_LBFGSTrainer::_progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
	    return reinterpret_cast<CRF_LBFGSTrainer *>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }


