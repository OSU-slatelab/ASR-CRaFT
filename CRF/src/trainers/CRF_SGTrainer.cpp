/*
 * CRF_SGTrainer.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_SGTrainer.h"

/*
 * CRF_SGTrainer constructor
 *
 * See superclass CRF_Trainer for details
 */
CRF_SGTrainer::CRF_SGTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str_mgr, char* wt_fname)
	: CRF_Trainer(crf_in, ftr_str_mgr, wt_fname)
{
	// Added by Ryan
	this->eta = 1.0;
	this->useAdagrad = false;
	this->eps = 1e-12;
}

// Added by Ryan
void CRF_SGTrainer::setNThreads(int n) {
	if (n <= 0) {
		cerr << "CRF_SGTrainer::setNThreads() Error: the number of threads has to be positive." << endl;
		exit(-1);
	}
	this->nThreads = n;
}

// Added by Ryan
void CRF_SGTrainer::setMinibatch(int m) {
	if (m <= 0) {
		cerr << "CRF_SGTrainer::setMinibatch() Error: minibatch size has to be positive." << endl;
		exit(-1);
	}
	this->minibatch = m;
}

// Added by Ryan
void CRF_SGTrainer::setEta(double eta) {
	this->eta = eta;
}

// Added by Ryan
void CRF_SGTrainer::setUseAdagrad(double useAdagrad) {
	this->useAdagrad = useAdagrad;
}

// Added by Ryan
/*
 *
 * The training function that is called outside.
 *
 */
void CRF_SGTrainer::train()
{
	if (this->minibatch == 1) {
		//sgtrain();
		sgtrainMinibatch();
	} else {
		sgtrainMinibatch();
	}
}

/*
 * CRF_SGTrainer::sgtrainMinibatch
 *
 * Performs stochastic gradient training with minibatch
 */
void CRF_SGTrainer::sgtrainMinibatch()
{
	double logLi = 0.0;
	double tmpLogLi;
	double tmp_Zx;
	double totLogLi = 0.0;

	// changed by Ryan
	//int iCounter = 0;
	int iCounter = this->crf_ptr->getInitIter();

	// Added by Ryan
	// rewind the training feature stream iCounter-1 times to make sure the
	// random sequence generator is consistent with the last iteration where
	// the previous training stopped at.
	for (int i = 0; i < iCounter; ++i) {
		this->ftr_strm_mgr->trn_stream->rewind();
		if (this->nThreads > 1) {
			this->ftr_strm_mgr->rewindAllChildrenTrn();
		}
	}

	cout << "Before loading feature streams ..." << endl;

	QNUInt32 uCounter = 0;
	// TODO: get all the children feature streams
	//CRF_FeatureStream *ftr_str=this->ftr_strm_mgr->trn_stream;
	/*
	CRF_FeatureStream **ftr_strs = new CRF_FeatureStream*[this->nThreads];
	if (this->nThreads == 1) {
		ftr_strs[0] = this->ftr_strm_mgr->trn_stream;
		if (!ftr_strs[0]) {
			cerr << "CRF_SGTrainer::sgtrainMinibatch() Error: feature stream is NULL for thread 0" << endl;
			exit(1);
		}
		cout << "After loading feature stream 0 ..." << endl;
	} else {
		for (int i = 0; i < this->nThreads; ++i) {
			CRF_FeatureStreamManager *child_ftr_strm_mgr = this->ftr_strm_mgr->getChild(i);
			if (!child_ftr_strm_mgr) {
				cerr << "CRF_SGTrainer::sgtrainMinibatch() Error: feature stream manager is NULL for thread " << i << endl;
				exit(1);
			}
			ftr_strs[i] = child_ftr_strm_mgr->trn_stream;
			if (!ftr_strs[i]) {
				cerr << "CRF_SGTrainer::sgtrainMinibatch() Error: feature stream is NULL for thread " << i << endl;
				exit(1);
			}
			cout << "After loading feature stream " << i << " ..." << endl;
		}
	}
	*/
	ofstream ofile;

	double* lambda = this->crf_ptr->getLambda();
	QNUInt32 lambdaLen = this->crf_ptr->getLambdaLen();
	//double* lambdaAcc = new double[lambdaLen];
	double* lambdaAcc = this->crf_ptr->getLambdaAcc();

	// Commented by Ryan:
	// This would be screwed if resuming from previous interrupted training,
	// since we have to read from previous variance file to resume this statistics.
	// But since this variable is not used, we are good for now.
	double* lambdaSqrAcc = new double[lambdaLen];

	double* lambdaAvg = new double[lambdaLen];
	double* lambdaVar = new double[lambdaLen];
	//int accCnt = 0;
	int accCnt = this->crf_ptr->getPresentations();
	double* grad = new double[lambdaLen];

	// Added by Ryan, for AdaGrad
	double* gradSqrAcc = this->crf_ptr->getGradSqrAcc();

	for (QNUInt32 i=0; i<lambdaLen; i++) {
		grad[i]=0.0;
		lambdaSqrAcc[i]=0.0;
		//lambdaAcc[i]=0.0;
	}

	// Added by Ryan
	for (QNUInt32 i=0; i<lambdaLen; i++) {
		if (accCnt == 0) {
			lambdaAvg[i]=0.0;
			lambdaVar[i]=0.0;
		} else {
			lambdaAvg[i]=(lambdaAcc[i]/(float)accCnt);
			lambdaVar[i]=(lambdaSqrAcc[i]/(float)accCnt)-(lambdaAvg[i]*lambdaAvg[i]);
		}
	}

	bool start=true;
	float invSquareVar=0.0;
	if (this->useGvar) {
		invSquareVar=1/this->gvar;
	}

	// Modified by Ryan, use the multithreaded version.
	//CRF_GradBuilder* gbuild = CRF_GradBuilder::create(this->crf_ptr,this->objective);
	//gbuild->setNodeList(new CRF_StateVector());
	/*
	CRF_GradAccumulator *gaccum =
			new CRF_Pthread_GradAccumulator(this->crf_ptr, this->useLogspace,
					this->crf_ptr->getFeatureMap()->getNumStates());
	*/
	CRF_Minibatch_GradAccumulator *gaccum = new CRF_Minibatch_GradAccumulator(
			this->crf_ptr, this->ftr_strm_mgr, this->nThreads);
	gaccum->setMinibatch(this->minibatch);
	gaccum->setUttReport(this->uttRpt);
	cout << "Utterance report interval: " << this->uttRpt << endl;
	cout << "Using Logspace training..." << endl;

	cout << "Before rewinding all feature streams ..." << endl;

	// TODO: rewind all feature streams
	//ftr_str->rewind();
	//QN_SegID segid = ftr_str->nextseg();
	/*
	for (int i = 0; i < this->nThreads; ++i) {
		ftr_strs[i]->rewind();
		cout << "After rewinding feature stream " << i << " ..." << endl;
	}
	*/
	gaccum->rewindAllAndNextSegs();

	time_t rawtime;

	// Added by Ryan, just for debugging
//	string testfname;
//	stringstream testss;
//	testss << this->weight_fname << ".i" << iCounter << ".out.beforetrain";
//	testss >> testfname;
//	this->crf_ptr->writeToFile(testfname.c_str());

	while (iCounter < this->maxIters) {

		// Added by Ryan, just for debugging
//		if (uCounter == 0) {
//			testss.clear();
//			testss << this->weight_fname << ".i" << iCounter << ".out.beforetrain";
//			testss >> testfname;
//			this->crf_ptr->writeToFile(testfname.c_str());
//			double ls = 0.0;
//			for (QNUInt32 i=0; i<lambdaLen; i++) {
//				ls += lambda[i];
//				cout << lambda[i] << "\t" << ls << endl;
//			}
//		}

		// Added by Ryan, just for debugging
//		double lambda_sum = 0.0, lambdaAcc_sum = 0.0, lambdaAvg_sum = 0.0;
//		for (QNUInt32 i=0; i<lambdaLen; i++) {
//			lambda_sum += lambda[i];
//			lambdaAcc_sum += lambdaAcc[i];
//			lambdaAvg_sum += lambdaAvg[i];
//		}
//		cout << " Before utt " << uCounter << " : ";
//		cout << " accCnt: " << accCnt;
//		cout << " lambda_sum: " << lambda_sum;
//		cout << " lambdaAcc_sum: " << lambdaAcc_sum;
//		cout << " lambdaAvg_sum: " << lambdaAvg_sum;
//		cout << " lambda[0]: " << lambda[0];
//		cout << " lambdaAcc[0]: " << lambdaAcc[0];
//		cout << " lambdaAvg[0]: " << lambdaAvg[0];
//		cout << " lambda[" << lambdaLen-1 << "]: " << lambda[lambdaLen-1];
//		cout << " lambdaAcc[" << lambdaLen-1 << "]: " << lambdaAcc[lambdaLen-1];
//		cout << " lambdaAvg[" << lambdaLen-1 << "]: " << lambdaAvg[lambdaLen-1];
//		cout << " grad[0]: " << grad[0];
//		cout << endl;

		if (start) {
			// Added the AdaGrad condition by Ryan
			if (useAdagrad) {
				cout << "Iteration: " << iCounter << " starting AdaGrad scaling factor (eta): " << this->eta << endl;
			} else {
				cout << "Iteration: " << iCounter << " starting LR: " << this->lr << endl;
			}
			start=false;
		}
		if (uCounter % this->uttRpt == 0) {
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time << " Beginning Utt: " << uCounter <<  endl;
			//delete time;
		}

		QNUInt32 uCounterInc = 0;
		bool isEndOfIter = false;
		tmp_Zx = uCounter%this->uttRpt; // Added for Ferr testing
		// TODO: use the multi-threaded version
		tmpLogLi = gaccum->accumulateGradient(grad, &tmp_Zx, &uCounterInc, &isEndOfIter);
		//tmpLogLi=gbuild->buildGradient(ftr_str,grad,&tmp_Zx);

		//Added by Ryan
		//Parameter tying.
		//this->crf_ptr->getFeatureMap()->tieGradient(grad, 11, 488);
//		this->crf_ptr->getFeatureMap()->tieGradient(grad);

		// TODO: accumulate uCounter, make the reporting block unconditional
		//tmpLogLi=logLi;
		logLi = tmpLogLi - tmp_Zx;
		totLogLi += logLi;
		uCounter += uCounterInc;
		// report information after each minibatch
		if (uCounter % this->uttRpt == 0) {
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time;
			cout << " Finished Utt: " << uCounter - 1;
			cout << " Batch-Avg Numerator: " << tmpLogLi / uCounterInc;
			cout << " Batch-Avg Zx: " << tmp_Zx / uCounterInc;
			cout << " Batch-Avg LogLi: " << logLi / uCounterInc;
			cout << " Iter-Avg LogLi: " << totLogLi / uCounter;

			// Added by Ryan, just for debugging
//			cout << " accCnt: " << accCnt;
//			cout << " lambda[0]: " << lambda[0];
//			cout << " lambdaAcc[0]: " << lambdaAcc[0];
//			cout << " lambdaAvg[0]: " << lambdaAvg[0];
//			cout << " grad[0]: " << grad[0];

			cout << endl;
		}
		//uCounter++;
		for (QNUInt32 i=0; i<lambdaLen; i++) {
			if (this->useGvar) {
				//grad[i]-=lambda[i]*invSquareVar;
				grad[i]-=grad[i]*invSquareVar;
			}

			// just for debugging
//			cout << "before update, lambda[" << i << "]=" << lambda[i] << endl;

			// Added the AdaGrad condition by Ryan
			if (useAdagrad) {
				//if (std::abs(gradSqrAcc[i]-0.0) < 1e-18) {
				//	cerr << "Error: zero accumulated gradient square sum: gradSqrAcc[" << i << "]=" << gradSqrAcc[i] << endl;
				//	exit(-1);
				//}
				gradSqrAcc[i] += grad[i]*grad[i];
				lambda[i] += this->eta / (std::sqrt(gradSqrAcc[i]) + this->eps) * grad[i];
			} else {
				lambda[i] += this->lr * grad[i];
			}
			lambdaAcc[i] += lambda[i];
			lambdaSqrAcc[i] += lambda[i] * lambda[i];
			grad[i] = 0.0;

			// just for debugging
//			cout << "after update, lambda[" << i << "]=" << lambda[i] << endl;
		}
		//accCnt++;
		accCnt += uCounterInc;

		//cout << "Utterance " << uCounter << ": Lambda[0]: " << lambda[0] << endl;
		//cout << "Utterance " << uCounter << ": Acc Lambda[0]: " << (lambdaAcc[0]/accCnt) << endl;
		//segid = ftr_str->nextseg();

		// TODO: check if EOS
		//if (segid == QN_SEGID_BAD) {
		if (isEndOfIter) {
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

			// commented out by Ryan: not outputting the variance to file for saving space
//			stringstream svar;
//			svar << this->weight_fname << ".i" << iCounter << ".var.out";
//			svar >> fname;
//			cout << "Writing Iteration " << iCounter << " weight variance to file " << fname << endl;
//			cout << "Writing after " << (float) accCnt << " samples" << endl;
//			chkwrite=this->crf_ptr->writeToFile(fname.c_str(),lambdaVar,lambdaLen);
//			if (!chkwrite) {
//				cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
//				exit(-1);
//			}

			if (useAdagrad) {
				stringstream sgradsqr;
				sgradsqr << this->weight_fname << ".i" << iCounter << ".gradSqrAcc.out";
				sgradsqr >> fname;
				cout << "Writing Iteration " << iCounter << " accumulated gradient square sum to file " << fname << endl;
				cout << "Writing after " << (float) accCnt << " samples" << endl;
				chkwrite=this->crf_ptr->writeToFile(fname.c_str(),gradSqrAcc,lambdaLen);
				if (!chkwrite) {
					cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
					exit(-1);
				}
			}

			// TODO: rewind all the thread feature streams
			//ftr_str->rewind();
			//ftr_str->nextseg();
			/*
			for (int i = 0; i < this->nThreads; ++i) {
				ftr_strs[i]->rewind();
			}
			*/
			gaccum->rewindAllAndNextSegs();

			touchDoneFileIter(iCounter);

			iCounter++;
			uCounter=0;
			totLogLi=0.0;
			start=true;

			// Added by Ryan, learning rate decay
			if (!useAdagrad) {
				this->lr *= this->lr_decay_rate;
				cout << "Learning rate is decayed by " << this->lr_decay_rate << " to " << this->lr << endl;
			}
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

	touchDoneFileFinal();

	// Commented out by Ryan: not outputting the variance to file for saving space
//	stringstream svar;
//	svar << this->weight_fname << ".var.out";
//	svar >> fname;
//	cout << "Writing Final Iteration weight variance to file " << fname << endl;
//	this->crf_ptr->writeToFile(fname.c_str(),lambdaVar,lambdaLen);

	delete[] lambdaAvg;
	delete[] lambdaVar;
	delete[] lambdaSqrAcc;
	delete[] grad;

	// Added by Ryan
	//delete[] ftr_strs;
	delete gaccum;
}

/*
 * CRF_SGTrainer::sgtrain
 *
 * Performs stochastic gradient training.
 */
void CRF_SGTrainer::sgtrain()
{
	double logLi = 0.0;
	double tmpLogLi;
	double tmp_Zx;
	double totLogLi = 0.0;

	// changed by Ryan
	//int iCounter = 0;
	int iCounter = this->crf_ptr->getInitIter();

	// Added by Ryan
	// rewind the training feature stream iCounter-1 times to make sure the
	// random sequence generator is consistent with the last iteration where
	// the previous training stopped at.
	for (int i = 0; i < iCounter; ++i) {
		this->ftr_strm_mgr->trn_stream->rewind();
		if (this->nThreads > 1) {
			this->ftr_strm_mgr->rewindAllChildrenTrn();
		}
	}

	int uCounter=0;
	CRF_FeatureStream *ftr_str = this->ftr_strm_mgr->trn_stream;
	ofstream ofile;

	CRF_GradBuilder* gbuild = CRF_GradBuilder::create(this->crf_ptr,this->objective);
	gbuild->setNodeList(new CRF_StateVector());
	cout << "Using Logspace training..." << endl;

	double* lambda = this->crf_ptr->getLambda();
	QNUInt32 lambdaLen = this->crf_ptr->getLambdaLen();
	//double* lambdaAcc = new double[lambdaLen];
	double* lambdaAcc = this->crf_ptr->getLambdaAcc();

	// Commented by Ryan:
	// This would be screwed if resuming from previous interrupted training,
	// since we have to read from previous variance file to resume this statistics.
	// But since this variable is not used, we are good for now.
	double* lambdaSqrAcc=new double[lambdaLen];

	double* lambdaAvg = new double[lambdaLen];
	double* lambdaVar = new double[lambdaLen];
	//int accCnt = 0;
	int accCnt = this->crf_ptr->getPresentations();
	double* grad = new double[lambdaLen];

	// Added by Ryan, for AdaGrad
	double* gradSqrAcc = this->crf_ptr->getGradSqrAcc();

	for (QNUInt32 i=0; i<lambdaLen; i++) {
		grad[i]=0.0;
		lambdaSqrAcc[i]=0.0;
		//lambdaAcc[i]=0.0;
	}

	// Added by Ryan
	for (QNUInt32 i=0; i<lambdaLen; i++) {
		if (accCnt == 0) {
			lambdaAvg[i]=0.0;
			lambdaVar[i]=0.0;
		} else {
			lambdaAvg[i]=(lambdaAcc[i]/(float)accCnt);
			lambdaVar[i]=(lambdaSqrAcc[i]/(float)accCnt)-(lambdaAvg[i]*lambdaAvg[i]);
		}
	}

	bool start=true;
	float invSquareVar=0.0;
	if (this->useGvar) {
		invSquareVar=1/this->gvar;
	}
	ftr_str->rewind();
	QN_SegID segid = ftr_str->nextseg();

    time_t rawtime;

	// Added by Ryan, just for debugging
//	string testfname;
//	stringstream testss;
//	testss << this->weight_fname << ".i" << iCounter << ".out.beforetrain";
//	testss >> testfname;
//	this->crf_ptr->writeToFile(testfname.c_str());

	while (iCounter<this->maxIters) {

		// Added by Ryan, just for debugging
//		if (uCounter == 0) {
//			testss.clear();
//			testss << this->weight_fname << ".i" << iCounter << ".out.beforetrain";
//			testss >> testfname;
//			this->crf_ptr->writeToFile(testfname.c_str());
//			double ls = 0.0;
//			for (QNUInt32 i=0; i<lambdaLen; i++) {
//				ls += lambda[i];
//				cout << lambda[i] << "\t" << ls << endl;
//			}
//		}

		// Added by Ryan, just for debugging
//		double lambda_sum = 0.0, lambdaAcc_sum = 0.0, lambdaAvg_sum = 0.0;
//		for (QNUInt32 i=0; i<lambdaLen; i++) {
//			lambda_sum += lambda[i];
//			lambdaAcc_sum += lambdaAcc[i];
//			lambdaAvg_sum += lambdaAvg[i];
//		}
//		cout << " Before utt " << uCounter << " : ";
//		cout << " accCnt: " << accCnt;
//		cout << " lambda_sum: " << lambda_sum;
//		cout << " lambdaAcc_sum: " << lambdaAcc_sum;
//		cout << " lambdaAvg_sum: " << lambdaAvg_sum;
//		cout << " lambda[0]: " << lambda[0];
//		cout << " lambdaAcc[0]: " << lambdaAcc[0];
//		cout << " lambdaAvg[0]: " << lambdaAvg[0];
//		cout << " lambda[" << lambdaLen-1 << "]: " << lambda[lambdaLen-1];
//		cout << " lambdaAcc[" << lambdaLen-1 << "]: " << lambdaAcc[lambdaLen-1];
//		cout << " lambdaAvg[" << lambdaLen-1 << "]: " << lambdaAvg[lambdaLen-1];
//		cout << " grad[0]: " << grad[0];
//		cout << endl;

		if (start) {
			// Added the AdaGrad condition by Ryan
			if (useAdagrad) {
				cout << "Iteration: " << iCounter << " starting AdaGrad scaling factor (eta): " << this->eta << endl;
			} else {
				cout << "Iteration: " << iCounter << " starting LR: " << this->lr << endl;
			}
			start=false;
		}
		if (uCounter % this->uttRpt == 0) {
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time << " Beginning Utt: " << uCounter << " SegID: " << segid <<  endl;
			//delete time;
		}

		tmp_Zx=uCounter%this->uttRpt; // Added for Ferr testing
		tmpLogLi=gbuild->buildGradient(ftr_str,grad,&tmp_Zx);

		//Added by Ryan
		//Parameter tying.
		//this->crf_ptr->getFeatureMap()->tieGradient(grad, 11, 488);
//		this->crf_ptr->getFeatureMap()->tieGradient(grad);

		//tmpLogLi=logLi;
		logLi=tmpLogLi - tmp_Zx;
		totLogLi += logLi;
		if (uCounter % this->uttRpt == 0) {
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time;
			cout << " Finished Utt: " << uCounter;
			cout << " Utt Numerator: " << tmpLogLi;
			cout << " Utt Zx: " << tmp_Zx;
			cout << " Utt LogLi: " << logLi;
			cout << " Iter-Avg LogLi: " << totLogLi/(uCounter+1);

			// Added by Ryan, just for debugging
//			cout << " accCnt: " << accCnt;
//			cout << " lambda[0]: " << lambda[0];
//			cout << " lambdaAcc[0]: " << lambdaAcc[0];
//			cout << " lambdaAvg[0]: " << lambdaAvg[0];
//			cout << " grad[0]: " << grad[0];

			cout << endl;
		}
		uCounter++;
		for (QNUInt32 i=0; i<lambdaLen; i++) {
			if (this->useGvar) {
				//grad[i]-=lambda[i]*invSquareVar;
				grad[i]-=grad[i]*invSquareVar;
			}

			// just for debugging
//			cout << "before update, lambda[" << i << "]=" << lambda[i] << endl;

			// Added the AdaGrad condition by Ryan
			if (useAdagrad) {
				//if (std::abs(gradSqrAcc[i]-0.0) < 1e-18) {
				//	cerr << "Error: zero accumulated gradient square sum: gradSqrAcc[" << i << "]=" << gradSqrAcc[i] << endl;
				//	exit(-1);
				//}
				gradSqrAcc[i] += grad[i]*grad[i];
				lambda[i] += this->eta / (std::sqrt(gradSqrAcc[i]) + this->eps) * grad[i];
			} else {
				lambda[i] += this->lr * grad[i];
			}
			lambdaAcc[i] += lambda[i];
			lambdaSqrAcc[i] += lambda[i] * lambda[i];
			grad[i] = 0.0;

			// just for debugging
//			cout << "after update, lambda[" << i << "]=" << lambda[i] << endl;
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

			// commented out by Ryan: not outputting the variance to file for saving space
//			stringstream svar;
//			svar << this->weight_fname << ".i" << iCounter << ".var.out";
//			svar >> fname;
//			cout << "Writing Iteration " << iCounter << " weight variance to file " << fname << endl;
//			cout << "Writing after " << (float) accCnt << " samples" << endl;
//			chkwrite=this->crf_ptr->writeToFile(fname.c_str(),lambdaVar,lambdaLen);
//			if (!chkwrite) {
//				cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
//				exit(-1);
//			}

			if (useAdagrad) {
				stringstream sgradsqr;
				sgradsqr << this->weight_fname << ".i" << iCounter << ".gradSqrAcc.out";
				sgradsqr >> fname;
				cout << "Writing Iteration " << iCounter << " accumulated gradient square sum to file " << fname << endl;
				cout << "Writing after " << (float) accCnt << " samples" << endl;
				chkwrite=this->crf_ptr->writeToFile(fname.c_str(),gradSqrAcc,lambdaLen);
				if (!chkwrite) {
					cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
					exit(-1);
				}
			}

			ftr_str->rewind();
			ftr_str->nextseg();

			touchDoneFileIter(iCounter);

			iCounter++;
			uCounter=0;
			totLogLi=0.0;
			start=true;

			// Added by Ryan, learning rate decay
			if (!useAdagrad) {
				this->lr *= this->lr_decay_rate;
				cout << "Learning rate is decayed by " << this->lr_decay_rate << " to " << this->lr << endl;
			}
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

	touchDoneFileFinal();

	// Commented out by Ryan: not outputting the variance to file for saving space
//	stringstream svar;
//	svar << this->weight_fname << ".var.out";
//	svar >> fname;
//	cout << "Writing Final Iteration weight variance to file " << fname << endl;
//	this->crf_ptr->writeToFile(fname.c_str(),lambdaVar,lambdaLen);

	delete[] lambdaAvg;
	delete[] lambdaVar;
	delete[] lambdaSqrAcc;
	delete[] grad;
}
