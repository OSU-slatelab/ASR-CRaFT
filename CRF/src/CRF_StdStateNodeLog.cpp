#include "CRF_StdStateNodeLog.h"

CRF_StdStateNodeLog::CRF_StdStateNodeLog(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf)
	: CRF_StdStateNode(fb, sizeof_fb, lab, crf)
{
	QNUInt32 nLabs=this->crf_ptr->getNLabs();
	this->logAddAcc = new double[nLabs];
}

CRF_StdStateNodeLog::~CRF_StdStateNodeLog()
{
	delete this->logAddAcc;
}

double CRF_StdStateNodeLog::computeTransMatrix()
{
	return this->computeTransMatrixLog();
}

double CRF_StdStateNodeLog::computeAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=0.0; // Reset alphaArray to zero before processing
	}

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->logAddAcc[0]=prev_alpha[0]+this->transMatrix[0+clab];
		double maxv=this->logAddAcc[0];
		for (QNUInt32 plab=1; plab<nLabs; plab++) {
	 		this->logAddAcc[plab]=prev_alpha[plab]+this->transMatrix[plab*nLabs+clab];
	 		if (this->logAddAcc[plab]>maxv) {
	 			maxv=logAddAcc[plab];
	 		}
	 	}
	 	try {
	 		this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,nLabs);
	 	}
	 	catch (exception &e) {
	 		string errstr="CRF_StdStateNodeLog::computeAlpha caught exception: "+string(e.what())+" while computing alpha";
			throw runtime_error(errstr);
			return(-1);
	 	}
	 	this->alphaArray[clab]+=this->stateArray[clab];
	 	//alpha_tot+=new_alpha[clab];
	 }


	return this->alphaScale;

}

double CRF_StdStateNodeLog::computeFirstAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;
	//for (QNUInt32 clab=0; clab<nLabs; clab++) {
	//	this->alphaArray[clab]=0.0; // Reset alphaArray to zero before processing
	//}

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
	 	this->alphaArray[clab]=this->stateArray[clab];
	 	//alpha_tot+=new_alpha[clab];
	 }


	return this->alphaScale;

}

double CRF_StdStateNodeLog::computeBeta(double* result_beta, double scale)
{
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]

	QNUInt32 nLabs = this->crf_ptr->getNLabs();
		//for (QNUInt32 clab=0; clab<nLabs; clab++) {
		//	this->betaArray[clab]=0.0;
		//}
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		//this->betaArray[clab]=0.0;
		this->tempBeta[clab]=this->betaArray[clab]+this->stateArray[clab];
		result_beta[clab]=0;
	}
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, this->num_labs, 1, this->num_labs, 1.0f, Mtrans, this->num_labs, this->tmp_beta, this->num_labs, 1.0f, new_beta,1);
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLabs, 1, nLabs, 1.0f, this->transMatrix, nLabs, this->tempBeta, nLabs, 1.0f, result_beta,1);

	for (QNUInt32 plab=0; plab<nLabs; plab++) {
		this->logAddAcc[0]=this->transMatrix[plab*nLabs+0]+this->tempBeta[0];
		double maxv=this->logAddAcc[0];
		for (QNUInt32 clab=1; clab<nLabs; clab++) {
			this->logAddAcc[clab]=this->transMatrix[plab*nLabs+clab]+this->tempBeta[clab];
			if (this->logAddAcc[clab]>maxv) {
				maxv=this->logAddAcc[clab];
			}
		}
		try {
			result_beta[plab]=logAdd(this->logAddAcc,maxv,nLabs);
		}
		catch (exception &e) {
	 		string errstr="CRF_StdStateNodeLog::computeBeta caught exception: "+string(e.what())+" while computing beta";
			throw runtime_error(errstr);
			return(-1);
		}
	}

	return this->alphaScale;
}

double* CRF_StdStateNodeLog::computeAlphaBeta(double Zx)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaBetaArray[clab]=this->alphaArray[clab]+this->betaArray[clab]-Zx;
	}
	return this->alphaBetaArray;
}


void CRF_StdStateNodeLog::setTailBeta()
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->betaArray[clab]=0.0;
	}
}

double CRF_StdStateNodeLog::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	// Reconsider this logic a bit - perhaps this should be moved back to the gradient builder and
	// a simpler set of "computeAlphaBeta" and "computeTransAlphaBeta" functions should be implemented
	// in the statenodes instead...

	// (Although if we do the above that makes it harder to plug in arbitrary state-paths, so perhaps
	//  we should just keep this the way it is).
	double logLi=0.0;
	double alpha_beta=0.0;
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
		alpha_beta_tot += alpha_beta;
		bool match=(clab==this->label);
		logLi+=this->crf_ptr->getFeatureMap()->computeExpFState(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab);
		if (prev_lab > nLabs) {
			// if prev_lab > nLabs, we're in the first label frame and there are no previous
			// transitions - skip the transition calculation in this case
			// but set the alpha_beta_trans_tot to 1.0 for the check below
			alpha_beta_trans_tot=1.0;
			// We STILL need to cycle through our "lc" values because of how the FeatureMap code works
			// Our feature map is expecting us to have lc point at the next input vector
			//double* tmp_ExpF = new double[this->crf_ptr->getLambdaLen()];
			for (QNUInt32 plab=0; plab<nLabs; plab++) {
				lc+=this->crf_ptr->getFeatureMap()->getNumTransFuncs(plab,clab);
				//this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,NULL,lc,tmp_ExpF,NULL,alpha_beta,false,plab,clab);
			}
			//delete tmp_ExpF;
		}
		else {
			// if prev_lab > nLabs, we're in the first label frame and there are no previous
			// transitions - skip the transition calculation in this case
			for (QNUInt32 plab=0; plab<nLabs; plab++) {
				QNUInt32 idx = plab*nLabs+clab;
				//alpha_beta=prev_alpha[plab]*this->transMatrix[idx]*this->stateArray[clab]*this->betaArray[clab]/Zx;
				alpha_beta=expE(prev_alpha[plab]+this->transMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
				alpha_beta_trans_tot+=alpha_beta;
				match=((clab==this->label)&&(plab==prev_lab));
				//logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab,plab);
				logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,plab,clab);
			}
		}
	}

	if ((alpha_beta_tot >1.1))  {
		string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot < 0.9) {
		string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.1) {
		string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot < 0.9) {
		string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;
}

double CRF_StdStateNodeLog::computeAlphaSum()
{
	//double Zx = 0.0;
	//QNUInt32 nLabs=this->crf_ptr->getNLabs();
	double Zx;
	try {
		Zx=logAdd(this->alphaArray,this->crf_ptr->getNLabs());
	}
	catch (exception& e) {
		string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: "+string(e.what());
		throw runtime_error(errstr);
	}
	return Zx;
}

double CRF_StdStateNodeLog::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	return this->transMatrix[prev_lab*nLabs+cur_lab]+this->stateArray[cur_lab];
}
