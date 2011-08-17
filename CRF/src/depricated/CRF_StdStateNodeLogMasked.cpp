/*
 * CRF_StdStateNodeLogMasked.cpp
 *
 *  Created on: Jul 27, 2009
 *      Author: morrijer
 */

#include "CRF_StdStateNodeLogMasked.h"

CRF_StdStateNodeLogMasked::CRF_StdStateNodeLogMasked(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf)
	:CRF_StdStateNodeLog(fb,sizeof_fb,lab,crf)
{
}

CRF_StdStateNodeLogMasked::~CRF_StdStateNodeLogMasked() {
	// TODO Auto-generated destructor stub
}

double CRF_StdStateNodeLogMasked::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	// Reconsider this logic a bit - perhaps this should be moved back to the gradient builder and
	// a simpler set of "computeAlphaBeta" and "computeTransAlphaBeta" functions should be implemented
	// in the statenodes instead...

	// (Although if we do the above that makes it harder to plug in arbitrary state-paths, so perhaps
	//  we should just keep this the way it is).
	double logLi=0.0;
	double alpha_beta=0.0;
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	//QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;
	QNUInt32 best_clab=0; double best_clab_score=this->alphaArray[0]+this->betaArray[0]-Zx;
	double tmp_score;
	for (QNUInt32 clab=1; clab<nLabs; clab++) {
		// We need to find the best label score
		tmp_score=this->alphaArray[clab]+this->betaArray[clab]-Zx;
		if (tmp_score>best_clab_score) {
			best_clab_score=tmp_score;
			best_clab=clab;
		}
	}
	QNUInt32 best_plab=0;
	double best_plab_score=prev_alpha[0]+this->transMatrix[0+best_clab]+this->stateArray[best_clab]+this->betaArray[best_clab]-Zx;
	if (prev_lab<=nLabs) {
		for (QNUInt32 plab=0; plab<nLabs; plab++) {
			// Next we find the best transition into the best clab
			// We only do this if we're not on the first frame of data
			QNUInt32 idx = plab*nLabs+best_clab;
			tmp_score=prev_alpha[plab]+this->transMatrix[idx]+this->stateArray[best_clab]+this->betaArray[best_clab]-Zx;
			if (tmp_score>best_plab_score) {
				best_plab_score=tmp_score;
				best_plab=plab;
			}
		}
	}
//	bool compute_clab_grad=!(this->label == best_clab);
	bool compute_plab_grad=!(this->label == best_clab && prev_lab==best_plab);
	bool compute_clab_grad=compute_plab_grad;  //Use the same criterion for state updates as transition updates
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
		alpha_beta_tot += alpha_beta;
		//bool match=(clab==this->label);
		//logLi+=this->crf_ptr->getFeatureMap()->computeExpFState(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab);
		logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,this->label,clab,compute_clab_grad);
		if (prev_lab > nLabs) {
			// if prev_lab > nLabs, we're in the first label frame and there are no previous
			// transitions - skip the transition calculation in this case
			// but set the alpha_beta_trans_tot to 1.0 for the check below
			alpha_beta_trans_tot=1.0;
			// We STILL need to cycle through our "lc" values because of how the FeatureMap code works
			// Our feature map is expecting us to have lc point at the next input vector
			//double* tmp_ExpF = new double[this->crf_ptr->getLambdaLen()];
			/*for (QNUInt32 plab=0; plab<nLabs; plab++) {
				lc+=this->crf_ptr->getFeatureMap()->getNumTransFuncs(plab,clab);
				//this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,NULL,lc,tmp_ExpF,NULL,alpha_beta,false,plab,clab);
			}*/
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
				compute_plab_grad=!(clab==best_clab && plab==best_plab);
				//match=((clab==this->label)&&(plab==prev_lab));
				//logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,plab,clab);
				logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab,compute_plab_grad);
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
