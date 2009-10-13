/*
 * CRF_StdStateNodeMasked.cpp
 *
 *  Created on: Jul 29, 2009
 *      Author: morrijer
 */

#include "CRF_StdStateNodeMasked.h"

CRF_StdStateNodeMasked::CRF_StdStateNodeMasked(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf)
	: CRF_StdStateNode(fb,sizeof_fb, lab,crf)
{
	// TODO Auto-generated constructor stub

}

CRF_StdStateNodeMasked::~CRF_StdStateNodeMasked() {
	// TODO Auto-generated destructor stub
}

double CRF_StdStateNodeMasked::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
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
		tmp_score=this->alphaArray[clab]*this->betaArray[clab]*this->alphaScale/Zx;
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
			tmp_score=prev_alpha[plab]*this->transMatrix[idx]*this->stateArray[best_clab]*this->betaArray[best_clab]/Zx;
			if (tmp_score>best_plab_score) {
				best_plab_score=tmp_score;
				best_plab=plab;
			}
		}
	}
	bool compute_clab_grad=!(this->label == best_clab);
	bool compute_plab_grad=!(this->label == best_clab && prev_lab==best_plab);
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		alpha_beta=this->alphaArray[clab]*this->betaArray[clab]*this->alphaScale/Zx;
		alpha_beta_tot += alpha_beta;
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
				alpha_beta=prev_alpha[plab]*this->transMatrix[idx]*this->stateArray[clab]*this->betaArray[clab]/Zx;
				alpha_beta_trans_tot+=alpha_beta;
				logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab,compute_plab_grad);
			}
		}
	}

	if ((alpha_beta_tot >1.1))  {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot < 0.9) {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.1) {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot < 0.9) {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;

}
