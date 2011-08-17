/*
 * CRF_StdNStateNodeMasked.cpp
 *
 *  Created on: Jul 29, 2009
 *      Author: morrijer
 */

#include "CRF_StdNStateNodeMasked.h"

CRF_StdNStateNodeMasked::CRF_StdNStateNodeMasked(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf)
	:CRF_StdNStateNode(fb,sizeof_fb,lab,crf)

{
	// TODO Auto-generated constructor stub

}

CRF_StdNStateNodeMasked::~CRF_StdNStateNodeMasked() {
	// TODO Auto-generated destructor stub
}

double CRF_StdNStateNodeMasked::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	double logLi=0.0;
	double alpha_beta=0.0;
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;

	QNUInt32 best_clab=0; double best_clab_score=this->alphaArray[0]*this->betaArray[0]*this->alphaScale/Zx;
	double tmp_score;
	for (QNUInt32 clab=1; clab<nLabs; clab++) {
		// We need to find the best label score
		tmp_score=this->alphaArray[clab]*this->betaArray[clab]*this->alphaScale/Zx;
		if (tmp_score>best_clab_score) {
			best_clab_score=tmp_score;
			best_clab=clab;
		}
	}
	// Start with the assumption that the best transition is a self transition
	QNUInt32 best_plab=best_clab;
	double best_plab_score=prev_alpha[best_clab]*this->diagTransMatrix[best_clab]*this->stateArray[best_clab]*this->betaArray[best_clab]/Zx;
	if (prev_lab<=nLabs) {
		if (best_clab % this->nStates == 0) {
		    for (QNUInt32 plab=0; plab<nLabs; plab++) {
			    QNUInt32 clab_idx=best_clab/this->nStates;
			    QNUInt32 idx=plab*nFullLabs+clab_idx;
			    QNUInt32 real_plab = plab*nStates+nStates-1;  // real_plab is the end state
			    tmp_score=prev_alpha[real_plab]*this->denseTransMatrix[idx]*this->stateArray[best_clab]*this->betaArray[best_clab]/Zx;
			    if (tmp_score>best_plab_score) {
			    	best_plab_score=tmp_score;
			    	best_plab=plab;
			    }
		    }
		}
		else {
			tmp_score=prev_alpha[best_clab-1]*this->offDiagTransMatrix[best_clab-1]*this->stateArray[best_clab]*this->betaArray[best_clab]/Zx;
			if (tmp_score>best_plab_score) {
				best_plab_score=tmp_score;
				best_plab=best_clab-1;
			}
		}
	}
	bool compute_clab_grad=!(this->label == best_clab);
	bool compute_plab_grad=!(this->label == best_clab && prev_lab==best_plab);


	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		alpha_beta=this->alphaArray[clab]*this->betaArray[clab]*this->alphaScale/Zx;
		/*cout << "CLAB: " << clab << " Alpha: " << this->alphaArray[clab] << " Beta: " << this->betaArray[clab];
		cout << " Scale: " << this->alphaScale << " Zx: " << Zx << endl;*/
		alpha_beta_tot += alpha_beta;
		//bool match=(clab==this->label);
		logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,this->label,clab,compute_clab_grad);
		//logLi+=this->crf_ptr->getFeatureMap()->computeExpFState(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab);
		// With the new transition matrices, we need to be careful.  The order of this computation needs to
		// match the computeTransitionMatrix code above

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
			// First we compute the self transitions
			alpha_beta=prev_alpha[clab]*this->diagTransMatrix[clab]*this->stateArray[clab]*this->betaArray[clab]/Zx;
			alpha_beta_trans_tot+=alpha_beta;
			//match=((clab==this->label)&&(clab==prev_lab));
			logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,clab,clab,compute_plab_grad);
			//logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab,clab);

			// Next we check to see if the clab state is an end state or not
			if (clab % this->nStates == 0) {
				// If so, we update the lambda values in a loop as above
				for (QNUInt32 plab=0; plab<nFullLabs; plab++) {
					// Our index into the dense transition matrix needs to be munged a bit
					QNUInt32 idx=plab*nFullLabs+clab/this->nStates;
					// And our previous label is actually the end state for the previous "label"
					QNUInt32 real_plab = plab*nStates+nStates-1;  // real_plab is the end state
					alpha_beta=prev_alpha[real_plab]*this->denseTransMatrix[idx]*this->stateArray[clab]*this->betaArray[clab]/Zx;
					alpha_beta_trans_tot+=alpha_beta;
					//match=((clab==this->label)&&(real_plab==prev_lab));
					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,real_plab,clab,compute_plab_grad);
					//logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,real_plab,clab);
				}
			}
			else {
				// If not, we only update the lambda values on the offDiagonal
				alpha_beta=prev_alpha[clab-1]*this->offDiagTransMatrix[clab-1]*this->stateArray[clab]*this->betaArray[clab]/Zx;
				alpha_beta_trans_tot+=alpha_beta;
				//match=((clab==this->label)&&(clab-1==prev_lab));
				logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,clab-1,clab,compute_plab_grad);
				//logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab-1,clab);
			}
		}
	}
	/*cout << "Alpha_beta tot: " << alpha_beta_tot << endl;
	cout << "Alpha_beta_trans tot: " << alpha_beta_trans_tot << endl;*/
	if ((alpha_beta_tot >1.1))  {
		string errstr="CRF_StdNStateNode::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot < 0.9) {
		string errstr="CRF_StdNStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.1) {
		string errstr="CRF_StdNStateNode::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot < 0.9) {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;

}

