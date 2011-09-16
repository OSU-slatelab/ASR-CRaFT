/*
 * CRF_StdSegStateNode.cpp
 *
 *  Created on: Sep 13, 2011
 *      Author: hey
 */

#include "CRF_StdSegStateNode.h"

CRF_StdSegStateNode::CRF_StdSegStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur)
	: CRF_StdStateNode(fb, sizeof_fb, lab, crf)
{
	this->labMaxDur = crf_ptr->getLabMaxDur();
	this->nodeLabMaxDur = nodeMaxDur;
	if (nodeLabMaxDur > labMaxDur)
	{
		string errstr="CRF_StdSegStateNode constructor: nodeLabMaxDur is larger than labMaxDur.";
		throw runtime_error(errstr);
	}
	if (nLabs % this->labMaxDur != 0) {
		string errstr="CRF_StdSegStateNode constructor: nLabs and labMaxDur do not correspond.";
		throw runtime_error(errstr);
	}
	this->nFullLabs = nLabs/this->labMaxDur;
	for (QNUInt32 clab = 0; clab < nLabs; clab++)
	{
		this->alphaArray[clab] = CRF_LogMath::LOG0;
		this->betaArray[clab] = CRF_LogMath::LOG0;
		this->alphaBetaArray[clab] = CRF_LogMath::LOG0;
	}

	// TODO: verify sizeof_fb since the feature buffer now contains the features for more than one segments.
	if (ftrBuf_size % this->nodeLabMaxDur != 0) {
		string errstr="CRF_StdSegStateNode constructor: ftrBuf_size and nodeLabMaxDur do not correspond.";
		throw runtime_error(errstr);
	}
	this->nFtrsPerSeg = ftrBuf_size / this->nodeLabMaxDur;
}

CRF_StdSegStateNode::~CRF_StdSegStateNode() {
	// TODO Auto-generated destructor stub
}

/*
 * CRF_StdSegStateNode::computeTransMatrix
 *
 * Computes the log of the state vector (stateArray) and the log of the transition matrix (transMatrix)
 *   and stores them as appropriate
 */
double CRF_StdSegStateNode::computeTransMatrix()
{
	double result=0.0;

	double* lambda = this->crf_ptr->getLambda();
	float* seg_ftr_buf = this->ftrBuf;
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		if (clab % this->labMaxDur == 0)
		{
			seg_ftr_buf = this->ftrBuf;
		}
		if (clab % this->labMaxDur >= this->nodeLabMaxDur)
		{
			continue;
		}
		this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,clab);
		for (QNUInt32 plab=0; plab<nLabs; plab++) {
			QNUInt32 idx=plab*nLabs+clab;
			this->transMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(seg_ftr_buf,lambda,plab,clab);
		}
		seg_ftr_buf += this->nFtrsPerSeg;
	}
	return result;
}

/*
 * CRF_StdSegStateNode::computeAlpha
 *
 * Input: prev_alpha - vector of alpha vectors of a few previous nodes
 *
 * Computes the alpha vector for the forward backward computation for this node.
 */
double CRF_StdSegStateNode::computeAlpha(double** prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		if (clab % this->labMaxDur >= this->nodeLabMaxDur)
		{
			continue;
		}
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
	 		string errstr="CRF_StdSegStateNode::computeAlpha() caught exception: "+string(e.what())+" while computing alpha";
	 		throw runtime_error(errstr);
			return(-1);
	 	}
	 	this->alphaArray[clab]+=this->stateArray[clab];
	 }

	return this->alphaScale;

}

/*
 * CRF_StdSegStateNode::computeFirstAlpha
 *
 * This version is to match the input type of computeAlpha(double** prev_alpha).
 *
 * Computes the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdSegStateNode::computeFirstAlpha(double** prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab+=this->labMaxDur) {
	 	this->alphaArray[clab]=this->stateArray[clab];
	 }
	return this->alphaScale;

}

/*
 * CRF_StdSegStateNode::computeBeta
 *
 * Inputs: result_beta - vector of vectors to store the result of beta vectors of a few previous nodes.
 *         scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Computes the beta vector for the node before this one and store it in result_beta
 */
double CRF_StdSegStateNode::computeBeta(double** result_beta, double scale)
{
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]

	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		if (clab % this->labMaxDur >= this->nodeLabMaxDur)
		{
			continue;
		}
		this->tempBeta[clab]=this->betaArray[clab]+this->stateArray[clab];
	}

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
			string errstr="CRF_StdSegStateNode::computeBeta() caught exception: "+string(e.what())+" while computing beta";
			throw runtime_error(errstr);
			return(-1);
		}
	}
	return this->alphaScale;
}

/*
 * CRF_StdSegStateNode::computeAlphaBeta
 *
 * Inputs: Zx - normalization constant
 *
 * Returns: array of state probabilities based on the precomputed alpha and beta for this node
 *
 */
double* CRF_StdSegStateNode::computeAlphaBeta(double Zx)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		if (clab % this->labMaxDur >= this->nodeLabMaxDur)
		{
			continue;
		}
		this->alphaBetaArray[clab]=this->alphaArray[clab]+this->betaArray[clab]-Zx;
	}
	return this->alphaBetaArray;
}

/*
 * CRF_StdSegStateNode::setTailBeta
 *
 * Sets the beta value in this node to the special case for the end of the sequence.
 */
void CRF_StdSegStateNode::setTailBeta()
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->betaArray[clab]=0.0;
	}
}

/*
 * CRF_StdSegStateNode::computeExpF
 *
 * Inputs: *ExpF - vector to store expected values of feature functions
 *         *grad - vector to store computed gradient values
 *         Zx - normalization constant
 *         **prev_alpha - vector of vectors of alpha values from a few previous nodes (for use in transition feature
 *            ExpF computation
 *         prev_lab - previous node label (transition feature ExpF computation)
 *
 * Returns:
 *
 * Computes gradient and expected values for features in this node and store them in *grad and
 *   *ExpF vectors respectively.  State features and transition features are computed in the same function.
 */
double CRF_StdSegStateNode::computeExpF(double* ExpF, double* grad, double Zx, double** prev_alpha, QNUInt32 prev_lab)
{
	double logLi=0.0;
	double alpha_beta=0.0;
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
		alpha_beta_tot += alpha_beta;
		bool match=(clab==this->label);
		logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,this->label,clab);
		if (prev_lab > nLabs) {
			// if prev_lab > nLabs, we're in the first label frame and there are no previous
			// transitions - skip the transition calculation in this case
			// but set the alpha_beta_trans_tot to 1.0 for the check below
			alpha_beta_trans_tot=1.0;
		}
		else {
			// Otherwise do the transition calculations
			for (QNUInt32 plab=0; plab<nLabs; plab++) {
				QNUInt32 idx = plab*nLabs+clab;
				alpha_beta=expE(prev_alpha[plab]+this->transMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
				alpha_beta_trans_tot+=alpha_beta;
				match=((clab==this->label)&&(plab==prev_lab));
				logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab);
			}
		}
	}

	//Added by Ryan, just for debugging
	//cout << "\tAlpha_beta_tot: " << alpha_beta_tot << "\tAlpha_beta_trans_tot: " << alpha_beta_trans_tot;

	if ((alpha_beta_tot >1.1))  {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot < 0.9) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.1) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot < 0.9) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;
}


// These functions do not work for CRF_StdSegStateNode

/*
 * CRF_StdSegStateNode::computeAlpha
 *
 * Input: prev_alpha - alpha vector of the previous node
 *
 * Computes the alpha vector for the forward backward computation for this node.
 */
double CRF_StdSegStateNode::computeAlpha(double* prev_alpha)
{
	string errstr="CRF_StdSegStateNode::computeAlpha(double* prev_alpha) threw exception: CRF_StdSegStateNode::computeAlpha(double** prev_alpha) should be used instead.";
	throw runtime_error(errstr);
}

/*
 * CRF_StdSegStateNode::computeFirstAlpha
 *
 * Computes the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdSegStateNode::computeFirstAlpha(double* prev_alpha)
{
	string errstr="CRF_StdSegStateNode::computeFirstAlpha(double* prev_alpha) threw exception: CRF_StdSegStateNode::computeFirstAlpha(double** prev_alpha) should be used instead.";
	throw runtime_error(errstr);
}

/*
 * CRF_StdSegStateNode::computeBeta
 *
 * Inputs: result_beta - vector to store the result of the computation
 *         scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Computes the beta vector for the node before this one and store it in result_beta
 */
double CRF_StdSegStateNode::computeBeta(double* result_beta, double scale)
{
	string errstr="CRF_StdSegStateNode::computeBeta(double* result_beta, double scale) threw exception: CRF_StdSegStateNode::computeBeta(double** result_beta, double scale) should be used instead.";
	throw runtime_error(errstr);
}

/*
 * CRF_StdSegStateNode::computeExpF
 *
 * Inputs: *ExpF - vector to store expected values of feature functions
 *         *grad - vector to store computed gradient values
 *         Zx - normalization constant
 *         *prev_alpha - vector of alpha values from the previous node (for use in transition feature
 *            ExpF computation
 *         prev_lab - previous node label (transition feature ExpF computation)
 *
 * Returns:
 *
 * Computes gradient and expected values for features in this node and store them in *grad and
 *   *ExpF vectors respectively.  State features and transition features are computed in the same function.
 */
double CRF_StdSegStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	string errstr="CRF_StdSegStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab) threw exception: CRF_StdSegStateNode::computeExpF(double* ExpF, double* grad, double Zx, double** prev_alpha, QNUInt32 prev_lab) should be used instead.";
	throw runtime_error(errstr);
}
