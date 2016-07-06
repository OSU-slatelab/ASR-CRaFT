/*
 * CRF_StdStateNode.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_StdStateNode.h"

/*
 * CRF_StdStateNode constructor
 *
 * Input: fb - feature buffer for features applicable for this node
 *        sizeof_fb - number of features in buffer fb
 *        lab - label for this node
 *        crf_in - pointer back to CRF model used by this node
 */
CRF_StdStateNode::CRF_StdStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf)
	: CRF_StateNode(fb, sizeof_fb, lab, crf)
{

	this->stateArray = new double[nLabs];
	this->transMatrix = new double[nLabs*nLabs];
	this->alphaArray = new double[nLabs];
	this->betaArray = new double[nLabs];
	this->alphaBetaArray = new double[nLabs];
	this->tempBeta = new double[nLabs];
	this->alphaSize = nLabs;
	this->alphaScale = 0.0;
	this->logAddAcc = new double[nLabs];
	this->alphaArrayAligned.assign(nLabs,CRF_LogMath::LOG0);
	this->betaArrayAligned.assign(nLabs,CRF_LogMath::LOG0);

}

/*
 * CRF_StdStateNode destructor
 */
CRF_StdStateNode::~CRF_StdStateNode()
{
	delete [] this->stateArray;
	delete [] this->transMatrix;
	delete [] this->alphaArray;
	delete [] this->betaArray;
	delete [] this->tempBeta;
	delete [] this->logAddAcc;

	// Ryan: why not delete alphaBetaArray?
}


/*
 * CRF_StateNode::computeTransMatrix
 *
 * Computes the log of the state vector (stateArray) and the log of the transition matrix (transMatrix)
 *   and stores them as appropriate
 */
double CRF_StdStateNode::computeTransMatrix()
{
	double result=0.0;

	double* lambda = this->crf_ptr->getLambda();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(this->ftrBuf,lambda,clab);

		for (QNUInt32 plab=0; plab<nLabs; plab++) {
			QNUInt32 idx=plab*nLabs+clab;
			this->transMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,plab,clab);
		}
	}
	return result;
}

/*
 * CRF_StdStateNode::computeAlpha
 *
 * Input: prev_alpha - alpha vector of the previous node
 *
 * Computes the alpha vector for the forward backward computation for this node.
 */
double CRF_StdStateNode::computeAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

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
	 		//changed by Ryan
	 		//string errstr="CRF_StdStateNodeLog::computeAlpha caught exception: "+string(e.what())+" while computing alpha";
	 		string errstr="CRF_StdStateNode::computeAlpha() caught exception: "+string(e.what())+" while computing alpha";
	 		throw runtime_error(errstr);
			return(-1);
	 	}
	 	this->alphaArray[clab]+=this->stateArray[clab];
	 }

	return this->alphaScale;

}

/*
 * CRF_StdStateNode::computeFirstAlpha
 *
 * Computes the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdStateNode::computeFirstAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
	 	this->alphaArray[clab]=this->stateArray[clab];
	}
	return this->alphaScale;

}

/*
 * CRF_StateNode::computeBeta
 *
 * Inputs: result_beta - vector to store the result of the computation
 *         scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Computes the beta vector for the node before this one and store it in result_beta
 */
double CRF_StdStateNode::computeBeta(double* result_beta, double scale)
{
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]

	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
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
			//changed by Ryan
	 		//string errstr="CRF_StdStateNodeLog::computeBeta caught exception: "+string(e.what())+" while computing beta";
	 		string errstr="CRF_StdStateNode::computeBeta() caught exception: "+string(e.what())+" while computing beta";
			throw runtime_error(errstr);
			return(-1);
		}
	}
	return this->alphaScale;
}

/*
 * CRF_StdStateNode::computeAlphaBeta
 *
 * Inputs: Zx - normalization constant
 *
 * Returns: array of state probabilities based on the precomputed alpha and beta for this node
 *
 */
double* CRF_StdStateNode::computeAlphaBeta(double Zx)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaBetaArray[clab]=this->alphaArray[clab]+this->betaArray[clab]-Zx;
	}
	return this->alphaBetaArray;
}

/*
 * CRF_StdStateNode::setTailBeta
 *
 * Sets the beta value in this node to the special case for the end of the sequence.
 */
void CRF_StdStateNode::setTailBeta()
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->betaArray[clab]=0.0;
	}
}

/*
 * CRF_StdStateNode::computeExpF
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
double CRF_StdStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
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
	if ((alpha_beta_tot >1.1))  {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot < 0.9) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.1) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot < 0.9) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;
}

/*
 * CRF_StdStateNode::computeAlphaSum
 *
 * Returns: Sum of the values in the alpha vector of this node
 *
 * Used to compute the normalization constant for the CRF.
 */
double CRF_StdStateNode::computeAlphaSum()
{
	double Zx;
	try {
		Zx=logAdd(this->alphaArray,nLabs);
	}
	catch (exception& e) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: "+string(e.what());
		string errstr="CRF_StdStateNode::computeAlphaSum() threw exception: "+string(e.what());
		throw runtime_error(errstr);
	}
	return Zx;
}

/*
 * CRF_StdStateNode::getTransValue
 *
 * Returns: Transition matrix value for transition prev_lab->cur_lab
 */
double CRF_StdStateNode::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	return this->transMatrix[prev_lab*nLabs+cur_lab];
}

/*
 * CRF_StdStateNode::getStateValue
 *
 * Returns: State vector value for state label cur_lab
 */
double CRF_StdStateNode::getStateValue(QNUInt32 cur_lab)
{
	return this->stateArray[cur_lab];
}

/*
 * CRF_StdStateNode::getFullTransValue
 *
 * Returns: Full value for transition prev_lab->cur_lab
 *   (i.e. transition value prev_lab->cur_lab plus the state value for cur_lab)
 */
double CRF_StdStateNode::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	return this->transMatrix[prev_lab*nLabs+cur_lab]+this->stateArray[cur_lab];


}

/*
 * CRF_StdStateNode::computeAlphaAlignedSum
 *
 * Returns: Sum of the values in the alphaAligned vector of this node
 *
 * Used to compute the "soft" normalization constant.
 */
double CRF_StdStateNode::computeAlphaAlignedSum()
{
	double Zx;
	try {
		Zx=logAdd(&(this->alphaArrayAligned),nLabs);
	}
	catch (exception& e) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: "+string(e.what());
		string errstr="CRF_StdStateNode::computeAlphaAlignedSum() threw exception: "+string(e.what());
		throw runtime_error(errstr);
	}
	return Zx;
}


