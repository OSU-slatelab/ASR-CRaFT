/*
 * CRF_StdNStateNode.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_StdNStateNode.h"

/*
 * CRF_StdNStateNode constructor
 *
 * see CRF_StateNode for details
 */
CRF_StdNStateNode::CRF_StdNStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf)
	: CRF_StateNode(fb, sizeof_fb, lab, crf)
{
	//QNUInt32 nLabs=this->crf_ptr->getNLabs();
	this->nStates=this->crf_ptr->getFeatureMap()->getNumStates();
	if (nLabs % this->nStates != 0) {
		string errstr="CRF_StdNStateNode constructor: nLabs and nStates do not correspond";
		throw runtime_error(errstr);
	}
	this->nFullLabs = nLabs/this->nStates;
	this->stateArray = new double[nLabs];
	this->denseTransMatrix = new double[nFullLabs*nFullLabs]; // Dense transition matrix
	this->diagTransMatrix = new double[nLabs];
	this->offDiagTransMatrix = new double[nLabs];
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
 * CRF_StdNStateNode destructor
 */
CRF_StdNStateNode::~CRF_StdNStateNode()
{
	delete [] this->diagTransMatrix;
	delete [] this->offDiagTransMatrix;
	delete [] this->denseTransMatrix;
	delete [] this->alphaArray;
	delete [] this->betaArray;
	delete [] this->alphaBetaArray;
	delete [] this->stateArray;
	delete [] this->tempBeta;
	delete [] this->logAddAcc;
}

/*
 * CRF_StdNStateNode::computeTransMatrix
 *
 * Computes the log of the state vector (stateArray).
 * computes the log of the transition matrix as three separate elements:
 *   diagTransMatrix (controlling the self transitions for all states
 *   offDiagTransMatrix (controlling transitions between internal states)
 *   denseTransMatrix (controlling transitions from exit states to entrance states)
 */
double CRF_StdNStateNode::computeTransMatrix()
{
	double result=0.0;
	QNUInt32 nLabs = this->crf_ptr->getNLabs();


	double* lambda = this->crf_ptr->getLambda();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {

		this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(this->ftrBuf,lambda,clab);

		// Here we need to work some magic.  All entries on the diagonal get their self transition assigned
		//this->diagTransMatrix[clab]=this->crf_ptr->getFeatureMap()->computeMij(this->ftrBuf,lambda,lc,clab,clab);
		this->diagTransMatrix[clab]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,clab,clab);

		// Besides the diagonal, we need to update the off-diagonal or the dense transition matrix
		// dense transition matrix is updated when the current label is a start state (e.g. when the
		// clab % nStates == 0
		if (clab % this->nStates == 0) {
			for (QNUInt32 plab=0; plab<nFullLabs; plab++) {
				// Our index into the dense transition matrix needs to be munged a bit
				QNUInt32 idx=plab*nFullLabs+clab/this->nStates;
				// And our previous label is actually the end state for the previous "label"
				QNUInt32 real_plab = plab*nStates+nStates-1;  // real_plab is the end state
				this->denseTransMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,real_plab,clab);
			}
		}
		else {
			// We're on the off-diagonal - clab is not a start  state and the only previous
			// transition must be from the immediate prior state (e.g. clab-1)
			this->offDiagTransMatrix[clab-1]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,clab-1,clab);
		}
	}
	return result;
}

/*
 * CRF_StdNStateNode::computeAlpha
 *
 * Input: prev_alpha - alpha vector of the previous node
 *
 * Computes the alpha vector for the forward backward computation for this node, restricted by
 *  the n-state topology.
 */
double CRF_StdNStateNode::computeAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		// Next we add in transitions from prior states
		// Compute the self transition - all labels get this
		this->alphaArray[clab]=prev_alpha[clab]+this->diagTransMatrix[clab];

		if (clab % this->nStates == 0) {
			// Here clab is a new start state, so all end state transitions to it must be computed
			QNUInt32 dense_clab = clab/this->nStates; //Used to index into dense transition matrix
			this->logAddAcc[0]=prev_alpha[nStates-1]+this->denseTransMatrix[dense_clab];
			double max=this->logAddAcc[0];
			for (QNUInt32 plab=1; plab<nFullLabs; plab++)
			{
				QNUInt32 real_prev = plab*nStates+nStates-1;
				QNUInt32 idx = plab*nFullLabs+dense_clab;
				this->logAddAcc[plab]=prev_alpha[real_prev]+this->denseTransMatrix[idx];
				if (this->logAddAcc[plab]>max) {
					max=this->logAddAcc[plab];
				}
			}
			double logSum=logAdd(this->logAddAcc,max,nFullLabs);
			this->alphaArray[clab]=logAdd(this->alphaArray[clab],logSum);
		}
		else {
			// Here clab is an interior state, so only transitions from the previous state need to be
			// accounted for
			this->alphaArray[clab]=logAdd(this->alphaArray[clab],prev_alpha[clab-1]+this->offDiagTransMatrix[clab-1]);
		}


		// And of course we multiply by the state feature array at the end
	 	this->alphaArray[clab]+=this->stateArray[clab];
	 }
	 this->alphaScale=0.0;

	return this->alphaScale;

}

/*
 * CRF_StdNStateNode::computeFirstAlpha
 *
 * Computes the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdNStateNode::computeFirstAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		// And of course we multiply by the state feature array at the end
		this->alphaArray[clab]=this->stateArray[clab];
	}
	this->alphaScale=0.0;

	return this->alphaScale;

}

/*
 * CRF_StdNStateNode::computeBeta
 *
 * Inputs: result_beta - vector to store the result of the computation
 *         scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Computes the beta vector for the node before this one and store it in result_beta, restricted by
 * n-state topology.
 */
double CRF_StdNStateNode::computeBeta(double* result_beta, double scale)
{
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]

	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->tempBeta[clab]=this->betaArray[clab]+this->stateArray[clab];
	}
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, this->num_labs, 1, this->num_labs, 1.0f, Mtrans, this->num_labs, this->tmp_beta, this->num_labs, 1.0f, new_beta,1);
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLabs, 1, nLabs, 1.0f, this->transMatrix, nLabs, this->tempBeta, nLabs, 1.0f, result_beta,1);

	for (QNUInt32 plab=0; plab<nLabs; plab++) {
		// Here the logic is the inverse of that used in the computeAlphas above
		// because the loop driver here is the previous label

		// First we add the self transition to the result
		result_beta[plab]=this->tempBeta[plab]+this->diagTransMatrix[plab];

		// Check to see if plab is an end state.  If it isn't then we only need to worry about adding in
		// plab+1

		if ( (plab+1) % this->nStates == 0) {
			QNUInt32 idx_plab = (plab+1)/this->nStates -1;
			this->logAddAcc[0]=this->denseTransMatrix[idx_plab*this->nFullLabs]+this->tempBeta[0];
			double max=this->logAddAcc[0];
			for (QNUInt32 clab=1; clab<this->nFullLabs; clab++) {
				QNUInt32 idx=idx_plab*this->nFullLabs+clab;
				QNUInt32 real_clab = clab*nStates;
				logAddAcc[clab]=this->denseTransMatrix[idx]+this->tempBeta[real_clab];
				if (this->logAddAcc[clab]>max) {
					max=this->logAddAcc[clab];
				}
			}
			double logSum=logAdd(this->logAddAcc,max,nFullLabs);
			result_beta[plab]=logAdd(result_beta[plab],logSum);
		}
		else {
			// plab isn't an end state - we only have to compute for the next state (e.g. plab+1)
			if (plab != nLabs-1) {
				// Sanity check, plab should never be the last label in the sequence (that should be
				// an end state).  Just in case
				// Remember - the offDiagTransMatrix is stored by clab index, which is always one
				// more than the corresponding plab index (so in this case, the index into the
				// ofDiagTransMatrix should just be plab, while the "current label" referenced in
				// tempBeta is indexed by plab+1
				result_beta[plab]=logAdd(result_beta[plab],this->offDiagTransMatrix[plab]+this->tempBeta[plab+1]);
			}
			else {
				string errstr="CRF_StdNStateNode::computeBeta() found odd error that needs to be checked";
				throw runtime_error(errstr);
			}
		}
	}

	return this->alphaScale;
}

/*
 * CRF_StdNStateNode::computeAlphaBeta
 *
 * Inputs: Zx - normalization constant
 *
 * Returns: array of state probabilities based on the precomputed alpha and beta for this node
 *
 */
double* CRF_StdNStateNode::computeAlphaBeta(double Zx)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaBetaArray[clab]=this->alphaArray[clab]+this->betaArray[clab]-Zx;
	}
	return this->alphaBetaArray;
}

/*
 * CRF_StdNStateNode::setTailBeta
 *
 * Sets the beta value in this node to the special case for the end of the sequence.
 */
void CRF_StdNStateNode::setTailBeta()
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->betaArray[clab]=0.0;
	}
}

/*
 * CRF_StdNStateNode::computeExpF
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
double CRF_StdNStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	double logLi=0.0;
	double alpha_beta=0.0;
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	//QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
		alpha_beta_tot += alpha_beta;
		logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,this->label,clab);
		// With the new transition matrices, we need to be careful.  The order of this computation needs to
		// match the computeTransitionMatrix code above
		if (prev_lab > nLabs) {
			// if prev_lab > nLabs, we're in the first label frame and there are no previous
			// transitions - skip the transition calculation in this case
			// but set the alpha_beta_trans_tot to 1.0 for the check below
			alpha_beta_trans_tot=1.0;
		}
		else {
			// First we compute the self transitions
			alpha_beta=expE(prev_alpha[clab]+this->diagTransMatrix[clab]+this->stateArray[clab]+this->betaArray[clab]-Zx);
			alpha_beta_trans_tot+=alpha_beta;
			logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,clab,clab);

			// Next we check to see if the clab state is an end state or not
			if (clab % this->nStates == 0) {
				QNUInt32 clab_idx=clab/this->nStates;
				// If so, we update the lambda values in a loop as above
				for (QNUInt32 plab=0; plab<nFullLabs; plab++) {
					// Our index into the dense transition matrix needs to be munged a bit
					QNUInt32 idx=plab*nFullLabs+clab_idx;
					// And our previous label is actually the end state for the previous "label"
					QNUInt32 real_plab = plab*nStates+nStates-1;  // real_plab is the end state
					alpha_beta=expE(prev_alpha[real_plab]+this->denseTransMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
					alpha_beta_trans_tot+=alpha_beta;
					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,real_plab,clab);
				}
			}
			else {
				// If not, we only update the lambda values on the offDiagonal
				alpha_beta=expE(prev_alpha[clab-1]+this->offDiagTransMatrix[clab-1]+this->stateArray[clab]+this->betaArray[clab]-Zx);
				alpha_beta_trans_tot+=alpha_beta;
				logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,clab-1,clab);
			}
		}
	}
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
		string errstr="CRF_StdNStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;

}
/*
 * CRF_StdNStateNode::computeAlphaSum
 *
 * Returns: Sum of the values in the alpha vector of this node
 *
 * Used to compute the normalization constant for the CRF.
 */
double CRF_StdNStateNode::computeAlphaSum()
{
	double Zx = 0.0;
	QNUInt32 nLabs=this->crf_ptr->getNLabs();

	// Added by Ryan, just for debugging
//	for (int i = 0; i < nLabs; i++)
//	{
//		cout << "alphaArray[" << i << "] = " << alphaArray[i] << endl;
//	}

	Zx=logAdd(alphaArray,nLabs);
	return Zx;
}


/*
 * CRF_StdStateNode::getTransValue
 *
 * Returns: Transition matrix value for transition prev_lab->cur_lab
 */
double CRF_StdNStateNode::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	if ((cur_lab % this->nStates == 0) && ((prev_lab+1) % this->nStates == 0)) {
		QNUInt32 cur_idx = cur_lab/this->nStates;
		QNUInt32 prev_idx = (prev_lab+1)/this->nStates -1;
		QNUInt32 idx=prev_idx*this->nFullLabs+cur_idx;
		return this->denseTransMatrix[idx];
	}
	else if (cur_lab == prev_lab) {
		return this->diagTransMatrix[cur_lab];
	}
	else if (cur_lab == prev_lab+1) {
		return this->offDiagTransMatrix[cur_lab-1];
	}
	else {
		// Changed by Ryan
		// In theory, it should return LOG0 instead of 0.
		// But we would report error here instead.
		//return 0;
		string errstr="CRF_StdNStateNode::getTransValue() threw exception: Invalid transition from previous lab (="
				+ stringify(prev_lab) + ") to current lab (=" + stringify(cur_lab) + ").";
		throw runtime_error(errstr);
	}
}

/*
 * CRF_StdNStateNode::getStateValue
 *
 * Returns: State vector value for state label cur_lab
 */
double CRF_StdNStateNode::getStateValue(QNUInt32 cur_lab)
{
	return this->stateArray[cur_lab];
}


/*
 * CRF_StdNStateNode::getFullTransValue
 *
 * Returns: Full value for transition prev_lab->cur_lab
 *   (i.e. transition value prev_lab->cur_lab plus the state value for cur_lab)
 */
double CRF_StdNStateNode::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	if ((cur_lab % this->nStates == 0) && ((prev_lab+1) % this->nStates == 0)) {
		QNUInt32 cur_idx = cur_lab/this->nStates;
		QNUInt32 prev_idx = (prev_lab+1)/this->nStates -1;
		QNUInt32 idx=prev_idx*this->nFullLabs+cur_idx;
		return this->denseTransMatrix[idx]+this->stateArray[cur_lab];
	}
	else if (cur_lab == prev_lab) {
		return this->diagTransMatrix[cur_lab]+this->stateArray[cur_lab];
	}
	else if (cur_lab == prev_lab+1) {
		return this->offDiagTransMatrix[cur_lab-1]+this->stateArray[cur_lab];
	}
	else {
		// Changed by Ryan
		// Instead of return LOG0, we would report error here.
		//return LOG0;
		string errstr="CRF_StdNStateNode::getFullTransValue() threw exception: Invalid transition from previous lab (="
				+ stringify(prev_lab) + ") to current lab (=" + stringify(cur_lab) + ").";
		throw runtime_error(errstr);
	}
}

/*
 * CRF_StdNStateNode::computeAlphaAlignedSum
 *
 * Returns: Sum of the values in the alphaAligned vector of this node
 *
 * Used to compute the "soft" normalization constant.
 */
double CRF_StdNStateNode::computeAlphaAlignedSum()
{
	double Zx;
	try {
		Zx=logAdd(&(this->alphaArrayAligned),nLabs);
	}
	catch (exception& e) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: "+string(e.what());
		string errstr="CRF_StdNStateNode::computeExpF() threw exception: "+string(e.what());
		throw runtime_error(errstr);
	}
	return Zx;
}

