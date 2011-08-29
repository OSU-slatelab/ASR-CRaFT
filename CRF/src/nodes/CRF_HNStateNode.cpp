#include "CRF_HNStateNode.h"

CRF_HNStateNode::CRF_HNStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf)
	: CRF_StateNode(fb, sizeof_fb, lab, crf)
{
	//QNUInt32 nLabs=this->crf_ptr->getNLabs();
	this->nStates=this->crf_ptr->getFeatureMap()->getNumStates();
	this->trueNLabs = nLabs*this->nStates;
	this->stateArray = new double[this->trueNLabs];
	this->denseTransMatrix = new double[nLabs*nLabs]; // Dense transition matrix
	this->diagTransMatrix = new double[this->trueNLabs];
	this->offDiagTransMatrix = new double[this->trueNLabs];
	this->alphaArray = new double[this->trueNLabs];
	this->betaArray = new double[this->trueNLabs];
	this->alphaBetaArray = new double[this->trueNLabs];
	this->tempBeta = new double[this->trueNLabs];
	this->alphaSize = this->trueNLabs;
	this->alphaScale = 0.0;
	this->logAddAcc = new double[this->trueNLabs];
	this->alphaArrayAligned.assign(this->trueNLabs,CRF_LogMath::LOG0);
	this->betaArrayAligned.assign(this->trueNLabs,CRF_LogMath::LOG0);
}

CRF_HNStateNode::~CRF_HNStateNode()
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


double CRF_HNStateNode::computeTransMatrix()
{
	double result=0.0;
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	double* lambda = this->crf_ptr->getLambda();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {

		for (QNUInt32 cstate=0; cstate<this->nStates; cstate++) {
			QNUInt32 tlab=clab*this->nStates+cstate;
			this->stateArray[tlab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(this->ftrBuf,lambda,clab);
			// All entries on the diagonal get their self transition assigned
			this->diagTransMatrix[tlab]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,tlab,tlab);
			if (cstate==0) {
				// Besides the diagonal, we need to update the off-diagonal or the dense transition matrix
				// dense transition matrix is updated when the current label is a start state
				// (e.g. cstate==0)
				for (QNUInt32 plab=0; plab<nLabs; plab++) {
					// Our index into the dense transition matrix needs to be munged a bit
					QNUInt32 idx=plab*nLabs+clab;
					// And our previous label is actually the end state for the previous "label"
					QNUInt32 real_plab = plab*nStates+nStates-1;  // real_plab is the end state
					this->denseTransMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,real_plab,tlab);
				}
			}
			else {
				// We're on the off-diagonal - clab is not a start  state and the only previous
				// transition must be from the immediate prior state (e.g. clab-1)
				this->offDiagTransMatrix[tlab-1]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,tlab-1,tlab);
			}
		}
	}
}

double CRF_HNStateNode::computeAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		for (QNUInt32 cstate=0; cstate<this->nStates; cstate++) {
			QNUInt32 tlab=clab*this->nStates+cstate;

			// Next we add in transitions from prior states
			// Compute the self transition - all labels get this
			this->alphaArray[tlab]=prev_alpha[tlab]+this->diagTransMatrix[tlab];
			if (cstate == 0) {
				this->logAddAcc[0]=prev_alpha[nStates-1]+this->denseTransMatrix[clab];
				double max=this->logAddAcc[0];
				for (QNUInt32 plab=1; plab<nLabs; plab++)
				{
					QNUInt32 real_prev = plab*nStates+nStates-1;
					QNUInt32 idx = plab*nLabs+clab;
					this->logAddAcc[plab]=prev_alpha[real_prev]+this->denseTransMatrix[idx];
					if (this->logAddAcc[plab]>max) {
						max=this->logAddAcc[plab];
					}
				}
				double logSum=logAdd(this->logAddAcc,max,nLabs);
				this->alphaArray[tlab]=logAdd(this->alphaArray[tlab],logSum);
			}
			else {
				// Here clab is an interior state, so only transitions from the previous state need to be
				// accounted for
				this->alphaArray[tlab]=logAdd(this->alphaArray[tlab],prev_alpha[tlab-1]+this->offDiagTransMatrix[tlab-1]);
			}

			// And of course we multiply by the state feature array at the end
			this->alphaArray[tlab]+=this->stateArray[tlab];
		}
	 }
	 this->alphaScale=0.0;

	return this->alphaScale;
}

double CRF_HNStateNode::computeAlphaAligned(vector<double>* prev_alpha, QNUInt32 prev_label, QNUInt32 next_label)
{
	/*
	 * The following code has been commented out to allow for
	 * compilation.  The code needs some revision to allow it
	 * to perform as intended - jjm 8/29/11
	if (prev_label != this->label) {
		// We're in the first state - only do processing for the initial state
		QNUInt32 tlab=this->label*nStates;
		// And we can only transition from the final state of prev_label
		QNUInt32 t_plab=prev_label*nStates+nStates-1;
		this->alphaArrayAligned[tlab]=prev_alpha[t_plab]+this->denseTransMatrix[prev_label*nLabs+tlab];
	}
	else if (next_label != this->label) {
		// We're in the final state - only do processing for the final state
		QNUInt32 tlab=this->label*nStates+nStates-1;
		// And we can make two transitions, first the transition from the previous state
		this->alphaArrayAligned[tlab]=prev_alpha[tlab-1]+this->offDiagTransMatrix[tlab-1];
		// And the self transition
		this->alphaArrayAligned[tlab]=logAdd(this->alphaArrayAligned[tlab],prev_alpha[tlab]+this->diagTransMatrix[tlab]);
	}
	else {
		// We're in some internal state, so we need to compute for all of the states
		for (QNUInt32 cstate=0; cstate<nStates; cstate++) {
			QNUInt32 tlab=clab*this->nStates+cstate;

			// Next we add in transitions from prior states
			// Compute the self transition - all labels get this
			this->alphaArrayAligned[tlab]=prev_alpha[tlab]+this->diagTransMatrix[tlab];
			if (cstate !=0) {
				// state 0 only gets the self transition
				// successor states also get predecessor state transitions
				this->alphaArrayAligned[tlab]=logAdd(this->alphaArrayAligned[tlab],prev_alpha[tlab-1]+this->offDiagTransMatrix[tlab-1]);
			}

		}

	}
	return this->alphaScale;
	*/
	return 0.0;
}

double CRF_HNStateNode::computeFirstAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	for (QNUInt32 clab=0; clab<nLabs*this->nStates; clab++) {
		// And of course we multiply by the state feature array at the end
		this->alphaArray[clab]=this->stateArray[clab];
	}
	this->alphaScale=0.0;

	return this->alphaScale;

}

double CRF_HNStateNode::computeBeta(double* result_beta, double scale)
{
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]

	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs*this->nStates; clab++) {
		this->tempBeta[clab]=this->betaArray[clab]+this->stateArray[clab];
	}
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, this->num_labs, 1, this->num_labs, 1.0f, Mtrans, this->num_labs, this->tmp_beta, this->num_labs, 1.0f, new_beta,1);
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLabs, 1, nLabs, 1.0f, this->transMatrix, nLabs, this->tempBeta, nLabs, 1.0f, result_beta,1);

	for (QNUInt32 plab=0; plab<nLabs; plab++) {
		// Here the logic is the inverse of that used in the computeAlphas above
		// because the loop driver here is the previous label
		for (QNUInt32 pstate=0; pstate<this->nStates; pstate++) {
			QNUInt32 tplab=plab+pstate;
			// First we add the self transition to the result
			result_beta[tplab]=this->tempBeta[tplab]+this->diagTransMatrix[tplab];
			// Check to see if plab is an end state.  If it isn't then we
			// only need to worry about adding in plab+1
			if ( pstate == this->nStates-1) {
				QNUInt32 idx_plab = (plab+1)/this->nStates -1;
				this->logAddAcc[0]=this->denseTransMatrix[plab*nLabs]+this->tempBeta[0];
				double max=this->logAddAcc[0];
				for (QNUInt32 clab=1; clab<this->nLabs; clab++) {
					QNUInt32 idx=plab*this->nLabs+clab;
					QNUInt32 real_clab = clab*nStates;
					logAddAcc[clab]=this->denseTransMatrix[idx]+this->tempBeta[real_clab];
					if (this->logAddAcc[clab]>max) {
						max=this->logAddAcc[clab];
					}
				}
				double logSum=logAdd(this->logAddAcc,max,nLabs);
				result_beta[tplab]=logAdd(result_beta[tplab],logSum);
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
					result_beta[tplab]=logAdd(result_beta[tplab],this->offDiagTransMatrix[tplab]+this->tempBeta[tplab+1]);
				}
				else {
					string errstr="CRF_HNStateNode::computeBeta() found odd error that needs to be checked";
					throw runtime_error(errstr);
				}
			}
		}
	}

	return this->alphaScale;
}

double* CRF_HNStateNode::computeAlphaBeta(double Zx)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs*nStates; clab++) {
		this->alphaBetaArray[clab]=this->alphaArray[clab]+this->betaArray[clab]-Zx;
	}
	return this->alphaBetaArray;
}


void CRF_HNStateNode::setTailBeta()
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs*nStates; clab++) {
		this->betaArray[clab]=0.0;
	}
}

double CRF_HNStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	double logLi=0.0;
	double alpha_beta=0.0;
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	//QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		for (QNUInt32 cstate=0; cstate<nStates; cstate++) {
			QNUInt32 tlab=clab*nStates+cstate;
			alpha_beta=expE(this->alphaArray[tlab]+this->betaArray[tlab]-Zx);
			alpha_beta_tot+=alpha_beta;
			// Below needs to be fixed to account for multi-state
			//logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,this->label,tlab);
			QNUInt32 t_label=this->label*nStates+cstate;
			logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,t_label,tlab);
			if (prev_lab > nLabs) {
				alpha_beta_trans_tot=1.0;
			}
			else {
				QNUInt32 t_plabel;

				// First we compute the self transitions
				alpha_beta=expE(prev_alpha[tlab]+this->diagTransMatrix[tlab]+this->stateArray[tlab]+this->betaArray[tlab]-Zx);
				alpha_beta_trans_tot+=alpha_beta;
				// Below needs to be fixed to account for multi-state
				logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,t_label,tlab,tlab);

				// Next we check to see if we're in a start state or not
				if (cstate == 0) {
					QNUInt32 clab_idx=clab/this->nStates;
					// If so, we update the lambda values in a loop as above
					for (QNUInt32 plab=0; plab<nLabs; plab++) {
						// Our index into the dense transition matrix needs to be munged a bit
						QNUInt32 idx=plab*nLabs+clab;
						// And our previous label is actually the end state for the previous "label"
						QNUInt32 real_plab = plab*nStates+nStates-1;  // real_plab is the end state
						alpha_beta=expE(prev_alpha[real_plab]+this->denseTransMatrix[idx]+this->stateArray[tlab]+this->betaArray[tlab]-Zx);
						alpha_beta_trans_tot+=alpha_beta;
						// Below needs to be fixed to account for multi-state
						logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,real_plab,clab);
					}
				}
				else {
					// If not, we only update the lambda values on the offDiagonal
					alpha_beta=expE(prev_alpha[tlab-1]+this->offDiagTransMatrix[tlab-1]+this->stateArray[tlab]+this->betaArray[tlab]-Zx);
					alpha_beta_trans_tot+=alpha_beta;
					// Below needs to be fixed to account for multi-state
					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,clab-1,clab);
				}
			}
		}
	}

	if ((alpha_beta_tot >1.1))  {
		string errstr="CRF_HNStateNode::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot < 0.9) {
		string errstr="CRF_HNStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.1) {
		string errstr="CRF_HNStateNode::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot < 0.9) {
		string errstr="CRF_HNStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;

}

double CRF_HNStateNode::computeAlphaSum()
{
	double Zx = 0.0;
	QNUInt32 nLabs=this->crf_ptr->getNLabs();
	Zx=logAdd(alphaArray,nLabs*nStates);
	return Zx;
}

double CRF_HNStateNode::computeAlphaAlignedSum()
{
	double Zx;
	try {
		Zx=logAdd(&(this->alphaArrayAligned),nLabs);
	}
	catch (exception& e) {
		string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: "+string(e.what());
		throw runtime_error(errstr);
	}
	return Zx;
}

double CRF_HNStateNode::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	if ((cur_lab % this->nStates == 0) && ((prev_lab+1) % this->nStates == 0)) {
		QNUInt32 cur_idx = cur_lab/this->nStates;
		QNUInt32 prev_idx = (prev_lab+1)/this->nStates -1;
		QNUInt32 idx=prev_idx*this->nLabs+cur_idx;
		return this->denseTransMatrix[idx];
	}
	else if (cur_lab == prev_lab) {
		return this->diagTransMatrix[cur_lab];
	}
	else if (cur_lab == prev_lab+1) {
		return this->offDiagTransMatrix[cur_lab-1];
	}
	else {
		return 0;
	}
}

double CRF_HNStateNode::getStateValue(QNUInt32 cur_lab)
{
	return this->stateArray[cur_lab];
}



double CRF_HNStateNode::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	if ((cur_lab % this->nStates == 0) && ((prev_lab+1) % this->nStates == 0)) {
		QNUInt32 cur_idx = cur_lab/this->nStates;
		QNUInt32 prev_idx = (prev_lab+1)/this->nStates -1;
		QNUInt32 idx=prev_idx*this->nLabs+cur_idx;
		return this->denseTransMatrix[idx]+this->stateArray[cur_lab];
	}
	else if (cur_lab == prev_lab) {
		return this->diagTransMatrix[cur_lab]+this->stateArray[cur_lab];
	}
	else if (cur_lab == prev_lab+1) {
		return this->offDiagTransMatrix[cur_lab-1]+this->stateArray[cur_lab];
	}
	else {
		return LOG0;
	}
}


/*double CRF_HNStateNode::computeTransMatrix()
{
	double result = this->computeTransMatrixLog();
	QNUInt32 nLabs=this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		try {
			this->stateArray[clab]=expE(this->stateArray[clab]);
		}
		catch (exception &e) {
			string errstr="CRF_HNStateNode::computeTransMatrix() caught exception "+string(e.what())+" while taking exp of stateArray cell"+stringify(clab);
			throw runtime_error(errstr);
		}
		try {
			this->diagTransMatrix[clab]=expE(this->diagTransMatrix[clab]);
		}
		catch (exception &e) {
			string errstr="CRF_HNStateNode::computeTransMatrix() caught exception "+string(e.what())+" while taking exp of diagTransMatrix cell "+stringify(clab);
			throw runtime_error(errstr);
		}

		if (clab % nStates == 0) {
			for (QNUInt32 plab=0; plab<nFullLabs; plab++) {
				try {
					QNUInt32 idx=plab*nFullLabs+clab/this->nStates;
					this->denseTransMatrix[idx]=expE(this->denseTransMatrix[idx]);
				}
				catch (exception &e) {
					string errstr="CRF_HNStateNode::computeTransMatrix() caught exception "+string(e.what())+" while taking exp of denseTransMatrix cell "+stringify(plab)+", "+stringify(clab);
					throw runtime_error(errstr);
				}
			}
		}
		else {
			try {
				this->offDiagTransMatrix[clab-1]=expE(this->offDiagTransMatrix[clab-1]);
			}
			catch (exception &e) {
				string errstr="CRF_HNStateNode::computeTransMatrix() caught exception "+string(e.what())+" while taking exp of offDiagTransMatrix cell "+stringify(clab);
				throw runtime_error(errstr);
			}
		}
	}
	return result;
}

double CRF_HNStateNode::computeAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=0.0; // Reset alphaArray to zero before processing
	}

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		// First we compute the self transition - all labels get this
		this->alphaArray[clab]=prev_alpha[clab]*this->diagTransMatrix[clab];
		// Next we add in transitions from prior states
		if (clab % this->nStates == 0) {
			// Here clab is a new start state, so all end state transitions to it must be computed
			QNUInt32 dense_clab = clab/this->nStates; //Used to index into dense transition matrix
			for (QNUInt32 plab=0; plab<nFullLabs; plab++)
			{
				QNUInt32 real_prev = plab*nStates+nStates-1;
				QNUInt32 idx = plab*nFullLabs+dense_clab;
				this->alphaArray[clab]+=prev_alpha[real_prev]*this->denseTransMatrix[idx];
			}
		}
		else {
			// Here clab is an interior state, so only transitions from the previous state need to be
			// accounted for
			this->alphaArray[clab]+=prev_alpha[clab-1]*this->offDiagTransMatrix[clab-1];
		}
		// And of course we multiply by the state feature array at the end
	 	this->alphaArray[clab]*=this->stateArray[clab];
	 	this->alphaScale += this->alphaArray[clab];
	 	//alpha_tot+=new_alpha[clab];
	 }

	if ( (this->alphaScale<1.0) && (this->alphaScale>-1.0)) this->alphaScale=1.0; //Don't inflate scores
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=this->alphaArray[clab]/this->alphaScale;
	}


	return this->alphaScale;

}

double CRF_HNStateNode::computeFirstAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=0.0; // Reset alphaArray to zero before processing
	}

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
	 	this->alphaArray[clab]=this->stateArray[clab];
	 	this->alphaScale += this->alphaArray[clab];
	 	//alpha_tot+=new_alpha[clab];
	 }

	if ( (this->alphaScale<1.0) && (this->alphaScale>-1.0)) this->alphaScale=1.0; //Don't inflate scores
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=this->alphaArray[clab]/this->alphaScale;
	}

	return this->alphaScale;

}


double CRF_HNStateNode::computeBeta(double* result_beta, double scale)
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
		this->tempBeta[clab]=this->betaArray[clab]*this->stateArray[clab]/scale;
		result_beta[clab]=0;
	}
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, this->num_labs, 1, this->num_labs, 1.0f, Mtrans, this->num_labs, this->tmp_beta, this->num_labs, 1.0f, new_beta,1);
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLabs, 1, nLabs, 1.0f, this->transMatrix, nLabs, this->tempBeta, nLabs, 1.0f, result_beta,1);

	for (QNUInt32 plab=0; plab<nLabs; plab++) {
		// Here the logic is the inverse of that used in the computeAlphas above
		// because the loop driver here is the previous label

		// First we add the self transition to the result
		result_beta[plab]=this->tempBeta[plab]*this->diagTransMatrix[plab];

		// Check to see if plab is an end state.  If it isn't then we only need to worry about adding in
		// plab+1
		if ( (plab+1) % this->nStates == 0) {
			// plab is an end state, we have to compute for all possible start states
			QNUInt32 idx_plab = (plab+1)/this->nStates -1;
			for (QNUInt32 clab=0; clab<this->nFullLabs; clab++) {
				QNUInt32 idx=idx_plab*this->nFullLabs+clab;
				QNUInt32 real_clab = clab*nStates;
				result_beta[plab]+=this->denseTransMatrix[idx]*this->tempBeta[real_clab];
			}
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
				result_beta[plab]+=this->offDiagTransMatrix[plab]*this->tempBeta[plab+1];
			}
			else {
				string errstr="CRF_HNStateNode::computeBeta() found odd error that needs to be checked";
				throw runtime_error(errstr);
			}
		}
	}

	return this->alphaScale;
}

double* CRF_HNStateNode::computeAlphaBeta(double Zx)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaBetaArray[clab]=this->alphaArray[clab]*this->betaArray[clab]/Zx;
	}
	return this->alphaBetaArray;
}


void CRF_HNStateNode::setTailBeta()
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->betaArray[clab]=1.0/this->alphaScale;
	}
}

double CRF_HNStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	double logLi=0.0;
	double alpha_beta=0.0;
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		alpha_beta=this->alphaArray[clab]*this->betaArray[clab]*this->alphaScale/Zx;
		alpha_beta_tot += alpha_beta;
		//bool match=(clab==this->label);
		logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,this->label,clab);
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
			//for (QNUInt32 plab=0; plab<nLabs; plab++) {
			//	lc+=this->crf_ptr->getFeatureMap()->getNumTransFuncs(plab,clab);
			//	this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,NULL,lc,tmp_ExpF,NULL,alpha_beta,false,plab,clab);
			//}
			//delete tmp_ExpF;
		}
		else {
			// First we compute the self transitions
			alpha_beta=prev_alpha[clab]*this->diagTransMatrix[clab]*this->stateArray[clab]*this->betaArray[clab]/Zx;
			alpha_beta_trans_tot+=alpha_beta;
			//match=((clab==this->label)&&(clab==prev_lab));
			logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,clab,clab);
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
					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,real_plab,clab);
					//logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,real_plab,clab);
				}
			}
			else {
				// If not, we only update the lambda values on the offDiagonal
				alpha_beta=prev_alpha[clab-1]*this->offDiagTransMatrix[clab-1]*this->stateArray[clab]*this->betaArray[clab]/Zx;
				alpha_beta_trans_tot+=alpha_beta;
				//match=((clab==this->label)&&(clab-1==prev_lab));
				logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(this->ftrBuf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,clab-1,clab);
				//logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab-1,clab);
			}
		}
	}
	if ((alpha_beta_tot >1.1))  {
		string errstr="CRF_HNStateNode::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot < 0.9) {
		string errstr="CRF_HNStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.1) {
		string errstr="CRF_HNStateNode::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot < 0.9) {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;

}

double CRF_HNStateNode::computeAlphaSum()
{
	double Zx = 0.0;
	QNUInt32 nLabs=this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		Zx+=this->alphaArray[clab];
	}
	return Zx;
}

double CRF_HNStateNode::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
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
		return 0;
	}
}

double CRF_HNStateNode::getStateValue(QNUInt32 cur_lab)
{
	return this->stateArray[cur_lab];
}

double CRF_HNStateNode::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	if ((cur_lab % this->nStates == 0) && ((prev_lab+1) % this->nStates == 0)) {
		QNUInt32 cur_idx = cur_lab/this->nStates;
		QNUInt32 prev_idx = (prev_lab+1)/this->nStates -1;
		QNUInt32 idx=prev_idx*this->nFullLabs+cur_idx;
		return this->denseTransMatrix[idx]*this->stateArray[cur_lab];
	}
	else if (cur_lab == prev_lab) {
		return this->diagTransMatrix[cur_lab]*this->stateArray[cur_lab];
	}
	else if (cur_lab == prev_lab+1) {
		return this->offDiagTransMatrix[cur_lab-1]*this->stateArray[cur_lab];
	}
	else {
		return 0;
	}
}
*/
