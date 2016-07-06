/*
 * CRF_StdSegNStateNode.cpp
 *
 *  Created on: Apr 29, 2012
 *      Author: Yanzhang (Ryan) He
 */

#include "CRF_StdSegNStateNode.h"

// TODO: remove the match between the ranges of two dimensions of the transition matrix, so to delete prevNode_nLabs and nextNode_nActualLabs.
/*
 * CRF_StdSegNStateNode constructor
 *
 * Input: fb - feature buffer for features applicable for this node
 *        sizeof_fb - number of features in buffer fb
 *        lab - label for this node
 *        crf_in - pointer back to CRF model used by this node
 *        nodeMaxDur - the maximum duration of labels for this node
 *        prevNode_nLabs - the number of all labels for previous node, for use in defining the size of the transition matrix.
 *        nextNode_nActualLabs - the number of actual labels (without duration) for next node, for use in beta calculation for current node.
 */
CRF_StdSegNStateNode::CRF_StdSegNStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs)
	: CRF_StdNStateNode(fb, sizeof_fb, lab, crf)
{
	this->labMaxDur = crf_ptr->getLabMaxDur();
	this->nodeLabMaxDur = nodeMaxDur;
	if (this->nodeLabMaxDur > this->labMaxDur)
	{
		char errstr[1024];
		sprintf(errstr, "CRF_StdSegNStateNode constructor caught exception: the maximum duration of labels for this node (%lu) is larger than the maximum possible duration of labels for all nodes (%lu).",(unsigned long)this->nodeLabMaxDur,(unsigned long)this->labMaxDur);
		throw runtime_error(errstr);
	}
	// commented out by Ryan, for CRF_StdSegNStateNode_WithoutDurLab
//	if (this->nLabs % this->labMaxDur != 0) {
//		char errstr[1024];
//		sprintf(errstr, "CRF_StdSegNStateNode constructor caught exception: the number of all labels (%lu) and the maximum duration of labels (%lu)do not correspond.", (unsigned long)this->nLabs, (unsigned long)this->labMaxDur);
//		throw runtime_error(errstr);
//	}
	this->nActualLabs = nLabs/this->labMaxDur;
	this->numAvailLabs = this->nActualLabs * this->nodeLabMaxDur;

	for (QNUInt32 clab = 0; clab < nLabs; clab++)
	{
		this->alphaArray[clab] = CRF_LogMath::LOG0;
		this->betaArray[clab] = CRF_LogMath::LOG0;
		this->alphaBetaArray[clab] = CRF_LogMath::LOG0;
	}

	// TODO: verify sizeof_fb since the feature buffer now contains the features for more than one segments.
	if (this->ftrBuf_size % this->nodeLabMaxDur != 0) {
		char errstr[1024];
		sprintf(errstr, "CRF_StdSegNStateNode constructor caught exception: the size of the feature buffer (%lu) and the maximum duration of labels for this node (%lu) do not correspond.", (unsigned long)this->ftrBuf_size, (unsigned long)this->nodeLabMaxDur);
		throw runtime_error(errstr);
	}
	this->nFtrsPerSeg = this->ftrBuf_size / this->nodeLabMaxDur;

	// initialize the number of previous nodes and next nodes.
	// they have to be set up before alpha and beta calculation.
	this->numPrevNodes = CRF_UINT32_MAX;
	this->numNextNodes = CRF_UINT32_MAX;

	this->prevNodeNLabs = prevNode_nLabs;
	this->nextNodeNActualLabs = nextNode_nActualLabs;

	// TODO: remove the match between the ranges of two dimensions of the transition matrix
	// transMatrix created in the base class has same size of horizontal and vertical indice.
	// but we need them to be different for this class.
//	delete [] this->transMatrix;
//	this->transMatrix = new double[this->prevNodeNLabs * this->nLabs];
}

/*
 * CRF_StdSegNStateNode destructor
 */
CRF_StdSegNStateNode::~CRF_StdSegNStateNode() {
	// transMatrix will be deleted in the base class.
}

/*
 * CRF_StdSegNStateNode::computeTransMatrix
 *
 * Computes the log of the state vector (stateArray) and the log of the transition matrix (transMatrix)
 *   and stores them as appropriate
 */
double CRF_StdSegNStateNode::computeTransMatrix()
{
	string errstr="CRF_StdSegNStateNode::computeTransMatrix() threw exception: 3-state segmental models for phone-duration labels are not implemented yet!";
	throw runtime_error(errstr);

//	checkNumPrevNodes();
//
//	double result=0.0;
//
//	double* lambda = this->crf_ptr->getLambda();
//	float* seg_ftr_buf = this->ftrBuf;
//	QNUInt32 clab = 0;
//	// Note: it should be this->numPrevNodes instead of this->nodeLabMaxDur as the number of iterations of the outer loop.
//	// It's because numPrevNodes is enough for transition calculation while nodeLabMaxDur might be larger than numPrevNodes for some nodes.
//	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
//	{
//		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,clab);
//
//			// TODO: to delete this chunk of code, which is from 1-state segmental model and only for reference
///*
////			//TODO: for (QNUInt32 plab=0; plab<nLabs_of_prevNode; plab++) {
////			//for (QNUInt32 plab=0; plab<nLabs; plab++) {
////			for (QNUInt32 plab = 0; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {
////
////				QNUInt32 idx = plab * this->nLabs + clab;
////				//TODO: design a feature map in which the transition matrix calculation can take different dimensions of plab (from prevNode) and clab (from current node).
////				this->transMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(seg_ftr_buf,lambda,plab,clab);
////
////				// just for debugging
//////				cout << "[" << plab << "]=" << this->transMatrix[idx] << " ";
////			}
//*/
//
//
//			// Here we need to work some magic.  All entries on the diagonal get their self transition assigned
//			//this->diagTransMatrix[clab]=this->crf_ptr->getFeatureMap()->computeMij(this->ftrBuf,lambda,lc,clab,clab);
//			this->diagTransMatrix[clab]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,clab,clab);
//
//			// Besides the diagonal, we need to update the off-diagonal or the dense transition matrix
//			// dense transition matrix is updated when the current label is a start state (e.g. when the
//			// clab % nStates == 0
//			if (clab % this->nStates == 0) {
//				for (QNUInt32 plab=0; plab<nFullLabs; plab++) {
//					// Our index into the dense transition matrix needs to be munged a bit
//					QNUInt32 idx=plab*nFullLabs+clab/this->nStates;
//					// And our previous label is actually the end state for the previous "label"
//					QNUInt32 real_plab = plab*nStates+nStates-1;  // real_plab is the end state
//					this->denseTransMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,real_plab,clab);
//				}
//			}
//			else {
//				// We're on the off-diagonal - clab is not a start  state and the only previous
//				// transition must be from the immediate prior state (e.g. clab-1)
//				this->offDiagTransMatrix[clab-1]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,clab-1,clab);
//			}
//
//			clab++;
//
//		}
//		seg_ftr_buf += this->nFtrsPerSeg;
//	}
//	// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
//	for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
//	{
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,clab);
//			clab++;
//		}
//		seg_ftr_buf += this->nFtrsPerSeg;
//	}
//	return result;
}

/*
 * CRF_StdSegNStateNode::computeAlpha
 *
 * Read alpha vectors of previous nodes directly from prevNode and store the result of the alpha vector in alphaArray.
 *
 * Compute the alpha vector for the forward backward computation for this node.
 */
double CRF_StdSegNStateNode::computeAlpha()
{
	string errstr="CRF_StdSegNStateNode::computeAlpha() threw exception: 3-state segmental models for phone-duration labels are not implemented yet!";
	throw runtime_error(errstr);

//	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
//	this->alphaScale=0.0;
//
//	checkNumPrevNodes();
//
//	QNUInt32 clab = 0;
//	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
//	{
//		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
//		double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			// Next we add in transitions from prior states
//			// Compute the self transition - all labels get this
//			this->alphaArray[clab] = prev_adj_seg_alpha[clab] + this->diagTransMatrix[clab];
//
//			// TODO: to delete this chunk of code, which is from 1-state segmental model and only for reference
///*
////			this->logAddAcc[0]=prev_adj_seg_alpha[0]+this->transMatrix[0+clab];
////
////			double maxv=this->logAddAcc[0];
////			//TODO: for (QNUInt32 plab = 1; plab < prevAdjacentSeg->getNLabs(); plab++) {  //full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
////			for (QNUInt32 plab = 1; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {   //faster implementation, not guaranteed to work for all classes of previous nodes.
////				//this->logAddAcc[plab]=prev_adj_seg_alpha[plab]+this->transMatrix[plab * prevAdjacentSeg->getNLabs() + clab];
////				this->logAddAcc[plab]=prev_adj_seg_alpha[plab]+this->transMatrix[plab * this->nLabs + clab];
////
////				if (this->logAddAcc[plab]>maxv) {
////					maxv=logAddAcc[plab];
////				}
////			}
////			try {
////				//TODO: this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,prevAdjacentSeg->getNLabs());        //full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
////				this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,prevAdjacentSeg->getNumAvailLabs());   //faster implementation, not guaranteed to work for all classes of previous nodes.
////			}
////			catch (exception &e) {
////				string errstr="CRF_StdSegNStateNode::computeAlpha() caught exception: "+string(e.what())+" while computing alpha";
////				throw runtime_error(errstr);
////				return(-1);
////			}
//*/
//
//
//			try {
//				if (clab % this->nStates == 0) {
//					// Here clab is a new start state, so all end state transitions to it must be computed
//					QNUInt32 dense_clab = clab/this->nStates; //Used to index into dense transition matrix
//					this->logAddAcc[0] = prev_adj_seg_alpha[nStates-1] + this->denseTransMatrix[dense_clab];
//					double max = this->logAddAcc[0];
//					for (QNUInt32 plab = 1; plab < nFullLabs; plab++)
//					{
//						QNUInt32 real_prev = plab * nStates + nStates - 1;
//						QNUInt32 idx = plab * nFullLabs + dense_clab;
//						this->logAddAcc[plab] = prev_adj_seg_alpha[real_prev] + this->denseTransMatrix[idx];
//						if (this->logAddAcc[plab] > max) {
//							max = this->logAddAcc[plab];
//						}
//					}
//					double logSum = logAdd(this->logAddAcc,max,nFullLabs);
//					this->alphaArray[clab] = logAdd(this->alphaArray[clab],logSum);
//				}
//				else {
//					// Here clab is an interior state, so only transitions from the previous state need to be
//					// accounted for
//					this->alphaArray[clab] = logAdd(this->alphaArray[clab],prev_adj_seg_alpha[clab-1]+this->offDiagTransMatrix[clab-1]);
//				}
//			}
//			catch (exception &e) {
//				string errstr = "CRF_StdSegNStateNode::computeAlpha() caught exception: " + string(e.what()) + " while computing alpha";
//				throw runtime_error(errstr);
//				return(-1);
//			}
//			this->alphaArray[clab]+=this->stateArray[clab];
//			clab++;
//		}
//	}
//	// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
//	for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
//	{
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			this->alphaArray[clab]=this->stateArray[clab];
//			clab++;
//		}
//	}
//
//	return this->alphaScale;

}


/*
 * CRF_StdSegNStateNode::computeFirstAlpha
 *
 * Compute the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdSegNStateNode::computeFirstAlpha()
{
	string errstr="CRF_StdSegNStateNode::computeFirstAlpha() threw exception: 3-state segmental models for phone-duration labels are not implemented yet!";
	throw runtime_error(errstr);

//	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
//	this->alphaScale = 0.0;
//
//	//numAvailLabs for the first node of the sequence is usually equal to nActualLabs (since nodeLabMaxDur==1).
//	for (QNUInt32 clab = 0; clab < this->numAvailLabs; clab++)
//	{
//		this->alphaArray[clab] = this->stateArray[clab];
//	}
//	return this->alphaScale;

}


/*
 * CRF_StdSegNStateNode::computeBeta
 *
 * Inputs: scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Read the beta vectors of next nodes directly from nextNode and store the result of the beta vector in betaArray.
 *
 * Compute the beta vector for the node before this one and store it in result_beta
 */
double CRF_StdSegNStateNode::computeBeta(double scale)
{
	string errstr="CRF_StdSegNStateNode::computeBeta() threw exception: 3-state segmental models for phone-duration labels are not implemented yet!";
	throw runtime_error(errstr);

//	// Logic desired:
//	//	* Compute beta_i[size of alpha[]+1] to be all 1s
//	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]
//
//	checkNumNextNodes();
//
//	// if numNextNodes == 0, this is the last node of the sequence.
//	// Sets the beta value in this node to the special case for the end of the sequence.
//	if (this->numNextNodes == 0)
//	{
//		setTailBeta();
//		return this->alphaScale;
//	}
//
//	QNUInt32 nextlab = 0;
//	for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
//	{
//		CRF_StateNode* nextAdjacentSeg = this->nextNodes[dur - 1];
//		double* next_adj_seg_beta = nextAdjacentSeg->getBeta();
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++) {
//			this->tempBeta[nextlab] = next_adj_seg_beta[nextlab] + nextAdjacentSeg->getStateValue(nextlab);
//			nextlab++;
//		}
//	}
//
//	double* internal_state_logAddAcc = new double[this->numNextNodes * 2];
//
//	for (QNUInt32 clab = 0; clab < this->numAvailLabs; clab++)
//	{
//		CRF_StateNode* nextAdjacentSeg = this->nextNodes[0];
//		//this->transMatrix[plab*nLabs+0]+this->tempBeta[0]
//
//		// calculate the first max value in the internal_state_logAddAcc array
//		internal_state_logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, clab) + this->tempBeta[clab];
//		double internal_state_maxv = internal_state_logAddAcc[0];
//		int internal_state_logAddAcc_id = 0;
//
//		// calculate the first max value in the logAddAcc array (for cross phone transitions)
//		double cross_phone_maxv;
//		if ((clab + 1) % this->nStates == 0) {
//			this->logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, 0) + this->tempBeta[0];
//			cross_phone_maxv = this->logAddAcc[0];
//		}
//
//		QNUInt32 nextlab = 0;
//
//		for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
//		{
//			nextAdjacentSeg = this->nextNodes[dur - 1];
//
//			// First we add the self transition to the result
//			internal_state_logAddAcc[internal_state_logAddAcc_id] = nextAdjacentSeg->getTransValue(clab, clab) + this->tempBeta[clab];
//			if (internal_state_logAddAcc[internal_state_logAddAcc_id] > internal_state_maxv) {
//				internal_state_maxv = internal_state_logAddAcc[internal_state_logAddAcc_id];
//			}
//			internal_state_logAddAcc_id++;
//
//			// Check to see if clab is an end state.
//			// If it isn't then we only need to worry about adding in clab+1.
//			if ((clab + 1) % this->nStates == 0) {
//				//TODO: It should be nActualLabs_of_nextNode instead of nActualLabs_of_thisNode as the number of iterations.
//				//for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++) {
//				for (QNUInt32 lab = 0; lab < this->nextNodeNActualLabs; lab++) {
//
//					//this->logAddAcc[clab]=this->transMatrix[plab*nLabs+clab]+this->tempBeta[clab];
//					this->logAddAcc[nextlab] = nextAdjacentSeg->getTransValue(clab, nextlab) + this->tempBeta[nextlab];
//
//					if (this->logAddAcc[nextlab] > cross_phone_maxv) {
//						cross_phone_maxv = logAddAcc[nextlab];
//					}
//
//					nextlab++;
//				}
//			}
//			else {
//				// clab isn't an end state - we only have to compute for the next state (i.e. clab+1)
//				if (clab != nLabs-1) {
//					// Sanity check, clab should never be the last label in the sequence (that should be
//					// an end state).  Just in case
//					internal_state_logAddAcc[internal_state_logAddAcc_id] = nextAdjacentSeg->getTransValue(clab, clab + 1) + this->tempBeta[clab + 1];
//					if (internal_state_logAddAcc[internal_state_logAddAcc_id] > internal_state_maxv) {
//						internal_state_maxv = internal_state_logAddAcc[internal_state_logAddAcc_id];
//					}
//					internal_state_logAddAcc_id++;
//				}
//				else {
//					string errstr="CRF_StdSegNStateNode::computeBeta() found odd error that needs to be checked";
//					throw runtime_error(errstr);
//				}
//			}
//		}
//		try {
//			this->betaArray[clab]=logAdd(internal_state_logAddAcc,internal_state_maxv,internal_state_logAddAcc_id);
//
//			if ((clab + 1) % this->nStates == 0) {
//				//TODO:It should be nActualLabs_of_nextNode instead of nActualLabs_of_thisNode.
//				//this->betaArray[clab]=logAdd(this->logAddAcc,maxv,this->nActualLabs * this->numNextNodes);
//				double logSum = logAdd(this->logAddAcc,cross_phone_maxv,this->nextNodeNActualLabs * this->numNextNodes);
//				this->betaArray[clab] = logAdd(this->betaArray[clab],logSum);
//			}
//		}
//		catch (exception &e) {
//			string errstr="CRF_StdSegNStateNode::computeBeta() caught exception: "+string(e.what())+" while computing beta";
//			throw runtime_error(errstr);
//			return(-1);
//		}
//	}
//
//	delete [] internal_state_logAddAcc;
//
//	return this->alphaScale;
}

/*
 * CRF_StdSegNStateNode::computeExpF
 *
 * Inputs: *ExpF - vector to store expected values of feature functions
 *         *grad - vector to store computed gradient values
 *         Zx - normalization constant
 *         prev_lab - previous node label (transition feature ExpF computation)
 *
 * Returns:
 *
 * Read alpha vectors of previous nodes directly from prevNode for use in transition feature ExpF computation.
 *
 * Compute gradient and expected values for features in this node and store them in *grad and
 *   *ExpF vectors respectively.  State features and transition features are computed in the same function.
 */
double CRF_StdSegNStateNode::computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 prev_lab)
{
	string errstr="CRF_StdSegNStateNode::computeExpF() threw exception: 3-state segmental models for phone-duration labels are not implemented yet!";
	throw runtime_error(errstr);
}

/*
 * CRF_StdSegNStateNode::computeAlphaBeta
 *
 * Inputs: Zx - normalization constant
 *
 * Returns: array of state probabilities based on the precomputed alpha and beta for this node
 *
 */
double* CRF_StdSegNStateNode::computeAlphaBeta(double Zx)
{
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab = 0; clab < this->numAvailLabs; clab++)
	{
		this->alphaBetaArray[clab]=this->alphaArray[clab]+this->betaArray[clab]-Zx;
	}
	return this->alphaBetaArray;
}

/*
 * CRF_StdSegNStateNode::computeAlphaSum
 *
 * Returns: Sum of the values in the alpha vector of this node
 *
 * Used to compute the normalization constant for the CRF.
 */
double CRF_StdSegNStateNode::computeAlphaSum()
{
	double Zx;
	try {
		//Zx=logAdd(this->alphaArray,this->nLabs);  // full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
		Zx=logAdd(this->alphaArray,this->numAvailLabs);
	}
	catch (exception& e) {
		string errstr="CRF_StdSegNStateNode::computeAlphaSum() threw exception: "+string(e.what());
		throw runtime_error(errstr);
	}
	return Zx;
}

/*
 * CRF_StdSegNStateNode::computeAlphaAlignedSum
 *
 * Returns: Sum of the values in the alphaAligned vector of this node
 *
 * Used to compute the "soft" normalization constant.
 */
double CRF_StdSegNStateNode::computeAlphaAlignedSum()
{
	double Zx;
	try {
		//Zx=logAdd(&(this->alphaArrayAligned),this->nLabs);  // full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
		Zx=logAdd(&(this->alphaArrayAligned),this->numAvailLabs);
	}
	catch (exception& e) {
		string errstr="CRF_StdSegNStateNode::computeAlphaAlignedSum() threw exception: "+string(e.what());
		throw runtime_error(errstr);
	}
	return Zx;
}

/*
 * CRF_StdSegNStateNode::getNActualLabs
 *
 * Returns: the number of actual labels (without duration).
 *
 * Accessor function for nActualLabs
 */
QNUInt32 CRF_StdSegNStateNode::getNActualLabs()
{
	return this->nActualLabs;
}

/*
 * CRF_StdSegNStateNode::checkNumPrevNodes
 *
 * check if numPrevNodes has been assigned and <= nodeLabMaxDur
 *
 */
bool CRF_StdSegNStateNode::checkNumPrevNodes()
{
	if (this->numPrevNodes == CRF_UINT32_MAX)
	{
		string errstr="CRF_StdSegNStateNode::checkNumPrevNodes() caught exception: the vector of previous nodes and the number of them have not been set yet.";
		throw runtime_error(errstr);
	}
	if (this->numPrevNodes > this->nodeLabMaxDur)
	{
		string errstr="CRF_StdSegNStateNode::checkNumPrevNodes() caught exception: the number of previous nodes is larger than the maximum label duration of the current node.";
		throw runtime_error(errstr);
	}
	return true;
}

/*
 * CRF_StdSegNStateNode::checkNumNextNodes
 *
 * check if numNextNodes has been assigned and <= labMaxDur
 *
 */
bool CRF_StdSegNStateNode::checkNumNextNodes()
{
	if (this->numNextNodes == CRF_UINT32_MAX)
	{
		string errstr="CRF_StdSegNStateNode::checkNumNextNodes() caught exception: the vector of following nodes and the number of them have not been set yet.";
		throw runtime_error(errstr);
	}
	if (this->numNextNodes > this->labMaxDur)
	{
		string errstr="CRF_StdSegNStateNode::checkNumNextNodes() caught exception: the number of following nodes is larger than the maximum label duration of all nodes.";
		throw runtime_error(errstr);
	}
	return true;
}

// Override these functions which do not work for CRF_StdSegNStateNode

/*
 * CRF_StdSegNStateNode::computeAlpha
 *
 * Input: prev_alpha - alpha vector of the previous node
 *
 * Computes the alpha vector for the forward backward computation for this node.
 */
double CRF_StdSegNStateNode::computeAlpha(double* prev_alpha)
{
	string errstr="CRF_StdSegNStateNode::computeAlpha(double* prev_alpha) threw exception: CRF_StdSegNStateNode::computeAlpha() should be used instead.";
	throw runtime_error(errstr);
}

/*
 * CRF_StdSegNStateNode::computeFirstAlpha
 *
 * Computes the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdSegNStateNode::computeFirstAlpha(double* prev_alpha)
{
	string errstr="CRF_StdSegNStateNode::computeFirstAlpha(double* prev_alpha) threw exception: CRF_StdSegNStateNode::computeFirstAlpha() should be used instead.";
	throw runtime_error(errstr);
}

/*
 * CRF_StdSegNStateNode::computeBeta
 *
 * Inputs: result_beta - vector to store the result of the computation
 *         scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Computes the beta vector for the node before this one and store it in result_beta
 */
double CRF_StdSegNStateNode::computeBeta(double* result_beta, double scale)
{
	string errstr="CRF_StdSegNStateNode::computeBeta(double* result_beta, double scale) threw exception: CRF_StdSegNStateNode::computeBeta(double scale) should be used instead.";
	throw runtime_error(errstr);
}

/*
 * CRF_StdSegNStateNode::computeExpF
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
double CRF_StdSegNStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	string errstr="CRF_StdSegNStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab) threw exception: CRF_StdSegNStateNode::computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 prev_lab) should be used instead.";
	throw runtime_error(errstr);
}
