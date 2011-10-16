/*
 * CRF_StdSegStateNode.cpp
 *
 *  Created on: Sep 13, 2011
 *      Author: hey
 */

#include "CRF_StdSegStateNode.h"

/*
 * CRF_StdSegStateNode constructor
 *
 * Input: fb - feature buffer for features applicable for this node
 *        sizeof_fb - number of features in buffer fb
 *        lab - label for this node
 *        crf_in - pointer back to CRF model used by this node
 *        nodeMaxDur - the maximum duration of labels for this node
 *        prevNode_nLabs - the number of all labels for previous node, for use in defining the size of the transition matrix.
 *        nextNode_nActualLabs - the number of actual labels (without duration) for next node, for use in beta calculation for current node.
 */
CRF_StdSegStateNode::CRF_StdSegStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs)
	: CRF_StdStateNode(fb, sizeof_fb, lab, crf)
{
	// just for debugging
//	cout << "Creating a CRF_StdSegStateNode node:" << endl;
//	cout << "sizeof_fb=" << sizeof_fb << endl;
//	cout << "lab=" << lab << endl;
//	cout << "nLabs=" << crf_ptr->getNLabs() << endl;
//	cout << "labMaxDur=" << crf_ptr->getLabMaxDur() << endl;
//	cout << "nodeLabMaxDur=" << nodeMaxDur << endl;
//	cout << "prevNode_nLabs=" << prevNode_nLabs << endl;
//	cout << "nextNode_nActualLabs=" << nextNode_nActualLabs << endl;

	this->labMaxDur = crf_ptr->getLabMaxDur();
	this->nodeLabMaxDur = nodeMaxDur;
	if (this->nodeLabMaxDur > this->labMaxDur)
	{
		char errstr[1024];
		sprintf(errstr, "CRF_StdSegStateNode constructor caught exception: the maximum duration of labels for this node (%lu) is larger than the maximum possible duration of labels for all nodes (%lu).",(unsigned long)this->nodeLabMaxDur,(unsigned long)this->labMaxDur);
		throw runtime_error(errstr);
	}
	// commented out by Ryan, for CRF_StdSegStateNode_WithoutDurLab
//	if (this->nLabs % this->labMaxDur != 0) {
//		char errstr[1024];
//		sprintf(errstr, "CRF_StdSegStateNode constructor caught exception: the number of all labels (%lu) and the maximum duration of labels (%lu)do not correspond.", (unsigned long)this->nLabs, (unsigned long)this->labMaxDur);
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
		sprintf(errstr, "CRF_StdSegStateNode constructor caught exception: the size of the feature buffer (%lu) and the maximum duration of labels for this node (%lu) do not correspond.", (unsigned long)this->ftrBuf_size, (unsigned long)this->nodeLabMaxDur);
		throw runtime_error(errstr);
	}
	this->nFtrsPerSeg = this->ftrBuf_size / this->nodeLabMaxDur;

	// initialize the number of previous nodes and next nodes.
	// they have to be set up before alpha and beta calculation.
	this->numPrevNodes = CRF_UINT32_MAX;
	this->numNextNodes = CRF_UINT32_MAX;

	this->prevNodeNLabs = prevNode_nLabs;
	this->nextNodeNActualLabs = nextNode_nActualLabs;

	// transMatrix created in the base class has same size of horizontal and vertical indice.
	// but we need them to be different for this class.
	delete [] this->transMatrix;
	this->transMatrix = new double[this->prevNodeNLabs * this->nLabs];

	// just for debugging
//	cout << "nActualLabs=" << this->nActualLabs << endl;
//	cout << "numAvailLabs=" << this->numAvailLabs << endl;
//	cout << "nFtrsPerSeg=" << this->nFtrsPerSeg << endl;
//	cout << "transMatrix_size=" << this->prevNodeNLabs << "*" << this->nLabs << endl;
}

/*
 * CRF_StdSegStateNode destructor
 */
CRF_StdSegStateNode::~CRF_StdSegStateNode() {
	// transMatrix will be deleted in the base class.
}

/*
 * CRF_StdSegStateNode::computeTransMatrix
 *
 * Computes the log of the state vector (stateArray) and the log of the transition matrix (transMatrix)
 *   and stores them as appropriate
 */
double CRF_StdSegStateNode::computeTransMatrix()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode::computeTransMatrix(): " << endl;

	checkNumPrevNodes();

	double result=0.0;

	double* lambda = this->crf_ptr->getLambda();
	float* seg_ftr_buf = this->ftrBuf;
	QNUInt32 clab = 0;
	// Note: it should be this->numPrevNodes instead of this->nodeLabMaxDur as the number of iterations of the outer loop.
	// It's because numPrevNodes is enough for transition calculation while nodeLabMaxDur might be larger than numPrevNodes for some nodes.
	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
	{
		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,clab);

			// just for debugging
//			cout << "dur=" << dur << ", lab=" << lab << ", clab=" << clab << endl << "TransMatrix: plab=";
//			cout << "clab=" << clab << endl << "TransMatrix: plab=";

			//TODO: for (QNUInt32 plab=0; plab<nLabs_of_prevNode; plab++) {
			//for (QNUInt32 plab=0; plab<nLabs; plab++) {
			for (QNUInt32 plab = 0; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {

				QNUInt32 idx = plab * this->nLabs + clab;
				//TODO: design a feature map in which the transition matrix calculation can take different dimensions of plab (from prevNode) and clab (from current node).
				this->transMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(seg_ftr_buf,lambda,plab,clab);

				// just for debugging
//				cout << "[" << plab << "]=" << this->transMatrix[idx] << " ";
			}

			// just for debugging
//			cout << endl;
//			cout << "StateArray[" << clab << "]=" << this->stateArray[clab] << endl;

			clab++;

		}
		seg_ftr_buf += this->nFtrsPerSeg;

		// just for debugging
//		cout << "seg_ftr_buf moved forward by nFtrsPerSeg(" << nFtrsPerSeg << ")" << endl;
	}
	// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
	for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
	{
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,clab);

			// just for debugging
//			cout << "dur=" << dur << ", lab=" << lab << ", clab=" << clab << endl << "plab=";
			//cout << "clab=" << clab << endl << "TransMatrix: plab=";
//			cout << endl;
//			cout << "StateArray[" << clab << "]=" << this->stateArray[clab] << endl;

			clab++;
		}
		seg_ftr_buf += this->nFtrsPerSeg;
	}
	return result;
}

/*
 * CRF_StdSegStateNode::computeAlpha
 *
 * Read alpha vectors of previous nodes directly from prevNode and store the result of the alpha vector in alphaArray.
 *
 * Stub function.
 * Should compute the alpha vector for the forward backward computation for this node.
 */
double CRF_StdSegStateNode::computeAlpha()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode::computeAlpha(): " << endl;

	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	checkNumPrevNodes();

	QNUInt32 clab = 0;
	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
	{
		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
		double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			this->logAddAcc[0]=prev_adj_seg_alpha[0]+this->transMatrix[0+clab];

			// just for debugging
//			cout << "dur=" << dur << ", lab=" << lab << ", clab=" << clab << endl << "alphaArray = logAdd(prev alpha + trans): plab=[0]=" << this->logAddAcc[0] << " ";

			double maxv=this->logAddAcc[0];
			//TODO: for (QNUInt32 plab = 1; plab < prevAdjacentSeg->getNLabs(); plab++) {  //full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
			for (QNUInt32 plab = 1; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {   //faster implementation, not guaranteed to work for all classes of previous nodes.
				//this->logAddAcc[plab]=prev_adj_seg_alpha[plab]+this->transMatrix[plab * prevAdjacentSeg->getNLabs() + clab];
				this->logAddAcc[plab]=prev_adj_seg_alpha[plab]+this->transMatrix[plab * this->nLabs + clab];

				// just for debugging
//				cout << "[" << plab << "]=" << this->logAddAcc[plab] << " ";

				if (this->logAddAcc[plab]>maxv) {
					maxv=logAddAcc[plab];
				}
			}
			try {
				//TODO: this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,prevAdjacentSeg->getNLabs());        //full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
				this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,prevAdjacentSeg->getNumAvailLabs());   //faster implementation, not guaranteed to work for all classes of previous nodes.
			}
			catch (exception &e) {
				string errstr="CRF_StdSegStateNode::computeAlpha() caught exception: "+string(e.what())+" while computing alpha";
				throw runtime_error(errstr);
				return(-1);
			}
			// just for debugging
//			cout << endl << "alphaArray[" << clab << "](" << this->alphaArray[clab] << ") + stateArray[" << clab << "](" << this->stateArray[clab] << ") =";

			this->alphaArray[clab]+=this->stateArray[clab];

			// just for debugging
//			cout << this->alphaArray[clab] << endl;

			clab++;
		}
	}
	// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
	for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
	{
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			this->alphaArray[clab]=this->stateArray[clab];

			// just for debugging
//			cout << "dur=" << dur << ", lab=" << lab << ", clab=" << clab << endl;
//			cout << "alphaArray[" << clab << "] added stateArray[" << clab << "]=" << this->stateArray[clab] << endl;

			clab++;
		}
	}

	return this->alphaScale;

}

/*
 * CRF_StateNode::computeFirstAlpha
 *
 * Stub function.
 * Should compute the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdSegStateNode::computeFirstAlpha()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode::computeFirstAlpha()" << endl;

	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	//nodeMaxLab for the first node of the sequence is usually equal to nActualLabs (since nodeLabMaxDur==1).
	for (QNUInt32 clab = 0; clab < this->numAvailLabs; clab++)
	{
		// just for debugging
//		cout << "alphaArray[" << clab << "]=stateArray[" << clab << "]=" << this->stateArray[clab] << " ";

		this->alphaArray[clab]=this->stateArray[clab];
	}

	// just for debugging
//	cout << endl;

	return this->alphaScale;

}

/*
 * CRF_StateNode::computeBeta
 *
 * Inputs: scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Read the beta vectors of next nodes directly from nextNode and store the result of the beta vector in betaArray.
 *
 * Stub function.
 * Should compute the beta vector for the node before this one and store it in result_beta
 */
double CRF_StdSegStateNode::computeBeta(double scale)
{
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]

//	QNUInt32 nLabs = this->crf_ptr->getNLabs();
//
//	for (QNUInt32 clab=0; clab<nLabs; clab++) {
//		this->tempBeta[clab]=this->betaArray[clab]+this->stateArray[clab];
//	}
//
//	for (QNUInt32 plab=0; plab<nLabs; plab++) {
//		this->logAddAcc[0]=this->transMatrix[plab*nLabs+0]+this->tempBeta[0];
//		double maxv=this->logAddAcc[0];
//		for (QNUInt32 clab=1; clab<nLabs; clab++) {
//			this->logAddAcc[clab]=this->transMatrix[plab*nLabs+clab]+this->tempBeta[clab];
//			if (this->logAddAcc[clab]>maxv) {
//				maxv=this->logAddAcc[clab];
//			}
//		}
//		try {
//			result_beta[plab]=logAdd(this->logAddAcc,maxv,nLabs);
//		}
//		catch (exception &e) {
//			string errstr="CRF_StdSegStateNode::computeBeta() caught exception: "+string(e.what())+" while computing beta";
//			throw runtime_error(errstr);
//			return(-1);
//		}
//	}

	// just for debugging
//	cout << "CRF_StdSegStateNode::computeBeta(): " << endl;

	checkNumNextNodes();

	// if numNextNodes == 0, this is the last node of the sequence.
	// Sets the beta value in this node to the special case for the end of the sequence.
	if (this->numNextNodes == 0)
	{
		setTailBeta();
		return this->alphaScale;
	}

	QNUInt32 nextlab = 0;
	for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
	{
		CRF_StateNode* nextAdjacentSeg = this->nextNodes[dur - 1];
		double* next_adj_seg_beta = nextAdjacentSeg->getBeta();
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++) {
			this->tempBeta[nextlab] = next_adj_seg_beta[nextlab] + nextAdjacentSeg->getStateValue(nextlab);

			// just for debugging
//			cout << "tempBeta[" << nextlab << "]=" << "next_seg_beta[" << nextlab << "](" << next_adj_seg_beta[nextlab] << ") + next_seg_state_value[" << nextlab << "](" << nextAdjacentSeg->getStateValue(nextlab) << ")=" << this->tempBeta[nextlab] << endl;

			nextlab++;
		}
	}

	for (QNUInt32 clab = 0; clab < this->numAvailLabs; clab++)
	{
		CRF_StateNode* nextAdjacentSeg = this->nextNodes[0];
		//this->transMatrix[plab*nLabs+0]+this->tempBeta[0]
		this->logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, 0) + this->tempBeta[0];
		double maxv=this->logAddAcc[0];
		QNUInt32 nextlab = 0;

		// just for debugging
//		cout << "clab=" << clab << ", betaArray = logAdd(trans + tempBeta):" << endl;

		for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
		{
			nextAdjacentSeg = this->nextNodes[dur - 1];
			//TODO: It should be nActualLabs_of_nextNode instead of nActualLabs_of_thisNode as the number of iterations.
			//for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++) {
			for (QNUInt32 lab = 0; lab < this->nextNodeNActualLabs; lab++) {

				//this->logAddAcc[clab]=this->transMatrix[plab*nLabs+clab]+this->tempBeta[clab];
				this->logAddAcc[nextlab] = nextAdjacentSeg->getTransValue(clab, nextlab) + this->tempBeta[nextlab];

				// just for debugging
//				cout << " nextlab=[" << nextlab << "]=" << "trans(" << nextAdjacentSeg->getTransValue(clab, nextlab) << ") + tempBeta(" << this->tempBeta[nextlab] << ")="<< this->logAddAcc[nextlab];

				if (this->logAddAcc[nextlab]>maxv) {
					maxv=logAddAcc[nextlab];
				}
				nextlab++;
			}
		}
		try {
			//TODO:It should be nActualLabs_of_nextNode instead of nActualLabs_of_thisNode.
			//this->betaArray[clab]=logAdd(this->logAddAcc,maxv,this->nActualLabs * this->numNextNodes);
			this->betaArray[clab]=logAdd(this->logAddAcc,maxv,this->nextNodeNActualLabs * this->numNextNodes);

			// just for debugging
//			cout << endl << "betaArray[" << clab << "] = logAdd(trans + tempBeta) = " << this->betaArray[clab] << endl;
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
	for (QNUInt32 clab = 0; clab < this->numAvailLabs; clab++)
	{
		this->alphaBetaArray[clab]=this->alphaArray[clab]+this->betaArray[clab]-Zx;
	}
	return this->alphaBetaArray;
}

/*
 * CRF_StateNode::computeExpF
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
 * Stub function.
 * Should compute gradient and expected values for features in this node and store them in *grad and
 *   *ExpF vectors respectively.  State features and transition features are computed in the same function.
 */
double CRF_StdSegStateNode::computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 prev_lab)
{
	// Added by Ryan, just for debugging
//	cout << "Very beginning!!" << endl;

	checkNumPrevNodes();

	double logLi=0.0;
	double alpha_beta=0.0;
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;
	float* seg_ftr_buf = this->ftrBuf;
	QNUInt32 clab = 0;
	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
	{
		// just for debugging
//		cout << "dur: " << dur << endl;

		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
		double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();

		// just for debugging
//		cout << "After getting prevAdjacentSeg." << endl;

		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			// just debugging
//			cout << "First phase:: numPrevNodes: " << numPrevNodes << ", dur: " << dur << ", lab: " << lab << ", clab: " << clab << endl;

			alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
			alpha_beta_tot += alpha_beta;
			bool match=(clab==this->label);

			// just for debugging
//			if (match)
//				cout << clab << " ";
//				cout << "state label match: " << clab << " " << endl;

			logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,this->label,clab);
//			//TODO: verify: this->nLabs or nLabs_of_prevNode (probably this)
//			if (prev_lab > this->nLabs) {
//				// if prev_lab > nLabs, we're in the first label frame and there are no previous
//				// transitions - skip the transition calculation in this case
//				// but set the alpha_beta_trans_tot to 1.0 for the check below
//				alpha_beta_trans_tot=1.0;
//			}
//			else {
				// Otherwise do the transition calculations
				//for (QNUInt32 plab=0; plab<nLabs; plab++) {
				for (QNUInt32 plab = 0; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {
					QNUInt32 idx = plab * this->nLabs + clab;
					alpha_beta=expE(prev_adj_seg_alpha[plab]+this->transMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
					alpha_beta_trans_tot+=alpha_beta;
					match=((clab==this->label)&&(plab==prev_lab));

					// just debugging
//					if (match)
//						cout << "(" << clab << "<-" <<  plab << ") ";
//						cout << "transition label match: " << "(" << clab << "<-" <<  plab << ") " << endl;

					// TODO: design a feature map in which the transition matrix calculation can take different dimensions of plab (from prevNode) and clab (from current node).
					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab);
				}
//			}
			clab++;
		}
		seg_ftr_buf += this->nFtrsPerSeg;
	}
	// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
	for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
	{
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			// just debugging
//			cout << "Second phase:: numPrevNodes: " << numPrevNodes << ", dur: " << dur << ", lab: " << lab << ", clab: " << clab << endl;

			alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
			alpha_beta_tot += alpha_beta;
			bool match=(clab==this->label);

			// just debugging
//			if (match)
//				cout << clab << " ";
//				cout << "state label match: " << clab << " " << endl;


			logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,this->label,clab);

			clab++;
		}
		seg_ftr_buf += this->nFtrsPerSeg;
	}

	if (this->numPrevNodes == 0)
	{
		// if numPrevNodes == 0, there are no previous transitions
		// - skip the transition calculation in this case
		// but set the alpha_beta_trans_tot to 1.0 for the check below
		alpha_beta_trans_tot=1.0;
	}

	//just for debugging
//	cout << "\tAlpha_beta_tot: " << alpha_beta_tot << "\tAlpha_beta_trans_tot: " << alpha_beta_trans_tot << endl;

	// alpha_beta_tot and alpha_beta_trans_tot are no longer equal to 1 but less than 1 except the ending node.
	if ((alpha_beta_tot >1.000001))  {
		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
//	else if (alpha_beta_tot < 0.999999) {
	else if (alpha_beta_tot < -0.000001) {
//		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Probability sums less than 0.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.000001) {
		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
//	else if (alpha_beta_trans_tot < 0.999999) {
	else if (alpha_beta_trans_tot < -0.000001) {
//		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Trans Probability sums less than 0.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;
}

/*
 * CRF_StdSegStateNode::computeAlphaSum
 *
 * Returns: Sum of the values in the alpha vector of this node
 *
 * Used to compute the normalization constant for the CRF.
 */
double CRF_StdSegStateNode::computeAlphaSum()
{
	double Zx;
	try {

		// just for debugging
//		cout << "before compute alpha sum" << endl;

		//Zx=logAdd(this->alphaArray,this->nLabs);  // full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
		Zx=logAdd(this->alphaArray,this->numAvailLabs);

		// just for debugging
//		cout << "Zx=logAdd(this->alphaArray, " << this->numAvailLabs << ")=" << Zx << endl;
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
 * CRF_StdSegStateNode::computeAlphaAlignedSum
 *
 * Returns: Sum of the values in the alphaAligned vector of this node
 *
 * Used to compute the "soft" normalization constant.
 */
double CRF_StdSegStateNode::computeAlphaAlignedSum()
{
	double Zx;
	try {
		//Zx=logAdd(&(this->alphaArrayAligned),this->nLabs);  // full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
		Zx=logAdd(&(this->alphaArrayAligned),this->numAvailLabs);
	}
	catch (exception& e) {
		//changed by Ryan
		//string errstr="CRF_StdStateNodeLog::computeExpF() threw exception: "+string(e.what());
		string errstr="CRF_StdStateNode::computeAlphaAlignedSum() threw exception: "+string(e.what());
		throw runtime_error(errstr);
	}
	return Zx;
}

/*
 * CRF_StdSegStateNode::getNActualLabs
 *
 * Returns: the number of actual labels (without duration).
 *
 * Accessor function for nActualLabs
 */
QNUInt32 CRF_StdSegStateNode::getNActualLabs()
{
	return this->nActualLabs;
}

/*
 * CRF_StdSegStateNode::checkNumPrevNodes
 *
 * check if numPrevNodes has been assigned and <= nodeLabMaxDur
 *
 */
bool CRF_StdSegStateNode::checkNumPrevNodes()
{
	if (this->numPrevNodes == CRF_UINT32_MAX)
	{
		string errstr="CRF_StdSegStateNode::checkNumPrevNodes() caught exception: the vector of previous nodes and the number of them have not been set yet.";
		throw runtime_error(errstr);
	}
	if (this->numPrevNodes > this->nodeLabMaxDur)
	{
		string errstr="CRF_StdSegStateNode::checkNumPrevNodes() caught exception: the number of previous nodes is larger than the maximum label duration of the current node.";
		throw runtime_error(errstr);
	}
	return true;
}

/*
 * CRF_StdSegStateNode::checkNumNextNodes
 *
 * check if numNextNodes has been assigned and <= labMaxDur
 *
 */
bool CRF_StdSegStateNode::checkNumNextNodes()
{
	if (this->numNextNodes == CRF_UINT32_MAX)
	{
		string errstr="CRF_StdSegStateNode::checkNumNextNodes() caught exception: the vector of following nodes and the number of them have not been set yet.";
		throw runtime_error(errstr);
	}
	if (this->numNextNodes > this->labMaxDur)
	{
		string errstr="CRF_StdSegStateNode::checkNumNextNodes() caught exception: the number of following nodes is larger than the maximum label duration of all nodes.";
		throw runtime_error(errstr);
	}
	return true;
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
	string errstr="CRF_StdSegStateNode::computeAlpha(double* prev_alpha) threw exception: CRF_StdSegStateNode::computeAlpha() should be used instead.";
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
	string errstr="CRF_StdSegStateNode::computeFirstAlpha(double* prev_alpha) threw exception: CRF_StdSegStateNode::computeFirstAlpha() should be used instead.";
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
	string errstr="CRF_StdSegStateNode::computeBeta(double* result_beta, double scale) threw exception: CRF_StdSegStateNode::computeBeta(double scale) should be used instead.";
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
	string errstr="CRF_StdSegStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab) threw exception: CRF_StdSegStateNode::computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 prev_lab) should be used instead.";
	throw runtime_error(errstr);
}
