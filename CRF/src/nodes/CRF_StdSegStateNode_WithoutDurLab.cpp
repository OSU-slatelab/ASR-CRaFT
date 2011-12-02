/*
 * CRF_StdSegStateNode_WithoutDurLab.cpp
 *
 *  Created on: Sep 29, 2011
 *      Author: hey
 */

#include "CRF_StdSegStateNode_WithoutDurLab.h"

CRF_StdSegStateNode_WithoutDurLab::CRF_StdSegStateNode_WithoutDurLab(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs)
: CRF_StdSegStateNode(fb, sizeof_fb, lab, crf, nodeMaxDur, prevNode_nLabs, nextNode_nActualLabs) {

	this->nActualLabs = nLabs;
	this->numAvailLabs = this->nActualLabs;

	this->alphaArray_WithDur = new double[this->nActualLabs * this->nodeLabMaxDur];
//	for (QNUInt32 id = 0; id < this->nActualLabs * this->nodeLabMaxDur; id++)
//	{
//		this->alphaArray_WithDur[id] = CRF_LogMath::LOG0;
//	}

	delete [] this->stateArray;
	this->stateArray = new double[this->nActualLabs * this->nodeLabMaxDur];
//	for (QNUInt32 id = 0; id < this->nActualLabs * this->nodeLabMaxDur; id++)
//	{
//		this->stateArray[id] = 0.0;
//	}

	delete [] this->transMatrix;
	this->transMatrix = new double[this->nActualLabs * this->nActualLabs * this->nodeLabMaxDur];
//	for (QNUInt32 id = 0; id < this->nActualLabs * this->nActualLabs * this->nodeLabMaxDur; id++)
//	{
//		this->transMatrix[id] = 0.0;
//	}

	delete [] this->tempBeta;
	this->tempBeta = new double[this->nActualLabs * this->labMaxDur];  // this tempBeta is beta(phone, dur, endtime)*stateValue(phone, dur, endtime) for current node, which is different from tempBeta in CRF_StdSegStateNode.
//	for (QNUInt32 id = 0; id < this->nActualLabs * this->labMaxDur; id++)
//	{
//		this->tempBeta[id] = CRF_LogMath::LOG0;
//	}

	delete [] this->logAddAcc;
	this->logAddAcc = new double[this->nActualLabs * this->labMaxDur];
//	for (QNUInt32 id = 0; id < this->nActualLabs * this->labMaxDur; id++)
//	{
//		this->logAddAcc[id] = CRF_LogMath::LOG0;
//	}

}

CRF_StdSegStateNode_WithoutDurLab::~CRF_StdSegStateNode_WithoutDurLab() {

	delete [] this->alphaArray_WithDur;

	// all other arrays are deleted in the base classes.
}

/*
 * CRF_StdSegStateNode_WithoutDurLab::computeTransMatrix
 *
 * Computes the log of the state vector (stateArray) and the log of the transition matrix (transMatrix)
 *   and stores them as appropriate
 */
double CRF_StdSegStateNode_WithoutDurLab::computeTransMatrix()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab::computeTransMatrix(): " << endl;

	checkNumPrevNodes();

	double result=0.0;

	double* lambda = this->crf_ptr->getLambda();
	float* seg_ftr_buf = this->ftrBuf;
//	QNUInt32 clab = 0;
	// Note: it should be this->numPrevNodes instead of this->nodeLabMaxDur as the number of iterations of the outer loop.
	// It's because numPrevNodes is enough for transition calculation while nodeLabMaxDur might be larger than numPrevNodes for some nodes.
	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
	{
		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
		double* transMatrixForCurDur = &(this->transMatrix[this->nActualLabs * this->nActualLabs * (dur - 1)]);
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			// pass lab to computeStateArrayValue() instead of clab
			//this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,clab);
			stateArrayForCurDur[lab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,lab);

			// just for debugging
//			cout << "lab=" << lab << ", dur=" << dur << endl << "TransMatrix: plab=";
			//cout << "clab=" << clab << endl << "TransMatrix: plab=";

			//TODO: for (QNUInt32 plab=0; plab<nLabs_of_prevNode; plab++) {
			//for (QNUInt32 plab=0; plab<nLabs; plab++) {
			for (QNUInt32 plab = 0; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {

				//QNUInt32 idx = plab * this->nLabs + clab;
				QNUInt32 idx = plab * this->nActualLabs + lab;
				//TODO: design a feature map in which the transition matrix calculation can take different dimensions of plab (from prevNode) and clab (from current node).
				//this->transMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(seg_ftr_buf,lambda,plab,clab);
				transMatrixForCurDur[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(seg_ftr_buf,lambda,plab,lab);

				// just for debugging
//				cout << "[" << plab << "]=" << this->transMatrix[idx] << " ";
			}

			// just for debugging
//			cout << endl;
//			cout << "StateArray[" << "lab=" << lab << ", dur=" << dur << "]=" << stateArrayForCurDur[lab] << endl;

//			clab++;

		}
		seg_ftr_buf += this->nFtrsPerSeg;

		// just for debugging
//		cout << "seg_ftr_buf moved forward by nFtrsPerSeg(" << nFtrsPerSeg << ")" << endl;
	}
	// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
	for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
	{
		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			//this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,clab);
			stateArrayForCurDur[lab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,lab);

			// just for debugging
//			cout << "lab=" << lab << ", dur=" << dur << endl << "TransMatrix: plab=";
			//cout << "clab=" << clab << endl << "TransMatrix: plab=";
//			cout << endl;
//			cout << "StateArray[" << "lab=" << lab << ", dur=" << dur << "]=" << stateArrayForCurDur[lab] << endl;

//			clab++;
		}
		seg_ftr_buf += this->nFtrsPerSeg;

		// just for debugging
//		cout << "seg_ftr_buf moved forward by nFtrsPerSeg(" << nFtrsPerSeg << ")" << endl;
	}
	return result;
}


/*
 * CRF_StdSegStateNode_WithoutDurLab::computeAlpha
 *
 * Read alpha vectors of previous nodes directly from prevNode and store the result of the alpha vector in alphaArray.
 *
 * Stub function.
 * Should compute the alpha vector for the forward backward computation for this node.
 */
double CRF_StdSegStateNode_WithoutDurLab::computeAlpha()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab::computeAlpha(): " << endl;

	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	checkNumPrevNodes();

//	double* tempAlphaArray = new double[this->nodeLabMaxDur];

//	QNUInt32 clab = 0;
	for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
	{
//		double maxv = LOG0;
		for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
		{
			CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
			double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
			double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
			double* transMatrixForCurDur = &(this->transMatrix[this->nActualLabs * this->nActualLabs * (dur - 1)]);

			QNUInt32 logAddID = 0;

			//this->logAddAcc[0]=prev_adj_seg_alpha[0]+this->transMatrix[0+clab];
			this->logAddAcc[0]=prev_adj_seg_alpha[0]+transMatrixForCurDur[0 * this->nActualLabs + lab];
			logAddID++;
			double maxv=this->logAddAcc[0];

			// just for debugging
//			cout << "lab=" << lab << ", dur=" << dur << endl << "alphaArray_WithDur = logAdd([prev alpha + trans], " << prevAdjacentSeg->getNumAvailLabs() << "): ";
//			cout << "plab=[0]=" << prev_adj_seg_alpha[0] <<
//					")+(" << transMatrixForCurDur[0 * this->nActualLabs + lab] <<
//					")=" << this->logAddAcc[0] << " ";

			//TODO: for (QNUInt32 plab = 1; plab < prevAdjacentSeg->getNLabs(); plab++) {  //full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
			for (QNUInt32 plab = 1; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {   //faster implementation, not guaranteed to work for all classes of previous nodes.
				//this->logAddAcc[plab]=prev_adj_seg_alpha[plab]+this->transMatrix[plab * prevAdjacentSeg->getNLabs() + clab];
				//this->logAddAcc[plab]=prev_adj_seg_alpha[plab]+this->transMatrix[plab * this->nLabs + clab];
				this->logAddAcc[logAddID]=prev_adj_seg_alpha[plab]+transMatrixForCurDur[plab * this->nActualLabs + lab];

				// just for debugging
//				cout << "plab=[" << plab << "]=(" << prev_adj_seg_alpha[plab] <<
//						")+(" << transMatrixForCurDur[plab * this->nActualLabs + lab] <<
//						")=" << this->logAddAcc[plab] << " ";

				if (this->logAddAcc[logAddID]>maxv) {
					maxv=logAddAcc[logAddID];
				}

				logAddID++;
			}

			try {
				//TODO: this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,prevAdjacentSeg->getNLabs());        //full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
				//this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,prevAdjacentSeg->getNumAvailLabs());   //faster implementation, not guaranteed to work for all classes of previous nodes.
				this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] = logAdd(this->logAddAcc,maxv,logAddID);   //faster implementation, not guaranteed to work for all classes of previous nodes.
			}
			catch (exception &e) {
				string errstr="CRF_StdSegStateNode_WithoutDurLab::computeAlpha() caught exception: "+string(e.what())+" while computing alpha";
				throw runtime_error(errstr);
				return(-1);
			}
			// just for debugging
//			cout << endl << "alphaArray_WithDur[lab=" << lab << ",dur=" << dur << "](" << this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] << ") + stateArray[lab=" << lab << ",dur=" << dur << "](" << stateArrayForCurDur[lab] << ") =";

			//this->alphaArray[clab]+=this->stateArray[clab];
			this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] += stateArrayForCurDur[lab];

			// just for debugging
//			cout << this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] << endl;

//			clab++;
		}
		// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
		for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
		{
			double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
			//this->alphaArray[clab]=this->stateArray[clab];
			this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] = stateArrayForCurDur[lab];

			// just for debugging
//			cout << "lab =" << lab << ", dur=" << dur << endl;
//			cout << "alphaArray_WithDur[lab=" << lab << ",dur=" << dur << "]=stateArray[lab=" << lab << ",dur=" << dur << "]=" << this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] << endl;

//			clab++;
		}

		try {
			this->alphaArray[lab] = logAdd(&(this->alphaArray_WithDur[lab*this->nodeLabMaxDur]),this->nodeLabMaxDur);
		}
		catch (exception &e) {
			string errstr="CRF_StdSegStateNode_WithoutDurLab::computeAlpha() caught exception: "+string(e.what())+" while computing alpha";
			throw runtime_error(errstr);
			return(-1);
		}

		// just for debugging
//		cout << "lab=" << lab << endl;
//		cout << "alphaArray[" << lab << "]=logAdd(alphaArray_WithDur[" << lab*this->nodeLabMaxDur << "], " << this->nodeLabMaxDur << ")=" << this->alphaArray[lab] << endl;
//		cout << "logAdd(" << this->alphaArray_WithDur[lab*this->nodeLabMaxDur] << ", " << this->alphaArray_WithDur[lab*this->nodeLabMaxDur + 1] << "), nodeLabMaxDur=" << this->nodeLabMaxDur << endl;
	}

//	delete [] tempAlphaArray;

	// just for debugging
//	cout << endl;

	return this->alphaScale;

}

/*
 * CRF_StdSegStateNode_WithoutDurLab::computeFirstAlpha
 *
 * Stub function.
 * Should compute the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdSegStateNode_WithoutDurLab::computeFirstAlpha()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab::computeFirstAlpha()" << endl;

	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	//nodeMaxLab for the first node of the sequence is usually equal to nActualLabs (since nodeLabMaxDur==1).
	for (QNUInt32 lab = 0; lab < this->numAvailLabs; lab++)
	{
		// just for debugging
//		cout << "alphaArray_WithDur[lab=" << lab << ",dur=" << 1 << "]=stateArray[" << lab << "]=" << this->stateArray[lab] << endl;

		//this->alphaArray[clab]+=this->stateArray[clab];
		this->alphaArray_WithDur[lab] = this->stateArray[lab];

		// just for debugging
//		cout << "alphaArray[" << lab << "]=stateArray[" << lab << "]=" << this->stateArray[lab] << endl;

		this->alphaArray[lab]=this->stateArray[lab];
	}

	// just for debugging
//	cout << endl;

	return this->alphaScale;

}


/*
 * CRF_StdSegStateNode_WithoutDurLab::computeBeta
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
double CRF_StdSegStateNode_WithoutDurLab::computeBeta(double scale)
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
//	cout << "CRF_StdSegStateNode_WithoutDurLab::computeBeta(): " << endl;

	checkNumNextNodes();

	// if numNextNodes == 0, this is the last node of the sequence.
	// Sets the beta value in this node to the special case for the end of the sequence.
	if (this->numNextNodes == 0)
	{
		setTailBeta();
		return this->alphaScale;
	}

//	QNUInt32 nextlab = 0;
//	for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
//	{
//		CRF_StateNode* nextAdjacentSeg = this->nextNodes[dur - 1];
//		double* next_adj_seg_beta = nextAdjacentSeg->getBeta();
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++) {
//			this->tempBeta[nextlab] = next_adj_seg_beta[nextlab] + nextAdjacentSeg->getStateValue(nextlab);
//
//			// just for debugging
////			cout << "tempBeta[" << nextlab << "]=" << "next_seg_beta[" << nextlab << "](" << next_adj_seg_beta[nextlab] << ") + next_seg_state_value[" << nextlab << "](" << nextAdjacentSeg->getStateValue(nextlab) << ")=" << this->tempBeta[nextlab] << endl;
//
//			nextlab++;
//		}
//	}

//	QNUInt32 nextlab = 0;
	for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
	{
		CRF_StateNode* nextAdjacentSeg = this->nextNodes[dur - 1];
		double* next_adj_seg_beta = nextAdjacentSeg->getBeta();

//		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
		double* tempBetaForCurDur = &(this->tempBeta[this->nActualLabs * (dur - 1)]);
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			//this->tempBeta[nextlab] = next_adj_seg_beta[nextlab] + nextAdjacentSeg->getStateValue(nextlab);
			tempBetaForCurDur[lab] = next_adj_seg_beta[lab] + nextAdjacentSeg->getStateValue(lab,dur);
//			this->tempBeta[this->nActualLabs * (dur - 1) + lab] = next_adj_seg_beta[lab] + nextAdjacentSeg->getStateValue(lab,dur);

			// just for debugging
//			cout << "tempBeta[lab=" << lab << ",dur=" << dur << "]=" << "next_adj_seg_beta[" << lab << "](" << next_adj_seg_beta[lab] << ") + next_seg_state_value[lab=" << lab << ",dur=" << dur << "](" << nextAdjacentSeg->getStateValue(lab,dur) << ")=" << tempBetaForCurDur[lab] << endl;

			//nextlab++;
		}
	}


	for (QNUInt32 clab = 0; clab < this->numAvailLabs; clab++)
	{
		//just for debugging
//		cout << "clab=" << clab << endl;

		CRF_StateNode* nextAdjacentSeg = this->nextNodes[0];
		//this->logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, 0) + this->tempBeta[0];
//		this->logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, 0, 1) + nextAdjacentSeg->getTempBeta(0, 1);
		this->logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, 0, 1) + this->tempBeta[0];
		double maxv=this->logAddAcc[0];
//		QNUInt32 nextlab = 0;
		QNUInt32 logAddID = 0;

		// just for debugging
//		cout << "clab=" << clab << ", betaArray = logAdd(trans + tempBeta):" << endl;

		for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
		{
			nextAdjacentSeg = this->nextNodes[dur - 1];
			double* tempBetaForCurDur = &(this->tempBeta[this->nActualLabs * (dur - 1)]);
			//TODO: It should be nActualLabs_of_nextNode instead of nActualLabs_of_thisNode as the number of iterations.
			for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++) {
//			for (QNUInt32 lab = 0; lab < this->nextNodeNActualLabs; lab++) {

				//this->logAddAcc[nextlab] = nextAdjacentSeg->getTransValue(clab, nextlab) + this->tempBeta[nextlab];
//				this->logAddAcc[logAddID] = nextAdjacentSeg->getTransValue(clab, lab, dur) + nextAdjacentSeg->getTempBeta(lab, dur);
				this->logAddAcc[logAddID] = nextAdjacentSeg->getTransValue(clab, lab, dur) + tempBetaForCurDur[lab];
//				this->logAddAcc[logAddID] = nextAdjacentSeg->getTransValue(clab, lab, dur) + this->tempBeta[this->nActualLabs * (dur - 1) + lab];

				// just for debugging
//				cout << " next lab[" << lab << "]_dur[" << dur << "]=" << "trans(" << nextAdjacentSeg->getTransValue(clab, lab, dur) << ") + tempBeta(" << nextAdjacentSeg->getTempBeta(lab, dur) << ")="<< this->logAddAcc[logAddID] << endl;
//				cout << " next lab[" << lab << "]_dur[" << dur << "]=" << "trans(" << nextAdjacentSeg->getTransValue(clab, lab, dur) << ") + tempBeta(" << tempBetaForCurDur[lab] << ")="<< this->logAddAcc[logAddID] << endl;

//				if (this->logAddAcc[nextlab]>maxv) {
//					maxv=logAddAcc[nextlab];
//				}
				if (this->logAddAcc[logAddID]>maxv) {
					maxv=logAddAcc[logAddID];
				}
//				nextlab++;
				logAddID++;
			}
		}
		try {
			//TODO:It should be nActualLabs_of_nextNode instead of nActualLabs_of_thisNode.
			//this->betaArray[clab]=logAdd(this->logAddAcc,maxv,this->nActualLabs * this->numNextNodes);
			//this->betaArray[clab]=logAdd(this->logAddAcc,maxv,this->nextNodeNActualLabs * this->numNextNodes);
			this->betaArray[clab]=logAdd(this->logAddAcc,maxv,logAddID);

			// just for debugging
//			cout << "betaArray[" << clab << "] = logAdd(trans + tempBeta) = " << this->betaArray[clab] << endl;
		}
		catch (exception &e) {
			string errstr="CRF_StdSegStateNode_WithoutDurLab::computeBeta() caught exception: "+string(e.what())+" while computing beta";
			throw runtime_error(errstr);
			return(-1);
		}
	}

	// calculate tempBeta on the current node for later use in the beta calculation of previous nodes.
//	for (QNUInt32 dur = 1; dur <= this->nodeLabMaxDur; dur++)
//	{
//		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
//		double* tempBetaForCurDur = &(this->tempBeta[this->nActualLabs * (dur - 1)]);
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			tempBetaForCurDur[lab] = this->betaArray[lab] + stateArrayForCurDur[lab];
//		}
//	}

	// just for debugging
//	cout << endl;

	return this->alphaScale;
}

///*
// * CRF_StdSegStateNode_WithoutDurLab::setTailBeta
// *
// * Sets the beta value and the *tempbeta* value in this node to the special case for the end of the sequence.
// */
//void CRF_StdSegStateNode_WithoutDurLab::setTailBeta()
//{
//	// just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab::setTailBeta()." << endl;
//
//	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
//	for (QNUInt32 clab=0; clab<nLabs; clab++) {
//		this->betaArray[clab]=0.0;
//
//		// just for debugging
//		cout << "this->betaArray[clab=" << clab << "]=0.0;" << endl;
//	}
//
//	for (QNUInt32 clab = 0; clab < this->nActualLabs; clab++)
//	{
//		// this->numNextNodes would be 0 for the last node of the sequence
//		//for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
//		for (QNUInt32 dur = 1; dur <= this->nodeLabMaxDur; dur++)
//		{
//			this->tempBeta[this->nActualLabs * (dur - 1) + clab] = 0.0;
//
//			// just for debugging
//			cout << "this->tempBeta[clab=" << clab << ", dur=" << dur << "]=0.0;" << endl;
//		}
//	}
//
//	// just for debugging
//	cout << "Tail betas have been set." << endl;
//
//	// just for debugging
//	cout << endl;
//}

/*
 * CRF_StdSegStateNode_WithoutDurLab::computeExpF
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
double CRF_StdSegStateNode_WithoutDurLab::computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 prev_lab)
{
	// Added by Ryan, just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab::computeExpF(): Very beginning!!" << endl;

	checkNumPrevNodes();

	// Added by Ryan, just for debugging
//	cout << "checkNumPrevNodes() passed." << endl;

	QNUInt32 actualLab = this->label;
	QNUInt32 labDur = CRF_LAB_BAD;
	QNUInt32 actualPLab = prev_lab;
	QNUInt32 plabDur = CRF_LAB_BAD;

	// Added by Ryan, just for debugging
//	cout << "Before label conversion:" << endl;
//	cout << "actualLab=" << actualLab << endl;
//	cout << "actualPLab=" << actualPLab << endl;

	if (actualLab != CRF_LAB_BAD)
	{
		if (actualLab >= this->nActualLabs * this->labMaxDur)
		{
			string errstr="CRF_StdSegStateNode_WithoutDurLab::computeExpF() threw exception: the label is larger than "+ stringify(this->nActualLabs * this->labMaxDur);
			throw runtime_error(errstr);
		}
		actualLab = this->label % this->nActualLabs;
		labDur = this->label / this->nActualLabs + 1;
	}
	if (actualPLab != CRF_LAB_BAD)
	{
		if (actualPLab >= this->nActualLabs * this->labMaxDur)
		{
			string errstr="CRF_StdSegStateNode_WithoutDurLab::computeExpF() threw exception: the previous label is larger than "+ stringify(this->nActualLabs * this->labMaxDur);
			throw runtime_error(errstr);
		}
		actualPLab = prev_lab % this->nActualLabs;
		plabDur = prev_lab / this->nActualLabs + 1;
	}

	// Added by Ryan, just for debugging
//	cout << "After label conversion:" << endl;
//	cout << "actualLab=" << actualLab << endl;
//	cout << "actualPLab=" << actualPLab << endl;

	double logLi=0.0;
	double alpha_beta=0.0;
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;
	float* seg_ftr_buf = this->ftrBuf;
//	QNUInt32 clab = 0;
	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
	{
		// just for debugging
//		cout << "dur: " << dur << endl;

		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
		double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
		double* transMatrixForCurDur = &(this->transMatrix[this->nActualLabs * this->nActualLabs * (dur - 1)]);

		// just for debugging
//		cout << "After getting prevAdjacentSeg." << endl;

		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			// just for debugging
//			cout << "First phase:: numPrevNodes: " << numPrevNodes << ", dur: " << dur << ", lab: " << lab << endl;

			//alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
			alpha_beta=expE(this->alphaArray_WithDur[lab * this->nodeLabMaxDur + dur - 1] + this->betaArray[lab] - Zx);
			alpha_beta_tot += alpha_beta;
			//bool match=(clab==this->label);
			bool match=(lab==actualLab && dur==labDur);

			// just for debugging
//			cout << "alpha_beta=" << alpha_beta << ", alpha_beta_tot=" << alpha_beta_tot << endl;

			// just for debugging
//			if (match)
				//cout << clab << " ";
//				cout << "state label match: actualLab=" << actualLab << ", duration=" << labDur << endl;

			//logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,this->label,clab);
			if (match)
			{
				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,actualLab,lab);
			} else {
				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,lab);
			}
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

					// just for debugging
//					cout << "plab: " << plab << ", ";

					//QNUInt32 idx = plab * this->nLabs + clab;
					QNUInt32 idx = plab * this->nActualLabs + lab;
					//alpha_beta=expE(prev_adj_seg_alpha[plab]+this->transMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
					alpha_beta=expE(prev_adj_seg_alpha[plab]+transMatrixForCurDur[idx]+stateArrayForCurDur[lab]+this->betaArray[lab]-Zx);
					alpha_beta_trans_tot+=alpha_beta;
					//match=((clab==this->label)&&(plab==prev_lab));
					match=(lab==actualLab && dur==labDur && plab==actualPLab);

					// just for debugging
//					cout << "alpha_beta=" << alpha_beta << ", alpha_beta_trans_tot=" << alpha_beta_trans_tot << endl;

					// just for debugging
//					if (match)
						//cout << "(" << clab << "<-" <<  plab << ") ";
//						cout << "transition label match: " << "(" << actualLab << "<-" <<  actualPLab << ") " << endl;

					// TODO: design a feature map in which the transition matrix calculation can take different size of plab (from prevNode) and clab (from current node).
					//logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab);
					if (match)
					{
						logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,actualPLab,actualLab,plab,lab);
					} else {
						logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,CRF_LAB_BAD,plab,lab);
					}
				}
//			}
//			clab++;
		}
		seg_ftr_buf += this->nFtrsPerSeg;
	}
	// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
	for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
	{
		// just for debugging
//		cout << "dur: " << dur << endl;

		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			// just debugging
//			cout << "Second phase:: numPrevNodes: " << numPrevNodes << ", dur: " << dur << ", lab: " << lab << endl;

			//alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
			alpha_beta=expE(this->alphaArray_WithDur[lab * this->nodeLabMaxDur + dur - 1] + this->betaArray[lab] - Zx);
			alpha_beta_tot += alpha_beta;
			//bool match=(clab==this->label);
			bool match=(lab==actualLab && dur==labDur);

			// just for debugging
//			cout << "alpha_beta=" << alpha_beta << ", alpha_beta_tot=" << alpha_beta_tot << endl;

			// just for debugging
//			if (match)
				//cout << clab << " ";
//				cout << "state label match: actualLab=" << actualLab << ", duration=" << labDur << endl;


			//logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,this->label,clab);
			if (match)
			{
				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,actualLab,lab);
			} else {
				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,lab);
			}

//			clab++;
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
		string errstr="CRF_StdSegStateNode_WithoutDurLab::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
//	else if (alpha_beta_tot < 0.999999) {
	else if (alpha_beta_tot < -0.000001) {
//		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		string errstr="CRF_StdSegStateNode_WithoutDurLab::computeExpF() threw exception: Probability sums less than 0.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.000001) {
		string errstr="CRF_StdSegStateNode_WithoutDurLab::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
//	else if (alpha_beta_trans_tot < 0.999999) {
	else if (alpha_beta_trans_tot < -0.000001) {
//		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		string errstr="CRF_StdSegStateNode_WithoutDurLab::computeExpF() threw exception: Trans Probability sums less than 0.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}

	// just for debugging
//	cout << endl;

	return logLi;
}

// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur)
 *
 * Correct getTransValue() function for CRF_StdSegStateNode_WithoutDurLab.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur)
{
	// just for debugging
//	cout << "getting transMatrix[" << this->nActualLabs << "*" << this->nActualLabs << "*" << (dur - 1) << "+" << prev_lab
//			<< "*" << this->nActualLabs << "+" << cur_lab << "]=transMatrix["
//			<< this->nActualLabs * this->nActualLabs * (dur - 1) + prev_lab * this->nActualLabs + cur_lab << "]="
//			//<< this->transMatrix[this->nActualLabs * this->nActualLabs * (dur - 1) + prev_lab * this->nActualLabs + cur_lab]
//			<< endl;

	return this->transMatrix[this->nActualLabs * this->nActualLabs * (dur - 1) + prev_lab * this->nActualLabs + cur_lab];
}

// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab::getStateValue(QNUInt32 cur_lab, QNUInt32 dur)
 *
 * Correct getStateValue() function for CRF_StdSegStateNode_WithoutDurLab.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab::getStateValue(QNUInt32 cur_lab, QNUInt32 dur)
{
	return this->stateArray[this->nActualLabs * (dur - 1) + cur_lab];
}

// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur)
 *
 * Correct getFullTransValue() function for CRF_StdSegStateNode_WithoutDurLab.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur)
{
	return getTransValue(prev_lab, cur_lab, dur) + getStateValue(cur_lab, dur);
}

// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab::getTempBeta(QNUInt32 cur_lab, QNUInt32 dur)
 *
 */
double CRF_StdSegStateNode_WithoutDurLab::getTempBeta(QNUInt32 cur_lab, QNUInt32 dur)
{
	return this->tempBeta[this->nActualLabs * (dur - 1) + cur_lab];
}

// Disable all these functions by overriding them with exception handling. Use their modified version below.
// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
 *
 * Disabled in CRF_StdSegStateNode_WithoutDurLab, use getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur) instead.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	string errstr="Error: use CRF_StdSegStateNode_WithoutDurLab::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur) instead.";
	throw runtime_error(errstr);
}

// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab::getStateValue(QNUInt32 cur_lab)
 *
 * Disabled in CRF_StdSegStateNode_WithoutDurLab, use getStateValue(QNUInt32 cur_lab, QNUInt32 dur) instead.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab::getStateValue(QNUInt32 cur_lab)
{
	string errstr="Error: use CRF_StdSegStateNode_WithoutDurLab::getStateValue(QNUInt32 cur_lab, QNUInt32 dur) instead.";
	throw runtime_error(errstr);
}

// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
 *
 * Disabled in CRF_StdSegStateNode_WithoutDurLab, use getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur) instead.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	string errstr="Error: use CRF_StdSegStateNode_WithoutDurLab::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur) instead.";
	throw runtime_error(errstr);
}

