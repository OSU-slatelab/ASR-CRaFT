/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: hey
 */

#include "CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr.h"
#include <omp.h>

CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs)
	: CRF_StdSegStateNode_WithoutDurLab(fb, sizeof_fb, lab, crf, nodeMaxDur, prevNode_nLabs, nextNode_nActualLabs)
{
	delete [] this->transMatrix;
	this->transMatrix = new double[this->nActualLabs * this->nActualLabs];
//	for (QNUInt32 id = 0; id < this->nActualLabs * this->nActualLabs; id++)
//	{
//      // Changed by Ryan. It should be CRF_LogMath::LOG0 instead of 0.0. TODO: verify.
//		//this->transMatrix[id] = 0.0;
//		this->transMatrix[id] = CRF_LogMath::LOG0;
//	}

	alphaPlusTrans = new double[this->nActualLabs];
	tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur = new double[this->nActualLabs];
}

CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::~CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr()
{
	delete [] this->alphaPlusTrans;
	delete [] this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur;
}

/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeTransMatrix
 *
 * Computes the log of the state vector (stateArray) and the log of the transition matrix (transMatrix)
 *   and stores them as appropriate
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeTransMatrix()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeTransMatrix(): " << endl;

	checkNumPrevNodes();

	double result=0.0;

	double* lambda = this->crf_ptr->getLambda();

	// fill in the transition matrix
	// Important Note: also need to fill in the transition matrix even for the first node,
	//                 since it will be used in computeExpF().
//	if (this->numPrevNodes > 0)
//	{
//		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - 1];
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			// just for debugging
//			cout << "lab=" << lab << endl << "TransMatrix: plab=";
			//cout << "clab=" << clab << endl << "TransMatrix: plab=";

			// assuming the number of labels in the previous node is equal to that in the current node
			//for (QNUInt32 plab = 0; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {
			for (QNUInt32 plab = 0; plab < this->nActualLabs; plab++) {

				//TODO: design a feature map in which the transition matrix calculation can take different dimensions of plab (from prevNode) and clab (from current node).

				//QNUInt32 idx = plab * this->nLabs + clab;
				//this->transMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(seg_ftr_buf,lambda,plab,clab);

//				QNUInt32 idx = plab * this->nActualLabs + lab;
//				transMatrixForCurDur[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(seg_ftr_buf,lambda,plab,lab);

				// in this class, transMatrix of the current node only stores
				// the transition features for duration 1. So we pass this->ftrBuf
				// as the feature vector to computeTransMatrixValue() since the
				// first set of features in this->ftrBuf is for duration 1.
				QNUInt32 idx = plab * this->nActualLabs + lab;
				this->transMatrix[idx]=this->crf_ptr->getFeatureMap()->computeTransMatrixValue(this->ftrBuf,lambda,plab,lab);

				// just for debugging
//				cout << "[" << plab << "]=" << this->transMatrix[idx] << " ";
			}

			// just for debugging
//			cout << endl;

//			clab++;
		}
//	}

	// fill in the state array
	float* seg_ftr_buf = this->ftrBuf;
//	QNUInt32 clab = 0;
	// Note: it should be this->numPrevNodes instead of this->nodeLabMaxDur as the number of iterations of the outer loop.
	// It's because numPrevNodes is enough for transition calculation while nodeLabMaxDur might be larger than numPrevNodes for some nodes.
	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
	{
//		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
//		double* transMatrixForCurDur = &(this->transMatrix[this->nActualLabs * this->nActualLabs * (dur - 1)]);
		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
		{
			// pass lab to computeStateArrayValue() instead of clab
			//this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,clab);
			stateArrayForCurDur[lab]=this->crf_ptr->getFeatureMap()->computeStateArrayValue(seg_ftr_buf,lambda,lab);

			// just for debugging
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
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha
 *
 * Read alpha vectors of previous nodes directly from prevNode and store the result of the alpha vector in alphaArray.
 *
 * Should compute the alpha vector for the forward backward computation for this node.
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha(): " << endl;

	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	checkNumPrevNodes();

//	double* tmpAlphaArray_gatherPreviousNodes = new double[this->nActualLabs];
//	double* tmpAlphaArray_gatherPreviousNodes_addTransValue = new double[this->nActualLabs];
	QNUInt32 logAddAccSize = this->nActualLabs;
	if (this->labMaxDur > logAddAccSize)
		logAddAccSize = this->labMaxDur;
	double* tmpLogAddAcc = new double[logAddAccSize];


	// add the transition values to the alphas in the previous node for future use

	// just for debugging
//	cout << "Add the transition values to the alphas in the previous node..." << endl;

	if (this->numPrevNodes >= 1)
	{
		CRF_StateNode* prevAdjacentSegForDurOne = this->prevNodes[this->numPrevNodes - 1];
		CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr* prevAdjacentSegForDurOne_dnmcast;
		try {
			prevAdjacentSegForDurOne_dnmcast = dynamic_cast<CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr*>(prevAdjacentSegForDurOne);
			if (prevAdjacentSegForDurOne_dnmcast == 0)
			{
				string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: the previous node is a base class object of CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr.";
				throw runtime_error(errstr);
			}
		}
		catch (exception &e) {
			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: "+string(e.what())+", while dynamic casting previous nodes.";
			throw runtime_error(errstr);
		}
		prevAdjacentSegForDurOne_dnmcast->computeAlphaPlusTrans(this->transMatrix);
	}


//	// gather (log-)alphaPlusTrans from previous nodes,
//	// (log-)summing (log-)alphaPlusTrans for the same previous label over different previous nodes
//
//	// just for debugging
//	cout << "gathering (log-)alphaPlusTrans from previous nodes..." << endl;

	// compute log-alphas

	// just for debugging
//	cout << "Computing log-alphas..." << endl;

	// should we use "plab < prevAdjacentSeg->getNumAvailLabs()" or "plab < this->nActualLabs"? -
	// since we assume the number of labels in the previous node is equal to that in the current node,
	// prevAdjacentSeg->getNumAvailLabs() == this->nActualLabs
	for (QNUInt32 clab = 0; clab < this->nActualLabs; clab++)
	{
		QNUInt32 logAddID = 0;
		for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
		{
			CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
//			double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
			CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr* prevAdjacentSeg_dnmcast;
			try {
				prevAdjacentSeg_dnmcast = dynamic_cast<CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr*>(prevAdjacentSeg);
				if (prevAdjacentSeg_dnmcast == 0)
				{
					string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: the previous node is a base class object of CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr.";
					throw runtime_error(errstr);
				}
			}
			catch (exception &e) {
				string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: "+string(e.what())+", while dynamic casting previous nodes.";
				throw runtime_error(errstr);
			}
//			tmpLogAddAcc[logAddID++] = prev_adj_seg_alpha[plab];
			tmpLogAddAcc[logAddID++] = prevAdjacentSeg_dnmcast->getAlphaPlusTrans(clab) +
					this->stateArray[this->nActualLabs * (dur - 1) + clab];
			this->alphaArray_WithDur[clab*this->nodeLabMaxDur + dur - 1] =
					prevAdjacentSeg_dnmcast->getAlphaPlusTrans(clab) +
					this->stateArray[this->nActualLabs * (dur - 1) + clab];

			// just for debugging
//			cout << "this->alphaArray_WithDur[clab=" << clab << ", dur=" << dur << "] = " <<
//					"prevAdjacentSeg_dnmcast(dur=" << dur << ")->getAlphaPlusTrans(clab=" << clab << ")(" << prevAdjacentSeg_dnmcast->getAlphaPlusTrans(clab) <<
//					") + stateArray[clab=" << clab << ", dur=" << dur << "](" << this->stateArray[this->nActualLabs * (dur - 1) + clab] <<
//					") = " << this->alphaArray_WithDur[clab*this->nodeLabMaxDur + dur - 1] << endl;
		}
		for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
		{
			tmpLogAddAcc[logAddID++] = this->stateArray[this->nActualLabs * (dur - 1) + clab];
			this->alphaArray_WithDur[clab*this->nodeLabMaxDur + dur - 1] =
					this->stateArray[this->nActualLabs * (dur - 1) + clab];

			// just for debugging
//			cout << "this->alphaArray_WithDur[clab=" << clab << ", dur=" << dur << "] = " <<
//					"stateArray[clab=" << clab << ", dur=" << dur << "](" << this->stateArray[this->nActualLabs * (dur - 1) + clab] <<
//					") = " << this->alphaArray_WithDur[clab*this->nodeLabMaxDur + dur - 1] << endl;
		}
		try {
//			tmpAlphaArray_gatherPreviousNodes[plab] = logAdd(tmpLogAddAcc,logAddID);
//
//			// just for debugging
//			cout << "tmpAlphaArray_gatherPreviousNodes[plab=" << plab << "]=" << tmpAlphaArray_gatherPreviousNodes[plab] << endl;

			this->alphaArray[clab] = logAdd(tmpLogAddAcc,logAddID);

			// just for debugging
//			cout << "alphaArray[clab=" << clab << "] = " << this->alphaArray[clab] << endl;
		}
		catch (exception &e) {
			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: "+string(e.what())+", while computing alpha";
			throw runtime_error(errstr);
			return(-1);
		}
	}

//	// add transition values (from current node to the next node) to the alphas for future use in next nodes.
//
//	// the transition matrix should have been the one in the next node, but since transition features are
//	// used but only the transition bias, the transition matrix in the current node would be the same as that
//	// in the next node.
//
//	// just for debugging
////	cout << "adding transition values to the computed log-alphas for future use in the next node..." << endl;
//
//	// assuming the number of labels in the next node is equal to that in the current node
//	for (QNUInt32 next_lab = 0; next_lab < this->nActualLabs; next_lab++)
//	{
//		QNUInt32 logAddID = 0;
//		for (QNUInt32 clab = 0; clab < this->nActualLabs; clab++)
//		{
//			tmpLogAddAcc[logAddID++] = this->alphaArray[clab] +
//					this->transMatrix[clab * this->nActualLabs + next_lab];
//
//			// just for debugging
////			cout << "alphaArray[clab=" << clab << "](" << this->alphaArray[clab] << ") + " <<
////					"transMatrix[clab=" << clab << ", next_lab=" << next_lab << "](" << this->transMatrix[clab * this->nActualLabs + next_lab] <<
////					") = " << tmpLogAddAcc[logAddID-1] << endl;
//		}
//		try {
//			this->alphaPlusTrans[next_lab] = logAdd(tmpLogAddAcc,logAddID);
//
//			// just for debugging
////			cout << "alphaPlusTrans[next_lab=" << next_lab << "]=" << this->alphaPlusTrans[next_lab] << endl;
//		}
//		catch (exception &e) {
//			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: "+string(e.what())+", while computing alpha";
//			throw runtime_error(errstr);
//			return(-1);
//		}
//	}


//	delete[] tmpAlphaArray_gatherPreviousNodes;
//	delete[] tmpAlphaArray_gatherPreviousNodes_addTransValue;
	delete[] tmpLogAddAcc;

	// just for debugging
//	cout << endl;

	return this->alphaScale;


//////////////////////////////////

//	// just for debugging
////	cout << "CRF_StdSegStateNode_WithoutDurLab::computeAlpha(): " << endl;
//
//	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
//	this->alphaScale=0.0;
//
//	checkNumPrevNodes();
//
////	double* tempAlphaArray = new double[this->nodeLabMaxDur];
//
////	QNUInt32 clab = 0;
//	for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//	{
////		double maxv = LOG0;
//		for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
//		{
//			CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
//			double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
//			double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
////			double* transMatrixForCurDur = &(this->transMatrix[this->nActualLabs * this->nActualLabs * (dur - 1)]);
//
//			QNUInt32 logAddID = 0;
//
//			//this->logAddAcc[0]=prev_adj_seg_alpha[0]+this->transMatrix[0+clab];
////			this->logAddAcc[0]=prev_adj_seg_alpha[0]+transMatrixForCurDur[0+lab];
//			this->logAddAcc[0]=prev_adj_seg_alpha[0]+this->transMatrix[0+lab];
//			logAddID++;
//			double maxv=this->logAddAcc[0];
//
//			// just for debugging
////			cout << "lab=" << lab << ", dur=" << dur << endl << "alphaArray_WithDur = logAdd([prev alpha + trans], " << this->nodeLabMaxDur << "): plab=[0]=" << this->logAddAcc[0] << " ";
//
//			//TODO: for (QNUInt32 plab = 1; plab < prevAdjacentSeg->getNLabs(); plab++) {  //full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
//			for (QNUInt32 plab = 1; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {   //faster implementation, not guaranteed to work for all classes of previous nodes.
//				//this->logAddAcc[plab]=prev_adj_seg_alpha[plab]+this->transMatrix[plab * prevAdjacentSeg->getNLabs() + clab];
//				//this->logAddAcc[plab]=prev_adj_seg_alpha[plab]+this->transMatrix[plab * this->nLabs + clab];
//				this->logAddAcc[logAddID]=prev_adj_seg_alpha[plab]+transMatrixForCurDur[plab * this->nActualLabs + lab];
//
//				// just for debugging
////				cout << "[" << plab << "]=" << this->logAddAcc[plab] << " ";
//
//				if (this->logAddAcc[logAddID]>maxv) {
//					maxv=logAddAcc[logAddID];
//				}
//
//				logAddID++;
//			}
//
//			try {
//				//TODO: this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,prevAdjacentSeg->getNLabs());        //full implementation. But not working now because logAdd() cannot do calculation on LOG0 yet.
//				//this->alphaArray[clab]=logAdd(this->logAddAcc,maxv,prevAdjacentSeg->getNumAvailLabs());   //faster implementation, not guaranteed to work for all classes of previous nodes.
//				this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] = logAdd(this->logAddAcc,maxv,logAddID);   //faster implementation, not guaranteed to work for all classes of previous nodes.
//			}
//			catch (exception &e) {
//				string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: "+string(e.what())+" while computing alpha";
//				throw runtime_error(errstr);
//				return(-1);
//			}
//			// just for debugging
////			cout << endl << "alphaArray_WithDur[lab=" << lab << ",dur=" << dur << "](" << this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] << ") + stateArray[lab=" << lab << ",dur=" << dur << "](" << stateArrayForCurDur[lab] << ") =";
//
//			//this->alphaArray[clab]+=this->stateArray[clab];
//			this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] += stateArrayForCurDur[lab];
//
//			// just for debugging
////			cout << this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] << endl;
//
////			clab++;
//		}
//		// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
//		for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
//		{
//			double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
//			//this->alphaArray[clab]=this->stateArray[clab];
//			this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] = stateArrayForCurDur[lab];
//
//			// just for debugging
////			cout << "lab =" << lab << ", dur=" << dur << endl;
////			cout << "alphaArray_WithDur[lab=" << lab << ",dur=" << dur << "]=stateArray[lab=" << lab << ",dur=" << dur << "]=" << this->alphaArray_WithDur[lab*this->nodeLabMaxDur + dur - 1] << endl;
//
////			clab++;
//		}
//
//		try {
//			this->alphaArray[lab] = logAdd(&(this->alphaArray_WithDur[lab*this->nodeLabMaxDur]),this->nodeLabMaxDur);
//		}
//		catch (exception &e) {
//			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: "+string(e.what())+" while computing alpha";
//			throw runtime_error(errstr);
//			return(-1);
//		}
//
//		// just for debugging
////		cout << "lab=" << lab << endl;
////		cout << "alphaArray[" << lab << "]=logAdd(alphaArray_WithDur[" << lab*this->nodeLabMaxDur << "], " << this->nodeLabMaxDur << ")=" << this->alphaArray[lab] << endl;
////		cout << "logAdd(" << this->alphaArray_WithDur[lab*this->nodeLabMaxDur] << ", " << this->alphaArray_WithDur[lab*this->nodeLabMaxDur + 1] << "), nodeLabMaxDur=" << this->nodeLabMaxDur << endl;
//	}
//
////	delete [] tempAlphaArray;
//
//	return this->alphaScale;

}

/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeFirstAlpha
 *
 * Should compute the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeFirstAlpha()
{
	// just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeFirstAlpha()" << endl;

	//QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;

	//nodeMaxLab for the first node of the sequence is usually equal to nActualLabs (since nodeLabMaxDur==1).
	for (QNUInt32 lab = 0; lab < this->numAvailLabs; lab++)
	{
		this->alphaArray_WithDur[lab] = this->stateArray[lab];

		// just for debugging
//		cout << "this->alphaArray_WithDur[lab=" << lab << ", dur=1] = " <<
//				"stateArray[lab=" << lab << ", dur=1](" << this->stateArray[lab] <<
//				") = " << this->alphaArray_WithDur[lab] << endl;

		this->alphaArray[lab]=this->stateArray[lab];

		// just for debugging
//		cout << "alphaArray[" << lab << "]=stateArray[" << lab << "]=" << this->alphaArray[lab] << endl;
	}


//	// add transition values (from current node to the next node) to the alphas for future use in next nodes.
//
//	// the transition matrix should have been the one in the next node, but since transition features are
//	// used but only the transition bias, the transition matrix in the current node would be the same as that
//	// in the next node.
//
//	// just for debugging
////	cout << "adding transition values to the computed log-alphas for future use in the next node..." << endl;
//
//	QNUInt32 logAddAccSize = this->nActualLabs;
//	if (this->labMaxDur > logAddAccSize)
//		logAddAccSize = this->labMaxDur;
//	double* tmpLogAddAcc = new double[logAddAccSize];
//
//	// assuming the number of labels in the next node is equal to that in the current node
//	for (QNUInt32 next_lab = 0; next_lab < this->nActualLabs; next_lab++)
//	{
//		QNUInt32 logAddID = 0;
//		for (QNUInt32 clab = 0; clab < this->nActualLabs; clab++)
//		{
//			tmpLogAddAcc[logAddID++] = this->alphaArray[clab] +
//					this->transMatrix[clab * this->nActualLabs + next_lab];
//
//			// just for debugging
////			cout << "alphaArray[clab=" << clab << "](" << this->alphaArray[clab] << ") + " <<
////					"transMatrix[clab=" << clab << ", next_lab=" << next_lab << "](" << this->transMatrix[clab * this->nActualLabs + next_lab] <<
////					") = " << tmpLogAddAcc[logAddID-1] << endl;
//		}
//		try {
//			this->alphaPlusTrans[next_lab] = logAdd(tmpLogAddAcc,logAddID);
//
//			// just for debugging
////			cout << "alphaPlusTrans[next_lab=" << next_lab << "]=" << this->alphaPlusTrans[next_lab] << endl;
//		}
//		catch (exception &e) {
//			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlpha() caught exception: "+string(e.what())+", while computing alpha";
//			throw runtime_error(errstr);
//			return(-1);
//		}
//	}
//
//	delete[] tmpLogAddAcc;

	// just for debugging
//	cout << endl;

	return this->alphaScale;

}

/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeBeta
 *
 * Inputs: scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Read the beta vectors of next nodes directly from nextNode and store the result of the beta vector in betaArray.
 *
 * Should compute the beta vector for the node before this one and store it in result_beta
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeBeta(double scale)
{
	// just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeBeta(): " << endl;

	checkNumNextNodes();

	// if numNextNodes == 0, this is the last node of the sequence.
	// Sets the beta value in this node to the special case for the end of the sequence.
	if (this->numNextNodes == 0)
	{
		setTailBeta();
		return this->alphaScale;
	}

	QNUInt32 logAddAccSize = this->nActualLabs;
	if (this->labMaxDur > logAddAccSize)
		logAddAccSize = this->labMaxDur;
	double* tmpLogAddAcc = new double[logAddAccSize];

	// add state value of next nodes to the log-betas of next nodes,
	// (log-)summing the state-value-added log-betas for the same next label over different next nodes (i.e. different next durations)

	// just for debugging
//	cout << "adding state value to the log-betas of next nodes..." << endl;

	// assuming the number of labels in the next node is equal to that in the current node
	for (QNUInt32 nextlab = 0; nextlab < this->nActualLabs; nextlab++)
	{
		QNUInt32 logAddID = 0;
		for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
		{
			CRF_StateNode* nextAdjacentSeg = this->nextNodes[dur - 1];
			double* next_adj_seg_beta = nextAdjacentSeg->getBeta();
			this->tempBeta[this->nActualLabs * (dur - 1) + nextlab] =
					nextAdjacentSeg->getStateValue(nextlab,dur) + next_adj_seg_beta[nextlab];
			tmpLogAddAcc[logAddID++] = this->tempBeta[this->nActualLabs * (dur - 1) + nextlab];

			// just for debugging
//			cout << "nextAdjacentSeg->getStateValue(nextlab=" << nextlab << ",dur=" << dur << ")(" << nextAdjacentSeg->getStateValue(nextlab,dur) <<
//					") + next_adj_seg_beta[nextlab=" << nextlab << "] = " << tmpLogAddAcc[logAddID-1] << endl;
		}
		try {
			this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[nextlab] = logAdd(tmpLogAddAcc,logAddID);

			// just for debugging
//			cout << "this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[nextlab="
//					<< nextlab << "] = " << this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[nextlab] << endl;
		}
		catch (exception &e) {
			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeBeta() caught exception: "+string(e.what())+", while computing beta";
			throw runtime_error(errstr);
			return(-1);
		}
	}

	// add transition values to the state-value-added log-betas

	// just for debugging
//	cout << "adding transition values to the state-value-added log-betas..." << endl;

	if (this->numNextNodes == 0)
	{
		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeBeta() caught exception: numNextNodes has to be > 0";
		throw runtime_error(errstr);
	}
	// transition matrice are the same for all the next nodes,
	// since for this class transition features for the next node are the same for all following labMaxDur nodes.
	CRF_StateNode* nextAdjacentSeg = this->nextNodes[0];
	for (QNUInt32 clab = 0; clab < this->nActualLabs; clab++)
	{
		QNUInt32 logAddID = 0;
		for (QNUInt32 nextlab = 0; nextlab < this->nActualLabs; nextlab++)
		{
			tmpLogAddAcc[logAddID++] = nextAdjacentSeg->getTransValue(clab, nextlab) +
					this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[nextlab];

			// just for debugging
//			cout << "nextAdjacentSeg->getTransValue(clab=" << clab << ", nextlab=" << nextlab <<
//					")(" << nextAdjacentSeg->getTransValue(clab, nextlab) <<
//					") + tmpBetaArray_gatherNextNodes_addStateValue[nextlab=" << nextlab <<
//					"](" << this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[nextlab] << ") = " <<
//					tmpLogAddAcc[logAddID-1] << endl;
		}
		try {
			this->betaArray[clab] = logAdd(tmpLogAddAcc,logAddID);

			// just for debugging
//			cout << "betaArray[clab=" << clab << "]=" << this->betaArray[clab] << endl;
		}
		catch (exception &e) {
			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeBeta() caught exception: "+string(e.what())+", while computing beta";
			throw runtime_error(errstr);
			return(-1);
		}
	}

	delete[] tmpLogAddAcc;

	// just for debugging
//	cout << endl;

	return this->alphaScale;

	//////////////////////////////////////////////////////////

//	// Logic desired:
//	//	* Compute beta_i[size of alpha[]+1] to be all 1s
//	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]
//
////	QNUInt32 nLabs = this->crf_ptr->getNLabs();
////
////	for (QNUInt32 clab=0; clab<nLabs; clab++) {
////		this->tempBeta[clab]=this->betaArray[clab]+this->stateArray[clab];
////	}
////
////	for (QNUInt32 plab=0; plab<nLabs; plab++) {
////		this->logAddAcc[0]=this->transMatrix[plab*nLabs+0]+this->tempBeta[0];
////		double maxv=this->logAddAcc[0];
////		for (QNUInt32 clab=1; clab<nLabs; clab++) {
////			this->logAddAcc[clab]=this->transMatrix[plab*nLabs+clab]+this->tempBeta[clab];
////			if (this->logAddAcc[clab]>maxv) {
////				maxv=this->logAddAcc[clab];
////			}
////		}
////		try {
////			result_beta[plab]=logAdd(this->logAddAcc,maxv,nLabs);
////		}
////		catch (exception &e) {
////			string errstr="CRF_StdSegStateNode::computeBeta() caught exception: "+string(e.what())+" while computing beta";
////			throw runtime_error(errstr);
////			return(-1);
////		}
////	}
//
//	// just for debugging
////	cout << "CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeBeta(): " << endl;
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
////	QNUInt32 nextlab = 0;
////	for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
////	{
////		CRF_StateNode* nextAdjacentSeg = this->nextNodes[dur - 1];
////		double* next_adj_seg_beta = nextAdjacentSeg->getBeta();
////		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++) {
////			this->tempBeta[nextlab] = next_adj_seg_beta[nextlab] + nextAdjacentSeg->getStateValue(nextlab);
////
////			// just for debugging
//////			cout << "tempBeta[" << nextlab << "]=" << "next_seg_beta[" << nextlab << "](" << next_adj_seg_beta[nextlab] << ") + next_seg_state_value[" << nextlab << "](" << nextAdjacentSeg->getStateValue(nextlab) << ")=" << this->tempBeta[nextlab] << endl;
////
////			nextlab++;
////		}
////	}
//
////	QNUInt32 nextlab = 0;
//	for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
//	{
//		CRF_StateNode* nextAdjacentSeg = this->nextNodes[dur - 1];
//		double* next_adj_seg_beta = nextAdjacentSeg->getBeta();
//
////		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
//		double* tempBetaForCurDur = &(this->tempBeta[this->nActualLabs * (dur - 1)]);
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			//this->tempBeta[nextlab] = next_adj_seg_beta[nextlab] + nextAdjacentSeg->getStateValue(nextlab);
//			tempBetaForCurDur[lab] = next_adj_seg_beta[lab] + nextAdjacentSeg->getStateValue(lab,dur);
////			this->tempBeta[this->nActualLabs * (dur - 1) + lab] = next_adj_seg_beta[lab] + nextAdjacentSeg->getStateValue(lab,dur);
//
//			// just for debugging
////			cout << "tempBeta[lab=" << lab << ",dur=" << dur << "]=" << "next_adj_seg_beta[" << lab << "](" << next_adj_seg_beta[lab] << ") + next_seg_state_value[lab=" << lab << ",dur=" << dur << "](" << nextAdjacentSeg->getStateValue(lab,dur) << ")=" << tempBetaForCurDur[lab] << endl;
//
//			//nextlab++;
//		}
//	}
//
//
//	for (QNUInt32 clab = 0; clab < this->numAvailLabs; clab++)
//	{
//		//just for debugging
////		cout << "clab=" << clab << endl;
//
//		CRF_StateNode* nextAdjacentSeg = this->nextNodes[0];
//		//this->logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, 0) + this->tempBeta[0];
////		this->logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, 0, 1) + nextAdjacentSeg->getTempBeta(0, 1);
//		this->logAddAcc[0] = nextAdjacentSeg->getTransValue(clab, 0, 1) + this->tempBeta[0];
//		double maxv=this->logAddAcc[0];
////		QNUInt32 nextlab = 0;
//		QNUInt32 logAddID = 0;
//
//		// just for debugging
////		cout << "clab=" << clab << ", betaArray = logAdd(trans + tempBeta):" << endl;
//
//		for (QNUInt32 dur = 1; dur <= this->numNextNodes; dur++)
//		{
//			nextAdjacentSeg = this->nextNodes[dur - 1];
//			double* tempBetaForCurDur = &(this->tempBeta[this->nActualLabs * (dur - 1)]);
//			//TODO: It should be nActualLabs_of_nextNode instead of nActualLabs_of_thisNode as the number of iterations.
//			for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++) {
////			for (QNUInt32 lab = 0; lab < this->nextNodeNActualLabs; lab++) {
//
//				//this->logAddAcc[nextlab] = nextAdjacentSeg->getTransValue(clab, nextlab) + this->tempBeta[nextlab];
////				this->logAddAcc[logAddID] = nextAdjacentSeg->getTransValue(clab, lab, dur) + nextAdjacentSeg->getTempBeta(lab, dur);
//				this->logAddAcc[logAddID] = nextAdjacentSeg->getTransValue(clab, lab, dur) + tempBetaForCurDur[lab];
////				this->logAddAcc[logAddID] = nextAdjacentSeg->getTransValue(clab, lab, dur) + this->tempBeta[this->nActualLabs * (dur - 1) + lab];
//
//				// just for debugging
////				cout << " next lab[" << lab << "]_dur[" << dur << "]=" << "trans(" << nextAdjacentSeg->getTransValue(clab, lab, dur) << ") + tempBeta(" << nextAdjacentSeg->getTempBeta(lab, dur) << ")="<< this->logAddAcc[logAddID] << endl;
////				cout << " next lab[" << lab << "]_dur[" << dur << "]=" << "trans(" << nextAdjacentSeg->getTransValue(clab, lab, dur) << ") + tempBeta(" << tempBetaForCurDur[lab] << ")="<< this->logAddAcc[logAddID] << endl;
//
////				if (this->logAddAcc[nextlab]>maxv) {
////					maxv=logAddAcc[nextlab];
////				}
//				if (this->logAddAcc[logAddID]>maxv) {
//					maxv=logAddAcc[logAddID];
//				}
////				nextlab++;
//				logAddID++;
//			}
//		}
//		try {
//			//TODO:It should be nActualLabs_of_nextNode instead of nActualLabs_of_thisNode.
//			//this->betaArray[clab]=logAdd(this->logAddAcc,maxv,this->nActualLabs * this->numNextNodes);
//			//this->betaArray[clab]=logAdd(this->logAddAcc,maxv,this->nextNodeNActualLabs * this->numNextNodes);
//			this->betaArray[clab]=logAdd(this->logAddAcc,maxv,logAddID);
//
//			// just for debugging
////			cout << "betaArray[" << clab << "] = logAdd(trans + tempBeta) = " << this->betaArray[clab] << endl;
//		}
//		catch (exception &e) {
//			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeBeta() caught exception: "+string(e.what())+" while computing beta";
//			throw runtime_error(errstr);
//			return(-1);
//		}
//	}
//
//	// calculate tempBeta on the current node for later use in the beta calculation of previous nodes.
////	for (QNUInt32 dur = 1; dur <= this->nodeLabMaxDur; dur++)
////	{
////		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
////		double* tempBetaForCurDur = &(this->tempBeta[this->nActualLabs * (dur - 1)]);
////		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
////		{
////			tempBetaForCurDur[lab] = this->betaArray[lab] + stateArrayForCurDur[lab];
////		}
////	}
//
//	return this->alphaScale;
}

/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF
 *
 * Inputs: *ExpF - vector to store expected values of feature functions
 *         *grad - vector to store computed gradient values
 *         Zx - normalization constant
 *         t_next_lab - next node true label (transition feature ExpF computation)
 *         [important note: this is different from other models which pass prev_lab instead of t_next_lab]
 *
 * Returns:
 *
 * Read alpha vectors of previous nodes directly from prevNode for use in transition feature ExpF computation.
 *
 * Should compute gradient and expected values for features in this node and store them in *grad and
 *   *ExpF vectors respectively.  State features and transition features are computed in the same function.
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 t_next_lab)
{
	// Added by Ryan, just for debugging
//	cout << "CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF(): Very beginning!!" << endl;

	checkNumPrevNodes();

	// Added by Ryan, just for debugging
//	cout << "checkNumPrevNodes() passed." << endl;

	QNUInt32 actualLab = this->label;
	QNUInt32 labDur = CRF_LAB_BAD;
	QNUInt32 actualNextLab = t_next_lab;
	QNUInt32 plabDur = CRF_LAB_BAD;

	// Added by Ryan, just for debugging
//	cout << "Before label conversion:" << endl;
//	cout << "actualLab=" << actualLab << endl;
//	cout << "actualNextLab=" << actualNextLab << endl;

	if (actualLab != CRF_LAB_BAD)
	{
		if (actualLab >= this->nActualLabs * this->labMaxDur)
		{
			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: the label is larger than "+ stringify(this->nActualLabs * this->labMaxDur);
			throw runtime_error(errstr);
		}
		actualLab = this->label % this->nActualLabs;
		labDur = this->label / this->nActualLabs + 1;
	}
	if (actualNextLab != CRF_LAB_BAD)
	{
		if (actualNextLab >= this->nActualLabs * this->labMaxDur)
		{
			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: the previous label is larger than "+ stringify(this->nActualLabs * this->labMaxDur);
			throw runtime_error(errstr);
		}
		actualNextLab = t_next_lab % this->nActualLabs;
		plabDur = t_next_lab / this->nActualLabs + 1;
	}

	// Added by Ryan, just for debugging
//	cout << "After label conversion:" << endl;
//	cout << "actualLab=" << actualLab << endl;
//	cout << "actualNextLab=" << actualNextLab << endl;

	double logLi=0.0;
	double alpha_beta=0.0;
	//QNUInt32 nLabs = this->crf_ptr->getNLabs();

	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;

	// compute expected state features

	// just for debugging
//	cout << "computing expected state features: " << endl;

	// Changed by Ryan, using OpenMP for parallel
	QNUInt32 lab, dur, next_lab, clab;
//	int tid, nthreads, chunk = 6;
//	nthreads = omp_get_num_threads();
//	//omp_set_num_threads(4);
//	//printf("%d\n", omp_get_num_threads());
///*** Spawn a parallel region explicitly scoping all variables ***/
//#pragma omp parallel shared(nthreads,chunk) private(tid,lab,dur,next_lab,clab)
//	{
//		//tid = omp_get_thread_num();
//		//if (tid == 0)
//		//{
//		//	nthreads = omp_get_num_threads();
//		//	printf("Starting to compute expected state features with %d threads\n", nthreads);
//		//}
//		//printf("Thread %d starting to compute expected state features...\n",tid);
//	#pragma omp for schedule (static, chunk)
	for (lab = 0; lab < this->nActualLabs; lab++)
	{
		//printf("Thread=%d did lab=%d\n", tid, lab);

		//float* seg_ftr_buf = this->ftrBuf;

		//#pragma omp for schedule (static, 2)
		for (dur = 1; dur <= this->nodeLabMaxDur; dur++)
		{
			float* seg_ftr_buf = this->ftrBuf + (dur - 1) * this->nFtrsPerSeg;

			// just for debugging
//			cout << "dur: " << dur << endl;

			// just for debugging
//			cout << "numPrevNodes: " << numPrevNodes << ", dur: " << dur << ", lab: " << lab << endl;

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

			//seg_ftr_buf += this->nFtrsPerSeg;
		}
	}
//		// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
//		for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
//		{
//			// just for debugging
//	//		cout << "dur: " << dur << endl;
//
//			// just debugging
////			cout << "Second phase:: numPrevNodes: " << numPrevNodes << ", dur: " << dur << ", lab: " << lab << endl;
//
//			//alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
//			alpha_beta=expE(this->alphaArray_WithDur[lab * this->nodeLabMaxDur + dur - 1] + this->betaArray[lab] - Zx);
//			alpha_beta_tot += alpha_beta;
//			//bool match=(clab==this->label);
//			bool match=(lab==actualLab && dur==labDur);
//
//			// just for debugging
////			cout << "alpha_beta=" << alpha_beta << ", alpha_beta_tot=" << alpha_beta_tot << endl;
//
//			// just for debugging
////			if (match)
//				//cout << clab << " ";
////				cout << "state label match: actualLab=" << actualLab << ", duration=" << labDur << endl;
//
//
//			//logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,this->label,clab);
//			if (match)
//			{
//				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,actualLab,lab);
//			} else {
//				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,lab);
//			}
//
////			clab++;
//			seg_ftr_buf += this->nFtrsPerSeg;
//		}
	//}  /*** End of parallel region ***/

	// compute expected transition features

	// The transition calculated here is DIFFERENT from other model:
	// 		instead of transition from previous nodes to the current node,
	// 		the transition in this model is from the current node to next nodes.
	// This way of encoding transition makes possible the integration of
	// 		preceding alphas which end at t and succeeding betas which start
	//		at t (where t is the frame time stamp of the current node) such that
	// 		the transition features can be efficiently calculated in N*N time,
	// 		rather than N*N*D.
	// The overall expected features remain the same as the traditional way eventually.

	// just for debugging
//	cout << "computing expected transition features: " << endl;

	// just for debugging
//	cout << "numPrevNodes: " << numPrevNodes << endl;
//	cout << "numNextNodes: " << numNextNodes << endl;

	if (this->numNextNodes > 0)
	{
		CRF_StateNode* nextAdjacentSeg = this->nextNodes[0];

		// Changed by Ryan, using OpenMP for parallel
//		//QNUInt32 next_lab, clab;
//		//int tid, nthreads, chunk = 6;
//		//omp_set_num_threads(4);
//		//printf("%d\n", omp_get_num_threads());
///*** Spawn a parallel region explicitly scoping all variables ***/
////#pragma omp parallel shared(nthreads,chunk) private(tid,next_lab,clab)
//		//{
//			//tid = omp_get_thread_num();
//			//if (tid == 0)
//			//{
//			//	nthreads = omp_get_num_threads();
//			//	printf("Starting to compute expected transition features with %d threads\n", nthreads);
//			//}
//			//printf("Thread %d starting to compute expected transition features...\n",tid);
//		#pragma omp for schedule (static, chunk)
		// assuming the number of labels in the next node is equal to that in the current node
		for (next_lab = 0; next_lab < this->nActualLabs; next_lab++)
		{
			//printf("Thread=%d did next_lab=%d\n", tid, next_lab);

			// just for debugging
//			cout << "next_lab: " << next_lab << endl;

			//#pragma omp for schedule (static, chunk)
			for (clab = 0; clab < this->nActualLabs; clab++)
			{
				// just for debugging
//				cout << "clab: " << clab << ", ";

				// Since no transition features are used but only the transition bias, the transition
				// matrix of the current node can be used as the exact proxy of that of the next node.

				//QNUInt32 idx = plab * this->nLabs + clab;
//				QNUInt32 idx = clab * this->nActualLabs + next_lab;
				//alpha_beta=expE(prev_adj_seg_alpha[plab]+this->transMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
//				alpha_beta=expE(this->alphaArray[clab] + this->transMatrix[idx] + this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[next_lab] - Zx);

				alpha_beta=expE(this->alphaArray[clab] + nextAdjacentSeg->getTransValue(clab, next_lab)
						+ this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[next_lab] - Zx);

				alpha_beta_trans_tot+=alpha_beta;
				//match=((clab==this->label)&&(plab==prev_lab));
				bool match=(clab==actualLab && next_lab==actualNextLab);

				// just for debugging
//				cout << "alphaArray[clab=" << clab << "](" << this->alphaArray[clab] <<
//						") + nextAdjacentSeg->getTransValue[clab=" << clab << ", next_lab=" << next_lab << "](" << nextAdjacentSeg->getTransValue(clab, next_lab) <<
//						") + tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[next_lab=" << next_lab << "](" <<
//						this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[next_lab] <<
//						") - Zx(" << Zx <<
//						") = " << this->alphaArray[clab] + nextAdjacentSeg->getTransValue(clab, next_lab) + this->tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur[next_lab] - Zx << endl;
//				cout << "alpha_beta=" << alpha_beta << ", alpha_beta_trans_tot=" << alpha_beta_trans_tot << endl;

				// just for debugging
//				if (match)
					//cout << "(" << clab << "<-" <<  plab << ") ";
//					cout << "transition label match: " << "(" << actualNextLab << "<-" <<  actualLab << ") " << endl;

				// TODO: design a feature map in which the transition matrix calculation can take different size of plab (from prevNode) and clab (from current node).
				//logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab);

				// in this class, computeTransExpF() only needs the transition features in
				// the node for duration 1. So we pass nextAdjacentSeg->getFtrBuffer()
				// as the feature vector to computeTransExpF() since the first set of
				// features in nextAdjacentSeg->getFtrBuffer() is for duration 1.
				if (match)
				{
					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(nextAdjacentSeg->getFtrBuffer(),lambda,ExpF,grad,alpha_beta,actualLab,actualNextLab,clab,next_lab);
				} else {
					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(nextAdjacentSeg->getFtrBuffer(),lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,CRF_LAB_BAD,clab,next_lab);
				}
			}
		}


//	if (this->numPrevNodes > 0)
//	{
//		double* prevAlphaPlusState = new double[this->nActualLabs];
//		double* tmpLogAddAcc = new double[this->nodeLabMaxDur];
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			// just for debugging
//			cout << "lab: " << lab << endl;
//
//			QNUInt32 logAddID = 0;
//			for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
//			{
//				CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
//				double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
//
//				// just for debugging
//				cout << "After getting prevAdjacentSeg." << endl;
//
//				tmpLogAddAcc[logAddID++] = prev_adj_seg_alpha[lab] +
//						this->stateArray[this->nActualLabs * (dur - 1) + lab];
//
//				// just for debugging
//				cout << "prev_adj_seg_alpha[lab=" << lab << ", dur=" << dur << "]("
//						<< prev_adj_seg_alpha[lab] <<
//						") + stateArray[lab=" << lab << ", dur=" << dur << "]("
//						<< this->stateArray[this->nActualLabs * (dur - 1) + lab] <<
//						") = " << tmpLogAddAcc[logAddID-1] << endl;
//			}
//			for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
//			{
//				tmpLogAddAcc[logAddID++] = this->stateArray[this->nActualLabs * (dur - 1) + lab];
//
//				// just for debugging
//				cout << "stateArray[lab=" << lab << ", dur=" << dur << "]("
//						<< this->stateArray[this->nActualLabs * (dur - 1) + lab] <<
//						") = " << tmpLogAddAcc[logAddID-1] << endl;
//			}
//			try {
//				prevAlphaPlusState[lab] = logAdd(tmpLogAddAcc,logAddID);
//
//				// just for debugging
//				cout << "prevAlphaPlusState[lab=" << lab << "] = " << prevAlphaPlusState[lab] << endl;
//			}
//			catch (exception &e) {
//				string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() caught exception: "+string(e.what())+", while computing expected features";
//				throw runtime_error(errstr);
//				return(-1);
//			}
//		}
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			// just for debugging
//			cout << "lab: " << lab << endl;
//
//			// assuming the number of labels in the previous node is equal to that in the current node
////			for (QNUInt32 plab = 0; plab < prevAdjacentSeg->getNumAvailLabs(); plab++)
//			for (QNUInt32 plab = 0; plab < this->nActualLabs; plab++)
//			{
//				// just for debugging
//				cout << "plab: " << plab << ", ";
//
//				//QNUInt32 idx = plab * this->nLabs + clab;
//				QNUInt32 idx = plab * this->nActualLabs + lab;
//				//alpha_beta=expE(prev_adj_seg_alpha[plab]+this->transMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
//				alpha_beta=expE(prevAlphaPlusState[plab] + this->transMatrix[idx]+ this->betaArray[lab] - Zx);
//				alpha_beta_trans_tot+=alpha_beta;
//				//match=((clab==this->label)&&(plab==prev_lab));
//				bool match=(lab==actualLab && plab==actualPLab);
//
//				// just for debugging
//				cout << "prevAlphaPlusState[plab=" << plab << "](" << prevAlphaPlusState[plab] <<
//						") + transMatrix[plab=" << plab << "](" << this->transMatrix[idx] <<
//						") + betaArray[lab=" << lab << "](" << this->betaArray[lab] <<
//						") - Zx(" << Zx <<
//						") = " << prevAlphaPlusState[plab] + this->transMatrix[idx]+ this->betaArray[lab] - Zx << endl;
//				cout << "alpha_beta=" << alpha_beta << ", alpha_beta_trans_tot=" << alpha_beta_trans_tot << endl;
//
//				// just for debugging
//				if (match)
//					//cout << "(" << clab << "<-" <<  plab << ") ";
//					cout << "transition label match: " << "(" << actualLab << "<-" <<  actualPLab << ") " << endl;
//
//				// TODO: design a feature map in which the transition matrix calculation can take different size of plab (from prevNode) and clab (from current node).
//				//logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab);
//
//				// since there is no transition feature but only transition bias,
//				// pass NULL as the feature vector to computeTransExpF()
//				if (match)
//				{
//					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(NULL,lambda,ExpF,grad,alpha_beta,actualPLab,actualLab,plab,lab);
//				} else {
//					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(NULL,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,CRF_LAB_BAD,plab,lab);
//				}
//			}
//		}
//
//		delete[] prevAlphaPlusState;
//		delete[] tmpLogAddAcc;
//
//
////		double* stateValueLogSumOverDur = new double[this->nActualLabs];
////		double* tmpLogAddAcc = new double[this->numPrevNodes];
////		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
////		{
////			// just for debugging
////			cout << "lab: " << lab << endl;
////
////			QNUInt32 logAddID = 0;
////			for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
////			{
////				// just for debugging
////				cout << "stateArray[lab=" << lab << ", dur=" << dur << "]="
////						<< this->stateArray[this->nActualLabs * (dur - 1) + lab] << endl;
////
////				tmpLogAddAcc[logAddID++] = this->stateArray[this->nActualLabs * (dur - 1) + lab];
////			}
////			try {
////				stateValueLogSumOverDur[lab] = logAdd(tmpLogAddAcc,logAddID);
////
////				// just for debugging
////				cout << "stateValueLogSumOverDur[" << lab << "]=" << stateValueLogSumOverDur[lab] << endl;
////			}
////			catch (exception &e) {
////				string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() caught exception: "+string(e.what())+", while computing expected features";
////				throw runtime_error(errstr);
////				return(-1);
////			}
////		}
////		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
////		{
////			// just for debugging
////			cout << "lab: " << lab << endl;
////
////			// transition matrice are the same for all the previous nodes,
////			// since only transition bias is used but not transition features
////			CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - 1];
////			double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
////			for (QNUInt32 plab = 0; plab < prevAdjacentSeg->getNumAvailLabs(); plab++)
////			{
////				// just for debugging
////				cout << "plab: " << plab << ", ";
////
////				//QNUInt32 idx = plab * this->nLabs + clab;
////				QNUInt32 idx = plab * this->nActualLabs + lab;
////				//alpha_beta=expE(prev_adj_seg_alpha[plab]+this->transMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
////				alpha_beta=expE(prev_adj_seg_alpha[plab]+this->transMatrix[idx]+stateValueLogSumOverDur[lab]+this->betaArray[lab]-Zx);
////				alpha_beta_trans_tot+=alpha_beta;
////				//match=((clab==this->label)&&(plab==prev_lab));
////				bool match=(lab==actualLab && plab==actualPLab);
////
////				// just for debugging
////				cout << "alpha_beta=" << alpha_beta << ", alpha_beta_trans_tot=" << alpha_beta_trans_tot << endl;
////
////				// just for debugging
////				if (match)
////					//cout << "(" << clab << "<-" <<  plab << ") ";
////					cout << "transition label match: " << "(" << actualLab << "<-" <<  actualPLab << ") " << endl;
////
////				// TODO: design a feature map in which the transition matrix calculation can take different size of plab (from prevNode) and clab (from current node).
////				//logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab);
////
////				// since there is no transition feature but only transition bias,
////				// pass NULL as the feature vector to computeTransExpF()
////				if (match)
////				{
////					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(NULL,lambda,ExpF,grad,alpha_beta,actualPLab,actualLab,plab,lab);
////				} else {
////					logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(NULL,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,CRF_LAB_BAD,plab,lab);
////				}
////			}
////		}
////		delete[] stateValueLogSumOverDur;
////		delete[] tmpLogAddAcc;
	} else {
		// numNextNodes == 0, which means it is the last node of the sequence and
		// there are no following transitions
		// - skip the transition calculation in this case
		// but set the alpha_beta_trans_tot to 1.0 for the check below
		alpha_beta_trans_tot=1.0;
	}
//	} /*** End of parallel region ***/


	//just for debugging
//	cout << "\tAlpha_beta_tot: " << alpha_beta_tot << "\tAlpha_beta_trans_tot: " << alpha_beta_trans_tot << endl;

	// alpha_beta_tot and alpha_beta_trans_tot are no longer equal to 1 but less than 1 except the ending node.
	if ((alpha_beta_tot >1.000001))  {
		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
//	else if (alpha_beta_tot < 0.999999) {
	else if (alpha_beta_tot < -0.000001) {
//		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: Probability sums less than 0.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.000001) {
		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
//	else if (alpha_beta_trans_tot < 0.999999) {
	else if (alpha_beta_trans_tot < -0.000001) {
//		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: Trans Probability sums less than 0.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}

	// for this model, due to the way we calculate transition alpha-betas,
	// we should have alpha_beta_tot = alpha_beta_trans_tot
	if (alpha_beta_tot - alpha_beta_trans_tot > 0.000001) {
		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: State Probability sums(" + stringify(alpha_beta_tot) + ") greater than Trans Probability sum(" + stringify(alpha_beta_trans_tot) + ") for over 0.000001.";
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot - alpha_beta_trans_tot < -0.000001) {
		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: State Probability sums(" + stringify(alpha_beta_tot) + ") less than Trans Probability sum(" + stringify(alpha_beta_trans_tot) + ") for over 0.000001.";
		throw runtime_error(errstr);
	}

	// just for debugging
//	cout << endl;

	return logLi;


	////////////////////////////////////////////

//	float* seg_ftr_buf = this->ftrBuf;
////	QNUInt32 clab = 0;
//	for (QNUInt32 dur = 1; dur <= this->numPrevNodes; dur++)
//	{
//		// just for debugging
////		cout << "dur: " << dur << endl;
//
//		CRF_StateNode* prevAdjacentSeg = this->prevNodes[this->numPrevNodes - dur];
//		double* prev_adj_seg_alpha = prevAdjacentSeg->getAlpha();
////		double* stateArrayForCurDur = &(this->stateArray[this->nActualLabs * (dur - 1)]);
////		double* transMatrixForCurDur = &(this->transMatrix[this->nActualLabs * this->nActualLabs * (dur - 1)]);
//
//		// just for debugging
////		cout << "After getting prevAdjacentSeg." << endl;
//
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			// just for debugging
////			cout << "First phase:: numPrevNodes: " << numPrevNodes << ", dur: " << dur << ", lab: " << lab << endl;
//
//			//alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
//			alpha_beta=expE(this->alphaArray_WithDur[lab * this->nodeLabMaxDur + dur - 1] + this->betaArray[lab] - Zx);
//			alpha_beta_tot += alpha_beta;
//			//bool match=(clab==this->label);
//			bool match=(lab==actualLab && dur==labDur);
//
//			// just for debugging
////			cout << "alpha_beta=" << alpha_beta << ", alpha_beta_tot=" << alpha_beta_tot << endl;
//
//			// just for debugging
////			if (match)
//				//cout << clab << " ";
////				cout << "state label match: actualLab=" << actualLab << ", duration=" << labDur << endl;
//
//			//logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,this->label,clab);
//			if (match)
//			{
//				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,actualLab,lab);
//			} else {
//				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,lab);
//			}
////			//TODO: verify: this->nLabs or nLabs_of_prevNode (probably this)
////			if (prev_lab > this->nLabs) {
////				// if prev_lab > nLabs, we're in the first label frame and there are no previous
////				// transitions - skip the transition calculation in this case
////				// but set the alpha_beta_trans_tot to 1.0 for the check below
////				alpha_beta_trans_tot=1.0;
////			}
////			else {
//				// Otherwise do the transition calculations
//				//for (QNUInt32 plab=0; plab<nLabs; plab++) {
//				for (QNUInt32 plab = 0; plab < prevAdjacentSeg->getNumAvailLabs(); plab++) {
//
//					// just for debugging
////					cout << "plab: " << plab << ", ";
//
//					//QNUInt32 idx = plab * this->nLabs + clab;
//					QNUInt32 idx = plab * this->nActualLabs + lab;
//					//alpha_beta=expE(prev_adj_seg_alpha[plab]+this->transMatrix[idx]+this->stateArray[clab]+this->betaArray[clab]-Zx);
//					alpha_beta=expE(prev_adj_seg_alpha[plab]+transMatrixForCurDur[idx]+stateArrayForCurDur[lab]+this->betaArray[lab]-Zx);
//					alpha_beta_trans_tot+=alpha_beta;
//					//match=((clab==this->label)&&(plab==prev_lab));
//					match=(lab==actualLab && dur==labDur && plab==actualPLab);
//
//					// just for debugging
////					cout << "alpha_beta=" << alpha_beta << ", alpha_beta_trans_tot=" << alpha_beta_trans_tot << endl;
//
//					// just for debugging
////					if (match)
//						//cout << "(" << clab << "<-" <<  plab << ") ";
////						cout << "transition label match: " << "(" << actualLab << "<-" <<  actualPLab << ") " << endl;
//
//					// TODO: design a feature map in which the transition matrix calculation can take different size of plab (from prevNode) and clab (from current node).
//					//logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,prev_lab,this->label,plab,clab);
//					if (match)
//					{
//						logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,actualPLab,actualLab,plab,lab);
//					} else {
//						logLi+=this->crf_ptr->getFeatureMap()->computeTransExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,CRF_LAB_BAD,plab,lab);
//					}
//				}
////			}
////			clab++;
//		}
//		seg_ftr_buf += this->nFtrsPerSeg;
//	}
//	// These are the cases when the current node serves as the beginning segment of the sequence, so there is no previous node.
//	for (QNUInt32 dur = this->numPrevNodes + 1; dur <= this->nodeLabMaxDur; dur++)
//	{
//		// just for debugging
////		cout << "dur: " << dur << endl;
//
//		for (QNUInt32 lab = 0; lab < this->nActualLabs; lab++)
//		{
//			// just debugging
////			cout << "Second phase:: numPrevNodes: " << numPrevNodes << ", dur: " << dur << ", lab: " << lab << endl;
//
//			//alpha_beta=expE(this->alphaArray[clab]+this->betaArray[clab]-Zx);
//			alpha_beta=expE(this->alphaArray_WithDur[lab * this->nodeLabMaxDur + dur - 1] + this->betaArray[lab] - Zx);
//			alpha_beta_tot += alpha_beta;
//			//bool match=(clab==this->label);
//			bool match=(lab==actualLab && dur==labDur);
//
//			// just for debugging
////			cout << "alpha_beta=" << alpha_beta << ", alpha_beta_tot=" << alpha_beta_tot << endl;
//
//			// just for debugging
////			if (match)
//				//cout << clab << " ";
////				cout << "state label match: actualLab=" << actualLab << ", duration=" << labDur << endl;
//
//
//			//logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,this->label,clab);
//			if (match)
//			{
//				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,actualLab,lab);
//			} else {
//				logLi+=this->crf_ptr->getFeatureMap()->computeStateExpF(seg_ftr_buf,lambda,ExpF,grad,alpha_beta,CRF_LAB_BAD,lab);
//			}
//
////			clab++;
//		}
//		seg_ftr_buf += this->nFtrsPerSeg;
//	}
//
//	if (this->numPrevNodes == 0)
//	{
//		// if numPrevNodes == 0, there are no previous transitions
//		// - skip the transition calculation in this case
//		// but set the alpha_beta_trans_tot to 1.0 for the check below
//		alpha_beta_trans_tot=1.0;
//	}
//
//	//just for debugging
////	cout << "\tAlpha_beta_tot: " << alpha_beta_tot << "\tAlpha_beta_trans_tot: " << alpha_beta_trans_tot << endl;
//
//	// alpha_beta_tot and alpha_beta_trans_tot are no longer equal to 1 but less than 1 except the ending node.
//	if ((alpha_beta_tot >1.000001))  {
//		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
//		throw runtime_error(errstr);
//	}
////	else if (alpha_beta_tot < 0.999999) {
//	else if (alpha_beta_tot < -0.000001) {
////		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
//		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: Probability sums less than 0.0 "+stringify(alpha_beta_tot);
//		throw runtime_error(errstr);
//	}
//	else if (alpha_beta_trans_tot > 1.000001) {
//		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
//		throw runtime_error(errstr);
//	}
////	else if (alpha_beta_trans_tot < 0.999999) {
//	else if (alpha_beta_trans_tot < -0.000001) {
////		string errstr="CRF_StdSegStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
//		string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeExpF() threw exception: Trans Probability sums less than 0.0 "+stringify(alpha_beta_trans_tot);
//		throw runtime_error(errstr);
//	}
//	return logLi;
}

// only specific to CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr
/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlphaPlusTrans(double* nextNodeTransMatrix)
 *
 * add transition values (from current node to the next node) to the alphas for future use in next nodes.
 *
 * nextNodeTransMatrix: the transition matrix of the next node
 *
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlphaPlusTrans(double* nextNodeTransMatrix)
{
	double result = 0.0;

	// add transition values (from current node to the next node) to the alphas for future use in next nodes.

	// just for debugging
//	cout << "adding transition values to the computed log-alphas for future use in the next node..." << endl;

	QNUInt32 logAddAccSize = this->nActualLabs;
	double* tmpLogAddAcc = new double[logAddAccSize];

	// assuming the number of labels in the next node is equal to that in the current node
	for (QNUInt32 next_lab = 0; next_lab < this->nActualLabs; next_lab++)
	{
		QNUInt32 logAddID = 0;
		for (QNUInt32 clab = 0; clab < this->nActualLabs; clab++)
		{
			tmpLogAddAcc[logAddID++] = this->alphaArray[clab] +
					nextNodeTransMatrix[clab * this->nActualLabs + next_lab];

			// just for debugging
//			cout << "alphaArray[clab=" << clab << "](" << this->alphaArray[clab] << ") + " <<
//					"nextNodeTransMatrix[clab=" << clab << ", next_lab=" << next_lab << "](" << nextNodeTransMatrix[clab * this->nActualLabs + next_lab] <<
//					") = " << tmpLogAddAcc[logAddID-1] << endl;
		}
		try {
			this->alphaPlusTrans[next_lab] = logAdd(tmpLogAddAcc,logAddID);

			// just for debugging
//			cout << "alphaPlusTrans[next_lab=" << next_lab << "]=" << this->alphaPlusTrans[next_lab] << endl;
		}
		catch (exception &e) {
			string errstr="CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::computeAlphaPlusTrans() caught exception: "+string(e.what())+", while computing alpha plus trans";
			throw runtime_error(errstr);
			return(-1);
		}
	}

	delete [] tmpLogAddAcc;

	return result;
}

// only specific to CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr
/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getAlphaPlusTrans(QNUInt32 next_lab)
 *
 * return the summation of [alpha(cur_lab) + trans_value(cur_lab, next_lab)]
 * over all cur_lab for the given next_lab
 *
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getAlphaPlusTrans(QNUInt32 next_lab)
{
	return this->alphaPlusTrans[next_lab];
}


// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
 *
 * Correct getTransValue() function for CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	return this->transMatrix[prev_lab * this->nActualLabs + cur_lab];
}

// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 cur_dur)
 *
 * Correct getFullTransValue() function for CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 cur_dur)
{
	return getTransValue(prev_lab, cur_lab) + getStateValue(cur_lab, cur_dur);
}

// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getTempBeta(QNUInt32 next_lab, QNUInt32 next_dur)
 *
 * Disabled. Not applicable in CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getTempBeta(QNUInt32 next_lab, QNUInt32 next_dur)
{
	string errstr="Error: getTempBeta() is not applicable in CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr.";
	throw runtime_error(errstr);
}

// Disable all these functions by overriding them with exception handling. Use their modified version below.
// Added by Ryan
/*
 * CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 cur_dur)
 *
 * Disabled in CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr, use getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab) instead.
 *
 */
double CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 cur_dur)
{
	string errstr="Error: use CRF_StdSegStateNode_WithoutDurLab_WithoutSegTransFtr::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab) instead.";
	throw runtime_error(errstr);
}
