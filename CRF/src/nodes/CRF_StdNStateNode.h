#ifndef CRF_STDNSTATENODE_H_
#define CRF_STDNSTATENODE_H_
/*
 * CRF_StdNStateNode.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_StdNStateNode
 */

#include "CRF_StateNode.h"
/*
 * class CRF_StdNStateNode
 *
 * Used in training and decoding processing.  Holds the features for a given state, the label (if during
 * training), the computed alpha and beta arrays, and the computed transition matrix.
 *
 * These nodes are for a restricted multi-state crf topology.  The number of states can vary, but must
 * be the same for all labels.  In addition, the states must be passed through sequentially with no
 * states skipped (i.e. a three state model must have a minimum of 3 frames of features, each assigned
 * to a different state).
 */
class CRF_StdNStateNode : public CRF_StateNode
{
protected:
	double* stateArray;
	double* denseTransMatrix;
	double* diagTransMatrix;
	double* offDiagTransMatrix;  // indexed by the previous state
	QNUInt32 nStates;
	QNUInt32 nFullLabs;  // nFullLabs = nLabs / nStates
	double* tempBeta;
	double* logAddAcc;
public:
	CRF_StdNStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdNStateNode();
	virtual double computeTransMatrix();
//	virtual double computeTransMatrixLog();
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
	virtual double* computeAlphaBeta(double Zx);
	virtual double computeAlphaSum();
	virtual double computeAlphaAlignedSum();
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	virtual double getStateValue(QNUInt32 cur_lab);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);

};

#endif /*CRF_STDNSTATENODE_H_*/
