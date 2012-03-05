#ifndef CRF_STATENODE_H_
#define CRF_STATENODE_H_
/*
 * CRF_StateNode.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_StateNode
 */

#include <string>
#include "fst/fstlib.h"
#include "../CRF.h"
#include "../CRF_Model.h"
#include <vector>

using namespace fst;
/*
 * class CRF_StateNode
 *
 * Used in training and decoding processing.  Holds the features for a given state, the label (if during
 * training), the computed alpha and beta arrays, and the computed transition matrix.
 *
 * This class is a base class for CRF nodes, and should not be used directly.  See the subclasses
 * CRF_StdStateNode and CRF_StdNStateNode for examples of implementation.
 */
class CRF_StateNode
{
protected:
	float* ftrBuf;
	QNUInt32 ftrBuf_size;
	QNUInt32 label;
	CRF_Model* crf_ptr;
	double* alphaArray;
	double* betaArray;
	double* alphaBetaArray;
	vector<double> alphaArrayAligned;
	vector<double> alphaArrayAlignedBase;
	vector<double> betaArrayAligned;
	vector<double> betaArrayAlignedBase;
	QNUInt32 alphaSize;
	double alphaScale;
	QNUInt32 nLabs;
	double* prevAlpha;

	// Added by Ryan
	CRF_StateNode** prevNodes;
	CRF_StateNode** nextNodes;
	QNUInt32 numPrevNodes;
	QNUInt32 numNextNodes;
	QNUInt32 numAvailLabs;

public:
	vector<uint> viterbiPhnIds;
	vector<int> viterbiPointers;

	// Added by Ryan
	vector<bool> isPhoneStartBoundary;
	vector<uint> viterbiDurs;

	CRF_StateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StateNode();
	virtual double computeTransMatrix();
	virtual double computeTransMatrixLog();
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual double* computeAlphaBeta(double Zx);
	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
	virtual double computeSoftExpF(double* ExpF, double* grad, double Zx, double soft_Zx, double* prev_alpha, vector<double>* prevAlphaAligned, bool firstFrame);
	virtual double computeAlphaSum();
	virtual double computeAlphaAlignedSum();
	virtual void reset(float *fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf_in);
	virtual double* getAlpha();
	virtual double* getPrevAlpha();
	virtual double* getBeta();
	virtual double* getAlphaBeta();
	virtual vector<double>* getAlphaAligned();
	virtual vector<double>* getBetaAligned();
	virtual vector<double>* getAlphaAlignedBase();
	virtual vector<double>* getBetaAlignedBase();
	virtual QNUInt32 getLabel();
	virtual double getAlphaScale();
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	virtual double getStateValue(QNUInt32 cur_lab);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	static CRF_StateNode* createStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual float *getFtrBuffer();
	virtual QNUInt32 getFtrBufferSize();

	// Added by Ryan
	static CRF_StateNode* createStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab,
			CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs,
			QNUInt32 nextNode_nActualLabs);
	virtual double computeAlpha();
	virtual double computeFirstAlpha();
	virtual double computeBeta(double scale=1.0);
	virtual double computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 prev_lab);
	virtual void setPrevNodes(CRF_StateNode** prev_nodes_ptr, QNUInt32 numPrevNodes);
	virtual void setNextNodes(CRF_StateNode** next_nodes_ptr, QNUInt32 numNextNodes);
	virtual CRF_StateNode** getPrevNodes();
	virtual CRF_StateNode** getNextNodes();
	virtual QNUInt32 getNLabs();
	virtual QNUInt32 getNumAvailLabs();
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur);
	virtual double getStateValue(QNUInt32 cur_lab, QNUInt32 dur);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 dur);
	virtual double getTempBeta(QNUInt32 cur_lab, QNUInt32 dur);
//	virtual void deleteFtrBuf();
};

#endif /*CRF_STATENODE_H_*/
