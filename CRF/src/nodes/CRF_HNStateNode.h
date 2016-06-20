#ifndef CRF_HNSTATENODE_H_
#define CRF_HNSTATENODE_H_

#include "CRF_StateNode.h"

class CRF_HNStateNode : public CRF_StateNode
{
protected:
	double* stateArray;
	double* denseTransMatrix;
	double* diagTransMatrix;
	double* offDiagTransMatrix;
	QNUInt32 nStates;
	QNUInt32 trueNLabs;
	double* tempBeta;
	double* logAddAcc;
public:
	CRF_HNStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_HNStateNode();
	virtual double computeTransMatrix();
//	virtual double computeTransMatrixLog();
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
	virtual double* computeAlphaBeta(double Zx);
	virtual double computeAlphaSum();
	virtual double computeAlphaAligned(double* prev_alpha, QNUInt32 prev_label, QNUInt32 next_label);
	virtual double computeAlphaAlignedSum();
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	virtual double getStateValue(QNUInt32 cur_lab);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);

};

#endif /*CRF_STDNSTATENODE_H_*/
