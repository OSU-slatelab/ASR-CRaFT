#ifndef CRF_STDSTATENODE_H_
#define CRF_STDSTATENODE_H_

#include "CRF_StateNode.h"

class CRF_StdStateNode : public CRF_StateNode
{
protected:
	double* stateArray;
	double* transMatrix;
	double* tempBeta;
public:
	CRF_StdStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdStateNode();
	virtual double computeTransMatrix();
	virtual double computeTransMatrixLog();
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
	virtual double* computeAlphaBeta(double Zx);
	virtual double computeAlphaSum();
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	virtual double getStateValue(QNUInt32 cur_lab);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
};

#endif /*CRF_STDSTATENODE_H_*/
