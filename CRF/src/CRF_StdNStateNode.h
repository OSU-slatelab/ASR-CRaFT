#ifndef CRF_STDNSTATENODE_H_
#define CRF_STDNSTATENODE_H_

#include "CRF_StateNode.h"

class CRF_StdNStateNode : public CRF_StateNode
{
protected:
	double* stateArray;
	double* denseTransMatrix;
	double* diagTransMatrix;
	double* offDiagTransMatrix;
	QNUInt32 nStates;
	QNUInt32 nFullLabs;
	double* tempBeta;
public:
	CRF_StdNStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdNStateNode();
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

#endif /*CRF_STDNSTATENODE_H_*/
