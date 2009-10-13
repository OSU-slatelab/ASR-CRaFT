#ifndef CRF_STDSTATENODELOG_H_
#define CRF_STDSTATENODELOG_H_

#include "CRF_StdStateNode.h"

class CRF_StdStateNodeLog : public CRF_StdStateNode
{
protected:
	double* logAddAcc;
public:
	CRF_StdStateNodeLog(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdStateNodeLog();
	virtual double computeTransMatrix();
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual double* computeAlphaBeta(double Zx);
	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
	virtual double computeAlphaSum();
	virtual double computeAlphaAlignedSum();
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
};

#endif /*CRF_STDSTATENODELOG_H_*/
