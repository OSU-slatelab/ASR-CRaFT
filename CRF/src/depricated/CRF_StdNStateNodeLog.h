#ifndef CRF_STDNSTATENODELOG_H_
#define CRF_STDNSTATENODELOG_H_

#include "CRF_StdNStateNode.h"

class CRF_StdNStateNodeLog : public CRF_StdNStateNode
{
protected:
	double* logAddAcc;
public:
	CRF_StdNStateNodeLog(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdNStateNodeLog();
	virtual double computeTransMatrix();
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
	virtual double* computeAlphaBeta(double Zx);
	virtual double computeAlphaSum();
	virtual double computeAlphaAlignedSum();
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
};

#endif /*CRF_STDNSTATENODELOG_H_*/
