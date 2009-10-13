#ifndef CRF_STDGRADBUILDER_H_
#define CRF_STDGRADBUILDER_H_

#include "CRF.h"
#include "CRF_Seq.h"
#include "CRF_GradBuilder.h"

using namespace std;
using namespace CRF_LogMath;

class CRF_StdGradBuilder : public CRF_GradBuilder
{
public:
	CRF_StdGradBuilder(CRF_Model* crf_in);
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad);
	virtual double computeTransMatrix(CRF_Seq* cur_seq);
	virtual double computeTransMatrixLog(CRF_Seq* cur_seq);
	virtual double computeAlpha(CRF_Seq* cur_seq, CRF_Seq* prev_seq);
	virtual double computeBeta(CRF_Seq* cur_seq, CRF_Seq* next_seq);
	virtual double computeExpF(CRF_Seq* cur_seq, double* grad, double Zx);
};

#endif /*CRF_STDGRADBUILDER_H_*/
