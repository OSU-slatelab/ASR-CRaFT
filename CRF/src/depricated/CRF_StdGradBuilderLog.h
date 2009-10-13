#ifndef CRF_STDGRADBUILDERLOG_H_
#define CRF_STDGRADBUILDERLOG_H_

#include "CRF.h"
#include "CRF_StdGradBuilder.h"

using namespace std;
using namespace CRF_LogMath;

class CRF_StdGradBuilderLog : public CRF_StdGradBuilder
{
protected:
	double* logAddAcc;
public:
	CRF_StdGradBuilderLog(CRF_Model* crf_in);
	//virtual ~CRF_StdGradBuilderLog();
	virtual double buildGradient(CRF_FeatureStream *ftr_strm, double* grad);
	virtual double computeAlphaLog(CRF_Seq* cur_seq);
	virtual double computeBetaLog(CRF_Seq* cur_seq);
	virtual double computeExpFLog(CRF_Seq* cur_seq, double* grad, double Zx);
};

#endif /*CRF_STDGRADBUILDERLOG_H_*/
