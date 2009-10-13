#ifndef CRF_NEWGRADBUILDERLOG_H_
#define CRF_NEWGRADBUILDERLOG_H_

#include "CRF_NewGradBuilder.h"
#include "CRF_StdStateNodeLog.h"
#include "CRF_StdStateVectorLog.h"

class CRF_NewGradBuilderLog : public CRF_NewGradBuilder
{
public:
	CRF_NewGradBuilderLog(CRF_Model* crf_in);
	virtual ~CRF_NewGradBuilderLog();
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
};

#endif /*CRF_NEWGRADBUILDERLOG_H_*/
