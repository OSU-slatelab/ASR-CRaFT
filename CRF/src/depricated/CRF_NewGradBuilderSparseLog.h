#ifndef CRF_NEWGRADBUILDERLOG_H_
#define CRF_NEWGRADBUILDERLOG_H_

#include "CRF_NewGradBuilder.h"
#include "CRF_StdStateNodeLog.h"
#include "CRF_StdStateVectorLog.h"

class CRF_NewGradBuilderSparseLog : public CRF_NewGradBuilder
{
public:
	CRF_NewGradBuilderSparseLog(CRF_Model* crf_in);
	virtual ~CRF_NewGradBuilderSparseLog();
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
};

#endif /*CRF_NEWGRADBUILDERLOG_H_*/
