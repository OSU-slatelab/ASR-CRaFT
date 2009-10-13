#ifndef CRF_LOCALPOSTERIORBUILDER_H_
#define CRF_LOCALPOSTERIORBUILDER_H_

#include "CRF.h"
#include "CRF_Model.h"
#include "CRF_StdGradBuilderLog.h"

class CRF_LocalPosteriorBuilder : public CRF_StdGradBuilderLog
{
public:
	CRF_LocalPosteriorBuilder(CRF_Model* crf_in);
	virtual ~CRF_LocalPosteriorBuilder();
	virtual CRF_Seq* buildFtrSeq(CRF_FeatureStream *ftr_strm);
};

#endif /*CRF_LOCALPOSTERIOR_H_*/
