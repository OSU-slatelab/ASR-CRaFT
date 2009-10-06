#ifndef CRF_NEWLOCALPOSTERIORBUILDER_H_
#define CRF_NEWLOCALPOSTERIORBUILDER_H_

#include "CRF_Model.h"
#include "CRF_StateVector.h"
#include "CRF_StdStateVectorLog.h"
#include "CRF_StdNStateVectorLog.h"
#include "CRF_FeatureStream.h"

class CRF_NewLocalPosteriorBuilder
{
protected:
	CRF_Model* crf;
	CRF_StateVector* nodeList;
	bool ownsNodeList;
	float* ftr_buf;
	QNUInt32* lab_buf;
	double* alpha_base;
	bool normalize;
public:
	CRF_NewLocalPosteriorBuilder(CRF_Model* crf_in, bool norm=true);
	virtual ~CRF_NewLocalPosteriorBuilder();
	virtual CRF_StateVector* buildFtrSeq(CRF_FeatureStream* ftr_strm);
	virtual void computeAlphaBeta(CRF_StateVector *nodeList);
	virtual CRF_StateVector* buildFtrSeqNState(CRF_FeatureStream* ftr_strm);
};

#endif /*CRF_NEWLOCALPOSTERIORBUILDER_H_*/
