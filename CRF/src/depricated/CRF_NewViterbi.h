#ifndef CRF_NEWVITERBI_H_
#define CRF_NEWVITERBI_H_
#include <vector>

#include "CRF.h"

#include "CRF_Model.h"
#include "CRF_FeatureStream.h"
#include "CRF_StdStateVectorLog.h"
#include "CRF_StdNStateVectorLog.h"
#include "CRF_LabelPath.h"

class CRF_NewViterbi
{
protected:
	CRF_StateVector* nodeList;
	CRF_Model* crf;
	CRF_FeatureStream* ftr_strm;
	float* ftr_buf;
	QNUInt32* lab_buf;
	QNUInt32 bunch_size;
	QNUInt32 num_ftrs;
	QNUInt32 num_labs;
public:
	CRF_NewViterbi(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_NewViterbi();
	virtual CRF_LabelPath* bestPath();
	virtual CRF_LabelPath* alignPath();
	virtual CRF_LabelPath* bestPathVec();
	virtual CRF_LabelPath* nStateBestPathVec();
};

#endif /*CRF_NEWVITERBI_H_*/
