#ifndef CRF_STDSTATEVECTORLOG_H_
#define CRF_STDSTATEVECTORLOG_H_

#include "CRF_StateVector.h"
#include "CRF_StdStateNodeLog.h"

class CRF_StdStateVectorLog : public CRF_StateVector
{

public:
	CRF_StdStateVectorLog();
	virtual ~CRF_StdStateVectorLog();
	virtual void set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in);
};

#endif /*CRF_STATEVECTORLOG_H_*/
