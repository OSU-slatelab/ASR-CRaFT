#ifndef CRF_STDNSTATEVECTORLOG_H_
#define CRF_STDNSTATEVECTORLOG_H_

#include "CRF_StdNStateVector.h"
#include "CRF_StdNStateNodeLog.h"

class CRF_StdNStateVectorLog : public CRF_StdNStateVector
{
public:
	CRF_StdNStateVectorLog();
	virtual ~CRF_StdNStateVectorLog();
	virtual void set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in);
};

#endif /*CRF_STDNSTATEVECTORLOG_H_*/
