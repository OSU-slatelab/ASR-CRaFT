#include "CRF_StdNStateVectorLog.h"

CRF_StdNStateVectorLog::CRF_StdNStateVectorLog()
{
}

CRF_StdNStateVectorLog::~CRF_StdNStateVectorLog()
{
}

void CRF_StdNStateVectorLog::set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in)
{
	if (idx >= this->size() ) {
		this->push_back(new CRF_StdNStateNodeLog(new_buf,num_ftrs,lab_buf,crf_in));
	}
	else {
		this->at(idx)->reset(new_buf,num_ftrs,lab_buf,crf_in);
	}
}
