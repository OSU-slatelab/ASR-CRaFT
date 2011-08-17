#ifndef CRF_STDNSTATEVECTOR_H_
#define CRF_STDNSTATEVECTOR_H_

#include "CRF_StateVector.h"
#include "CRF_StdNStateNode.h"

class CRF_StdNStateVector : public CRF_StateVector
{
public:
	CRF_StdNStateVector();
	virtual ~CRF_StdNStateVector();
	virtual void set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in);
};

#endif /*CRF_STDNSTATEVECTOR_H_*/
