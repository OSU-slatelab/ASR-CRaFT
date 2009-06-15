#ifndef CRF_STDSTATEVECTOR_H_
#define CRF_STDSTATEVECTOR_H_

#include "CRF_StateVector.h"
#include "CRF_StdStateNode.h"

class CRF_StdStateVector : public CRF_StateVector
{
public:
	CRF_StdStateVector();
	virtual ~CRF_StdStateVector();
	virtual void set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in);
};

#endif /*CRF_STDSTATEVECTOR_H_*/
