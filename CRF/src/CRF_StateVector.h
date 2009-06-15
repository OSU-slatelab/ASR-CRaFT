#ifndef CRF_STATEVECTOR_H_
#define CRF_STATEVECTOR_H_
#include <vector>

#include "CRF.h"
#include "CRF_Model.h"
#include "CRF_StateNode.h"

class CRF_StateVector : public vector <CRF_StateNode*>
{
private:
	QNUInt32 nodeCount;
public:
	CRF_StateVector();
	virtual ~CRF_StateVector();
	virtual void set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in);
	virtual void setNodeCount(QNUInt32 cnt);
	virtual QNUInt32 getNodeCount();
};

#endif /*CRF_STATEVECTOR_H_*/
