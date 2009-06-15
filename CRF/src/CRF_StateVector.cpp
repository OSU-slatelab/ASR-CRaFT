#include "CRF_StateVector.h"

CRF_StateVector::CRF_StateVector()
{
}

CRF_StateVector::~CRF_StateVector()
{
}

void CRF_StateVector::set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in)
{
}

void CRF_StateVector::setNodeCount(QNUInt32 cnt)
{
	this->nodeCount=cnt;
}

QNUInt32 CRF_StateVector::getNodeCount()
{
	return this->nodeCount;
}