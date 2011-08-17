#include "CRF_LabelPath.h"

CRF_LabelPath::CRF_LabelPath(QNUInt32 length)
	: len(length)
{
	this->path = new QNUInt32[length];
}

CRF_LabelPath::~CRF_LabelPath()
{
	delete[] this->path;
}

QNUInt32* CRF_LabelPath::getLabelPath()
{
	return this->path;
}

QNUInt32 CRF_LabelPath::getLength()
{
	return this->len;
}
