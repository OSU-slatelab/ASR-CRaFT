#ifndef CRF_LABELPATH_H_
#define CRF_LABELPATH_H_

#include "CRF.h"

class CRF_LabelPath
{
private:
	QNUInt32* path;
	QNUInt32 len;
public:
	CRF_LabelPath(QNUInt32 length);
	virtual ~CRF_LabelPath();
	virtual QNUInt32* getLabelPath();
	virtual QNUInt32 getLength();
};

#endif /*CRF_LABELPATH_H_*/
