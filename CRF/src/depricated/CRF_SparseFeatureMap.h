#ifndef CRF_SPARSEFEATUREMAP_H_
#define CRF_SPARSEFEATUREMAP_H_

#include <vector>

#include "CRF_FeatureMap.h"

class CRF_SparseFeatureMap : public CRF_FeatureMap
{
protected:
	vector< vector<QNUInt32> > stateVector;
public:
	CRF_SparseFeatureMap();
	virtual ~CRF_SparseFeatureMap();
};

#endif /*CRF_SPARSEFEATUREMAP_H_*/
