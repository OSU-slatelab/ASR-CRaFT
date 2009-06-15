#ifndef CRF_STDTRANSFEATUREMAP_H_
#define CRF_STDTRANSFEATUREMAP_H_

#include "CRF_StdFeatureMap.h"

class CRF_StdTransFeatureMap : public CRF_StdFeatureMap
{
protected:
	QNUInt32 transFidxStart;
	QNUInt32 transFidxEnd;
public:
	CRF_StdTransFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas);
	virtual ~CRF_StdTransFeatureMap();
	virtual double computeMij(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 plab, QNUInt32 clab);
	virtual double computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab);
	virtual void setTransFtrStart(QNUInt32 st);
	virtual void setTransFtrEnd(QNUInt32 end);
	virtual QNUInt32 init();	
};

#endif /*CRF_STDTRANSFEATUREMAP_H_*/
