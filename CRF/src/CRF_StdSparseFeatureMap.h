#ifndef CRF_STDSPARSEFEATUREMAP_H_
#define CRF_STDSPARSEFEATUREMAP_H_

#include "CRF_FeatureMap.h"
#include "CRF_StdFeatureMap.h"

class CRF_StdSparseFeatureMap : public CRF_StdFeatureMap
{
public:
	CRF_StdSparseFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas);
	virtual ~CRF_StdSparseFeatureMap();
	virtual double computeRi(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 clab);
	virtual double computeMij(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 plab, QNUInt32 clab);
	virtual double computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab);
	virtual double computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab);
	virtual double computeExpFState(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab);
	virtual double computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab, QNUInt32 plab);
	virtual double computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab);
	virtual double computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab);
};

#endif /*CRF_STDFEATUREMAP_H_*/
