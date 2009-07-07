#ifndef CRF_FEATUREMAP_H_
#define CRF_FEATUREMAP_H_

#include "CRF.h"

class CRF_FeatureMap
{
protected:
	QNUInt32 numFtrFuncs;
	QNUInt32 numStates;
public:
	CRF_FeatureMap();
	CRF_FeatureMap(QNUInt32 nlabs, QNUInt32 nfeas);
	virtual ~CRF_FeatureMap();
	virtual double computeRi(float* ftr_buf, double* lambda, QNUInt32 &lc, QNUInt32 clab);
	virtual double computeMij(float* ftr_buf, double* lambda, QNUInt32 &lc, QNUInt32 plab, QNUInt32 clab);
	virtual double computeExpFState(float* ftr_buf, double* lambda, QNUInt32 &lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab);
	virtual double computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32 &lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab);
	virtual QNUInt32 getNumFtrFuncs();
	virtual QNUInt32 getNumStates();
	virtual QNUInt32 getNumStateFuncs(QNUInt32 clab);
	virtual QNUInt32 getNumTransFuncs(QNUInt32 plab, QNUInt32 clab);
	virtual void setNumStates(QNUInt32 ns);
	virtual QNUInt32 recalc();
	virtual string getMapDescriptor(QNUInt32 lambdaNum);
};

#endif /*CRF_FEATUREMAP_H_*/
