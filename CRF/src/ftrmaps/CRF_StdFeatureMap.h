#ifndef CRF_STDFEATUREMAP_H_
#define CRF_STDFEATUREMAP_H_

#include "CRF_FeatureMap.h"

class CRF_StdFeatureMap : public CRF_FeatureMap
{
protected:
	/*QNUInt32 numLabs;
	QNUInt32 numFeas;
	QNUInt32 numStates;
	bool useStateBias;
	bool useTransBias;
	bool useStateFtrs;
	bool useTransFtrs;
	QNUInt32 stateFidxStart;
	QNUInt32 stateFidxEnd;
	QNUInt32 transFidxStart;
	QNUInt32 transFidxEnd;
	double stateBiasVal;
	double transBiasVal;*/

	QNUInt32 numStateFuncs;
	QNUInt32 numTransFuncs;
	QNUInt32 numActualLabels;
	QNUInt32 transMult;
	QNUInt32* stateFeatureIdxCache;
	QNUInt32* transFeatureIdxCache;

	virtual QNUInt32 computeStateFeatureIdx(QNUInt32 clab, QNUInt32 fno=0);
	virtual QNUInt32 computeTransFeatureIdx(QNUInt32 clab, QNUInt32 plab, QNUInt32 fno=0);

public:
	CRF_StdFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas);
	CRF_StdFeatureMap(CRF_FeatureMap_config* cnf);
	virtual ~CRF_StdFeatureMap();

	/*virtual double computeRi(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 clab);
	virtual double computeMij(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 plab, QNUInt32 clab);
	virtual double computeExpFState(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab);
	virtual double computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab, QNUInt32 plab);
	*/

	virtual double computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab);
	virtual double computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab);
	virtual double computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab, bool compute_grad=true);
	virtual double computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad=true);

	/*virtual void setStateFtrRange(QNUInt32 st, QNUInt32 end);
	virtual void setTransFtrRange(QNUInt32 st, QNUInt32 end);
	virtual void setNumStates(QNUInt32 ns);
	virtual void setUseStateBias(bool useState);
	virtual void setUseTransBias(bool useTrans);
	virtual void setUseStateFtrs(bool useState);
	virtual void setUseTransFtrs(bool useTrans);
	virtual void setStateBiasVal(double stateBias);
	virtual void setTransBiasVal(double transBias);*/

	virtual QNUInt32 getNumStateFuncs(QNUInt32 clab);
	virtual QNUInt32 getNumTransFuncs(QNUInt32 plab, QNUInt32 clab);
	virtual QNUInt32 getNumStates();
	virtual QNUInt32 getStateFeatureIdx(QNUInt32 clab, QNUInt32 fno=0);
	virtual QNUInt32 getTransFeatureIdx(QNUInt32 clab, QNUInt32 plab, QNUInt32 fno=0);
	virtual QNUInt32 getStateBiasIdx(QNUInt32 clab);
	virtual QNUInt32 getTransBiasIdx(QNUInt32 plab, QNUInt32 clab);
	virtual QNUInt32 recalc();
	virtual string getMapDescriptor(QNUInt32 lambdaNum);
	virtual void accumulateFeatures(float *ftr_buf, double *accumulator, QNUInt32 lab);
};

#endif /*CRF_STDFEATUREMAP_H_*/
