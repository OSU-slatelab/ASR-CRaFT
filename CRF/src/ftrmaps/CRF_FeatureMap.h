#ifndef CRF_FEATUREMAP_H_
#define CRF_FEATUREMAP_H_

#include "../CRF.h"

#ifndef QN_UINT32_MAX
#define QN_UINT32_MAX (0xffffffff)
#endif

struct CRF_FeatureMap_config {
	ftrmaptype map_type;
	QNUInt32 numLabs;
	QNUInt32 numFeas;
	QNUInt32 numStates;
	bool useStateFtrs;
	QNUInt32 stateFidxStart;
	QNUInt32 stateFidxEnd;
	bool useTransFtrs;
	QNUInt32 transFidxStart;
	QNUInt32 transFidxEnd;
	bool useStateBias;
	bool useTransBias;
	double stateBiasVal;
	double transBiasVal;
};

class CRF_FeatureMap
{
protected:
	CRF_FeatureMap_config* config;
	QNUInt32 numFtrFuncs;
public:
	CRF_FeatureMap(QNUInt32 nlabs, QNUInt32 nfeas);
	CRF_FeatureMap(CRF_FeatureMap_config* cnf);
	virtual ~CRF_FeatureMap();
	static CRF_FeatureMap* createFeatureMap(CRF_FeatureMap_config* cnf);

	/*virtual double computeRi(float* ftr_buf, double* lambda, QNUInt32 &lc, QNUInt32 clab);
	virtual double computeMij(float* ftr_buf, double* lambda, QNUInt32 &lc, QNUInt32 plab, QNUInt32 clab);
	virtual double computeExpFState(float* ftr_buf, double* lambda, QNUInt32 &lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab);
	virtual double computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32 &lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab);
	*/

	virtual double computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab);
	virtual double computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab);
	virtual double computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab, bool compute_grad=true);
	virtual double computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad=true);
	virtual QNUInt32 getNumFtrFuncs();
	virtual QNUInt32 getNumStates();
	virtual QNUInt32 getNumStateFuncs(QNUInt32 clab);
	virtual QNUInt32 getNumTransFuncs(QNUInt32 plab, QNUInt32 clab);
	virtual QNUInt32 getStateFeatureIdx(QNUInt32 clab, QNUInt32 fno);
	virtual QNUInt32 getTransFeatureIdx(QNUInt32 plab, QNUInt32 clab, QNUInt32 fno);
	virtual QNUInt32 getStateBiasIdx(QNUInt32 clab);
	virtual QNUInt32 getTransBiasIdx(QNUInt32 plab, QNUInt32 clab);
	virtual void setNumStates(QNUInt32 ns);
	virtual QNUInt32 recalc();
	virtual string getMapDescriptor(QNUInt32 lambdaNum);
	virtual void accumulateFeatures(float *ftr_buf, double *accumulator, QNUInt32 lab);
};

#endif /*CRF_FEATUREMAP_H_*/