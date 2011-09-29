#ifndef CRF_FEATUREMAP_H_
#define CRF_FEATUREMAP_H_
/*
 * CRF_FeatureMap.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_FeatureMap_config and CRF_FeatureMap
 */

#include "../CRF.h"

#ifndef QN_UINT32_MAX
#define QN_UINT32_MAX (0xffffffff)
#endif
/*
 * struct CRF_FeatureMap_config
 *
 * Used by the CRF_FeatureMap object to control settings.  Generally used in the
 * constructor for CRF_FeatureMap.
 *
 */
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

	//Added by Ryan, for parameter tying and segmental CRFs
	//the maximum duration if labels are phone-duration combination
	QNUInt32 maxDur;
	//the start index of duration features (binary coded, one-hot features) if any
	QNUInt32 durFtrStart;
	//the number of actual labels (without duration)
	QNUInt32 nActualLabs;
};
/*
 * class CRF_FeatureMap
 *
 * Controls the interaction of features taken from the feature stream and the
 * CRF_Model itself.  Determines the association between weights in the lambda
 * vector of the CRF and the feature functions.
 *
 * This object should be treated as merely an interface class for an actual
 * FeatureMap object and should never be instantiated on its own.  The factory
 * function "createFeatureMap" should be updated to account for new FeatureMap
 * objects as they are created.
 *
 */
class CRF_FeatureMap
{
protected:
	CRF_FeatureMap_config* config;
	QNUInt32 numFtrFuncs;

	//Added by Ryan
	virtual void tieGradForSingleParam(double* grad, QNUInt32 numParam, QNUInt32 start, QNUInt32 step);

public:
	CRF_FeatureMap(QNUInt32 nlabs, QNUInt32 nfeas);
	CRF_FeatureMap(CRF_FeatureMap_config* cnf);
	virtual ~CRF_FeatureMap();
	static CRF_FeatureMap* createFeatureMap(CRF_FeatureMap_config* cnf);

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

	//Added by Ryan
	//virtual void tieGradient(double* grad, QNUInt32 maxDur, QNUInt32 durFtrStart);
	virtual void tieGradient(double* grad);
};

#endif /*CRF_FEATUREMAP_H_*/
