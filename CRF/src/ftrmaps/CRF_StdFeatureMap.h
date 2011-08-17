#ifndef CRF_STDFEATUREMAP_H_
#define CRF_STDFEATUREMAP_H_
/*
 * CRF_StdFeatureMap.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */

#include "CRF_FeatureMap.h"

/*
 * class CRF_StdFeatureMap
 * subclass of CRF_FeatureMap
 *
 * Implements a very basic "standard" FeatureMap for CRF processing.  This
 * map assumes that every label will be associated with every input feature.
 * The association between features and labels or transitions can be turned on
 * and off, but only at a high level (i.e. a CRF can be created with no transition
 * features, or no state features, but it cannot be created where only one label
 * has state features).
 *
 * Even with these restrictions, there is flexibility.  States can have feature
 * values or not.  States can have bias values or not.  Transitions can have
 * feature values or not.  Transitions can have bias values or not.  All of
 * these options are controlled via the configuration struct CRF_FeatureMap_config.
 *
 */
class CRF_StdFeatureMap : public CRF_FeatureMap
{
protected:


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

	virtual double computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab);
	virtual double computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab);
	virtual double computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab, bool compute_grad=true);
	virtual double computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad=true);

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
