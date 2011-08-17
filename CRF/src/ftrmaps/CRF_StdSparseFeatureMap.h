#ifndef CRF_STDSPARSEFEATUREMAP_H_
#define CRF_STDSPARSEFEATUREMAP_H_
/*
 * CRF_StdFeatureMap.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */

#include "CRF_FeatureMap.h"
#include "CRF_StdFeatureMap.h"
/*
 * class CRF_StdSparseFeatureMap
 * subclass of CRF_StdFeatureMap
 *
 * Implements a very basic "sparse" FeatureMap for CRF processing.  This
 * map assumes that every label will be associated with every input feature.
 * The association between features and labels or transitions can be turned on
 * and off, but only at a high level (i.e. a CRF can be created with no transition
 * features, or no state features, but it cannot be created where only one label
 * has state features).
 *
 * What makes this map different than the "standard" map (that it is a subclass of)
 * is that this map allows the input feature stream to be sparse.  A sparse stream
 * has each feature defined as a pair of floats - the first element of the pair is
 * the feature id number, and the second is the actual value.
 *
 * Sparseness does not extend to assignment of features to labels, and as with the
 * "standard" this map requires that all labels receive all features in the input
 * stream.
 *
 */
class CRF_StdSparseFeatureMap : public CRF_StdFeatureMap
{
public:
	CRF_StdSparseFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas);
	CRF_StdSparseFeatureMap(CRF_FeatureMap_config* cnf);
	virtual ~CRF_StdSparseFeatureMap();

	virtual double computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab);
	virtual double computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab);
	virtual double computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab, bool compute_grad=true);
	virtual double computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad=true);
	virtual void accumulateFeatures(float *ftr_buf, double *accumulator, QNUInt32 lab);
};

#endif /*CRF_STDFEATUREMAP_H_*/
