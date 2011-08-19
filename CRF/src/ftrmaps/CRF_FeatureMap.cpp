/*
 * CRF_FeatureMap.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_FeatureMap.h"
#include "CRF_StdFeatureMap.h"
#include "CRF_StdSparseFeatureMap.h"

/*
 * CRF_FeatureMap constructor
 *
 * Input: nlabs - number of possible labels in the CRF
 *        nfeas - number of features per label in the CRF
 *
 */
CRF_FeatureMap::CRF_FeatureMap(QNUInt32 nlabs, QNUInt32 nfeas)
{
}

/*
 * CRF_FeatureMap constructor
 *
 * Input: *cnf - struct containing feature mapping information
 *
 * Note:  This constructor is preferred for new code.
 *
 */
CRF_FeatureMap::CRF_FeatureMap(CRF_FeatureMap_config* cnf)
  :config(cnf)
{
	this->numFtrFuncs=0;
}

/*
 * CRF_FeatureMap destructor
 *
 *
 */
CRF_FeatureMap::~CRF_FeatureMap()
{
}

/*
 * CRF_FeatureMap::createFeatureMap
 *
 * Input: *cnf - struct containing feature mapping information
 *
 * Static factory function that takes a FeatureMap configuration object as input
 * and returns a new FeatureMap object.
 * Update this function whenever a new FeatureMap object is derived.
 */
CRF_FeatureMap* CRF_FeatureMap::createFeatureMap(CRF_FeatureMap_config* cnf)
{
	CRF_FeatureMap* new_map;
	switch (cnf->map_type)
	{
	case STDSPARSE:
		new_map = new CRF_StdSparseFeatureMap(cnf);
		break;
	case STDSPARSETRANS:
		new_map = new CRF_StdSparseFeatureMap(cnf);
		break;
	case STDTRANS:
		new_map = new CRF_StdFeatureMap(cnf);
		break;
	default:
		// STDSTATE is the default
		new_map = new CRF_StdFeatureMap(cnf);
	}
	return new_map;
}

/*
 * CRF_FeatureMap::computeStateArrayValue
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        clab - label to compute the CRF value for
 *
 * Stub function.  When implemented in a subclass, this function should take the
 * features in ftr_buf and return the state value for the label clab given the
 * CRF lambda values in lambda.
 */
double CRF_FeatureMap::computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab)
{
	return 0.0;
}

/*
 * CRF_FeatureMap::computeTransMatrixValue
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        plab - the previous label to compute the CRF value from
 *        clab - the current label to compute the CRF value for
 *
 * Stub function.  When implemented in a subclass, this function should take the
 * features in ftr_buf and return the transition value for the label pair
 * plab -> clab given the CRF lambda values in lambda.
 */
double CRF_FeatureMap::computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab)
{
	return 0.0;
}

/*
 * CRF_FeatureMap::computeStateExpFValue
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        *ExpF - vector to store expected values for each state feature
 *        *grad - vector used to store gradient values
 *        alpha_beta - gamma value of the current label (clab)
 *        t_clab - true label (used in training)
 *        clab - label to compute Expected state feature values for
 *        compute_grad - flag to control whether grad vector is updated or not
 *
 * Stub function.  When implemented in a subclass, this function should take the
 * features in ftr_buf and fill the appropriate values in the ExpF vector and the
 * gradient vector using the value computed in alpha_beta.  Gradient vector is updated
 * when the true label (t_clab) is equal to the label being examined (clab).  ExpF
 * vector is updated regardless.
 */
double CRF_FeatureMap::computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab,bool compute_grad)
{
	return 0.0;
}

/*
 * CRF_FeatureMap::computeTransExpFValue
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        *ExpF - vector to store expected values for each state feature
 *        *grad - vector used to store gradient values
 *        alpha_beta - gamma value of the current transition (plab->clab)
 *        t_plab - true previous label (used in training)
 *        t_clab - true label (used in training)
 *        plab - previous label used for computing expected transition values
 *        clab - label to compute expected transition values
 *        compute_grad - flag to control whether grad vector is updated or not
 *
 * Stub function.  As computeStateExpF above, except that the expected value and
 * gradient are computed for transition pairs instead of for single labels.
 */
double CRF_FeatureMap::computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad)
{
	return 0.0;
}

/*
 * CRF_FeatureMap::getNumStateFuncs
 *
 * Input: clab - label under examination
 *
 * Stub function.  When implemented in a subclass, should return the number of
 * state functions defined for the label identified by the parameter clab.
 */
QNUInt32 CRF_FeatureMap::getNumStateFuncs(QNUInt32 clab)
{
	return 0;
}

/*
 * CRF_FeatureMap::getNumTransFuncs
 *
 * Input: plab - previous label under examination
 *        clab - label under examination
 *
 * Stub function.  When implemented in a subclass, should return the number of
 * transition functions defined for the transition pair plab->clab.
 */
QNUInt32 CRF_FeatureMap::getNumTransFuncs(QNUInt32 plab, QNUInt32 clab)
{
	return 0;
}

/*
 * CRF_FeatureMap::getStateFeatureIdx
 *
 * Input: clab - label under examination
 *        fno - feature number under examination
 *
 * Stub function.  When implemented in a subclass, should return the index into
 * the lambda vector for the input feature identified by fno and the label
 * identified by clab.
 */
QNUInt32 CRF_FeatureMap::getStateFeatureIdx(QNUInt32 clab, QNUInt32 fno)
{
	return QN_UINT32_MAX;
}

/*
 * CRF_FeatureMap::getStateBiasIdx
 *
 * Input: clab - label under examination
 *
 * Stub function.  When implemented in a subclass, should return the index into
 * the lambda vector for the state bias feature for the label identified by clab.
 */
QNUInt32 CRF_FeatureMap::getStateBiasIdx(QNUInt32 clab)
{
	return QN_UINT32_MAX;
}

/*
 * CRF_FeatureMap::getTransFeatureIdx
 *
 * Input: plab - previous label under examination
 *        clab - current label under examination
 *        fno - feature number under examination
 *
 * Stub function.  When implemented in a subclass, should return the index into
 * the lambda vector for the transition tuple of the previous label (plab), the
 * current label (clab) and the input feature indexed as fno.
 */
QNUInt32 CRF_FeatureMap::getTransFeatureIdx(QNUInt32 plab, QNUInt32 clab, QNUInt32 fno)
{
	return QN_UINT32_MAX;
}

/*
 * CRF_FeatureMap::getTransBiasIdx
 *
 * Input: plab - previous label under examination
 *        clab - current label under examination
 *        fno - feature number under examination
 *
 * Stub function.  When implemented in a subclass, should return the index into
 * the lambda vector for the transition tuple of the previous label (plab), the
 * current label (clab) and the input feature indexed as fno.
 */
QNUInt32 CRF_FeatureMap::getTransBiasIdx(QNUInt32 clab, QNUInt32 plab)
{
	return QN_UINT32_MAX;
}

/*
 * CRF_FeatureMap::getNumFtrFuncs
 *
 * Accessor function to return the total number of feature functions in the
 * feature map.
 */
QNUInt32 CRF_FeatureMap::getNumFtrFuncs()
{
	return this->numFtrFuncs;
}

/*
 * CRF_FeatureMap::getNumStates
 *
 * Accessor function to return the state topology of the crf.
 */
QNUInt32 CRF_FeatureMap::getNumStates()
{
	return this->config->numStates;
}

/*
 * CRF_FeatureMap::setNumStates
 *
 * Mutator function to set the state topology of the crf
 */

void CRF_FeatureMap::setNumStates(QNUInt32 ns)
{
	this->config->numStates=ns;
}

/*
 * CRF_FeatureMap::recalc
 *
 * Stub function.  Used to recompute internal values after using a mutator
 * function.
 */
QNUInt32 CRF_FeatureMap::recalc()
{
	return 0;
}

/*
 * CRF_FeatureMap::getMapDescriptor
 *
 * Stub function.  When implemented, should return a string that describes the
 * feature map for the input lambda index lambdaNum.
 */
string CRF_FeatureMap::getMapDescriptor(QNUInt32 lambdaNum)
{
	return string("");
}

/*
 * CRF_FeatureMap::accumulateFeatures
 *
 * Stub function.  When implemented, fill the vector *accumulator with features from
 * *ftr_buf for the label lab.  Used in training the CRF via AIS.
 */
void CRF_FeatureMap::accumulateFeatures(float *ftr_buf,double *accumulator,QNUInt32 lab) {
	return;
}

/*
 * Added by Ryan
 *
 * CRF_FeatureMap::tieGradient
 *
 * Stub function. Tie gradients for tied parameters.
 *
 */
//void CRF_FeatureMap::tieGradient(double* grad, QNUInt32 maxDur, QNUInt32 durFtrStart)
void CRF_FeatureMap::tieGradient(double* grad)
{
	return;
}

/*
 * Added by Ryan
 *
 * CRF_FeatureMap::tieGradForSingleParam
 *
 * Stub function.
 * The subroutine for tying gradients of a series of parameters to a single parameter.
 *
 */
void CRF_FeatureMap::tieGradForSingleParam(double* grad, QNUInt32 numParam, QNUInt32 start, QNUInt32 step)
{
	return;
}
