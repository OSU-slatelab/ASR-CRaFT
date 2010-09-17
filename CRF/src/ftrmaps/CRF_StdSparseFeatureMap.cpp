/*
 * CRF_StdSparseFeatureMap.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_StdSparseFeatureMap.h"

/*
 * CRF_StdSparseFeatureMap constructor
 *
 * Input: nlabs - number of possible labels in the CRF
 *        nfeas - number of features per label in the CRF
 *
 */
CRF_StdSparseFeatureMap::CRF_StdSparseFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas)
	: CRF_StdFeatureMap(nlabs,nfeas)
{
	this->recalc();
}

/*
 * CRF_StdSparseFeatureMap constructor
 *
 * Input: *cnf - struct containing feature mapping information
 *
 * Note:  This constructor is preferred for new code.
 *
 */
CRF_StdSparseFeatureMap::CRF_StdSparseFeatureMap(CRF_FeatureMap_config* cnf)
	: CRF_StdFeatureMap(cnf)
{
	this->recalc();
}

/*
 * CRF_StdSparseFeatureMap destructor
 *
 *
 */
CRF_StdSparseFeatureMap::~CRF_StdSparseFeatureMap()
{
}

/*
 * CRF_StdSparseFeatureMap::computeStateArrayValue
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        clab - label to compute the CRF value for
 *
 * Returns: computed state value for the label clab given the set of features (ftr_buf)
 *          and the current CRF labmda values (lambda)
 *
 * Contains logic to read the elements of ftr_buf as index/value pairs.
 */
double CRF_StdSparseFeatureMap::computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab)
{
	double stateValue=0.0;
	QNUInt32 lc=this->stateFeatureIdxCache[clab];
	//QNUInt32 tmp_lc;
	if (config->useStateFtrs) {
		for (QNUInt32 fidx=0; fidx<config->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= config->stateFidxStart) && (new_idx <= config->stateFidxEnd)) {
				//tmp_lc=lc+new_idx;
				//cout << "S: Computing for new_idx: " << new_idx << " value " << ftr_buf[fidx+1] << " tmp_lc " << tmp_lc << endl;
				stateValue+=ftr_buf[fidx+1]*lambda[lc+new_idx];
			}
		}
	}
	if (config->useStateBias) {
		//tmp_lc=lc+this->numStateFuncs-1;
		//cout << "SB: Computing for new_idx: XX" << " value " << 1 << " tmp_lc " << tmp_lc << " C: " << clab << endl;
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs-1; //Advance to end of state functions, then backup
		stateValue+=lambda[lc+this->numStateFuncs-1];
	}
	return stateValue;
}

/*
 * CRF_StdSparseFeatureMap::computeTransMatrixValue
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        plab - the previous label to compute the CRF value from
 *        clab - the current label to compute the CRF value for
 *
 * Returns: computed state value for the transition between labels plab and clab
 *           given the set of features (ftr_buf) and the current CRF labmda values (lambda)
 *
 * Contains logic to read the elements of ftr_buf as index/value pairs.
  */
double CRF_StdSparseFeatureMap::computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab)
{
	// Here, lc is the start of the trans feature weights for the combination of plab,clab
	double transMatrixValue=0.0;
	QNUInt32 lc=this->transFeatureIdxCache[plab*config->numLabs+clab];
	if (config->useTransFtrs) {
		for (QNUInt32 fidx=0; fidx<config->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= config->transFidxStart) && (new_idx<=config->transFidxEnd)) {
				//cout << "T: Computing for new_idx: " << new_idx << " value " << ftr_buf[fidx+1] << " tmp_lc " << tmp_lc << endl;
				//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+plab*this->numTransFuncs+new_idx;
				transMatrixValue+=ftr_buf[fidx+1]*lambda[lc+new_idx];
				//lc++;
			}
		}
	}
	if (config->useTransBias) {
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+(plab+1)*this->numTransFuncs-1;
		//cout << "TB: Computing for new_idx: XX"  << " value " << 1 << " tmp_lc " << tmp_lc << " P: " << plab << " C: "<< clab << endl;
		transMatrixValue+=lambda[lc+this->numTransFuncs-1];
		//lc++;
	}
	return transMatrixValue;
}

/*
 * CRF_StdSparseFeatureMap::computeStateExpF
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
 * Returns: log likelihood of the state features
 *
 * Fills the ExpF vector and gradient vector with their values using the pre-computed
 * alpha_beta (probability of being in the state with the label clab).
 *
 * Contains logic to read the elements of ftr_buf as index/value pairs.
 *
 */
double CRF_StdSparseFeatureMap::computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab,bool compute_grad)
{
	double logLi=0.0;
	QNUInt32 lc = this->stateFeatureIdxCache[clab];
	QNUInt32 tmp_lc;
	if (config->useStateFtrs) {
		for (QNUInt32 fidx=0; fidx<config->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= config->stateFidxStart) && (new_idx <= config->stateFidxEnd)) {
				tmp_lc=lc+new_idx;
				ExpF[tmp_lc]+=alpha_beta*ftr_buf[fidx+1];
				if (compute_grad && (t_clab == clab)) {
					grad[tmp_lc]+=ftr_buf[fidx+1];
					logLi += lambda[tmp_lc]*ftr_buf[fidx+1];
				}
			}
		}
	}
	if (config->useStateBias) {
		tmp_lc=lc+this->numStateFuncs-1;
		ExpF[tmp_lc]+=alpha_beta;
		if (compute_grad && (t_clab==clab)) {
			grad[tmp_lc]+=1;
			logLi += lambda[tmp_lc];
		}
	}
	return logLi;
}

/*
 * CRF_StdSparseFeatureMap::computeTransExpF
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
 * Returns: log likelihood of the state features
 *
 * Fills the ExpF vector and gradient vector with their values using the pre-computed
 * alpha_beta (probability of being in the transition between plab and clab).
 *
 * Contains logic to read the elements of ftr_buf as index/value pairs.
 */
double CRF_StdSparseFeatureMap::computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad)
{
	double logLi=0.0;
	QNUInt32 lc = this->transFeatureIdxCache[plab*config->numLabs+clab];
	QNUInt32 tmp_lc;
	if (config->useTransFtrs) {
		for (QNUInt32 fidx=0; fidx<config->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= config->transFidxStart) && (new_idx<=config->transFidxEnd)) {
				tmp_lc=lc+new_idx;
				ExpF[tmp_lc]+=alpha_beta*ftr_buf[fidx+1];
				if (compute_grad && (clab==t_clab) && (plab==t_plab)) {
					grad[tmp_lc]+=ftr_buf[fidx+1];
					logLi += lambda[tmp_lc]*ftr_buf[fidx+1];
				}
				//lc++;
			}
		}
	}
	if (config->useTransBias) {
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+(plab+1)*this->numTransFuncs-1;
		tmp_lc=lc+this->numTransFuncs-1;
		ExpF[tmp_lc] += alpha_beta;
		if (compute_grad && (clab==t_clab) && (plab==t_plab)) {
			if (compute_grad) { grad[tmp_lc]+=1;}
			logLi+=lambda[tmp_lc];
		}
	}
	return logLi;
}

/*
 * CRF_StdSparseFeatureMap::accumulateFeatures
 *
 * Fills the vector *accumulator with features from  *ftr_buf for the label clab.
 *   Used in training the CRF via AIS.
 */
void CRF_StdSparseFeatureMap::accumulateFeatures(float *ftr_buf,double *accumulator,QNUInt32 clab) {
	int offset=clab*(config->stateFidxEnd-config->stateFidxStart);

	for (QNUInt32 fidx=0; fidx<config->numFeas; fidx=fidx+2)
	{
		QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
		if ((new_idx >= config->stateFidxStart) && (new_idx <= config->stateFidxEnd)) {
			accumulator[offset+new_idx]+=ftr_buf[fidx+1];
		}
	}
}
