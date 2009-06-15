#include "CRF_FeatureMap.h"

CRF_FeatureMap::CRF_FeatureMap()
{
	this->numStates=1;
}

CRF_FeatureMap::CRF_FeatureMap(QNUInt32 nlabs, QNUInt32 nfeas)
{
	this->numStates=1;
}


CRF_FeatureMap::~CRF_FeatureMap()
{
}

double CRF_FeatureMap::computeRi(float* ftr_buf, double* lambda, QNUInt32 &lc, QNUInt32 clab)
{
	return 0.0;
}

double CRF_FeatureMap::computeMij(float* ftr_buf, double* lambda, QNUInt32 &lc, QNUInt32 plab, QNUInt32 clab)
{
	return 0.0;
}

double CRF_FeatureMap::computeExpFState(float* ftr_buf, double* lambda, QNUInt32 &lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab)
{
	return 0.0;
}

double CRF_FeatureMap::computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32 &lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab)
{
	return 0.0;
}

QNUInt32 CRF_FeatureMap::getNumStateFuncs(QNUInt32 clab)
{
	return 0;
}

QNUInt32 CRF_FeatureMap::getNumTransFuncs(QNUInt32 plab, QNUInt32 clab)
{
	return 0;
}


QNUInt32 CRF_FeatureMap::getNumFtrFuncs()
{
	return this->numFtrFuncs;
}

QNUInt32 CRF_FeatureMap::getNumStates()
{
	return this->numStates;
}

void CRF_FeatureMap::setNumStates(QNUInt32 ns)
{
	this->numStates=ns;
}

QNUInt32 CRF_FeatureMap::recalc()
{
	return 0;
}
