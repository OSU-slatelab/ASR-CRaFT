#include "CRF_FeatureMap.h"
#include "CRF_StdFeatureMap.h"
#include "CRF_StdSparseFeatureMap.h"

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

double CRF_FeatureMap::computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab)
{
	return 0.0;
}

double CRF_FeatureMap::computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab)
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

double CRF_FeatureMap::computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab)
{
	return 0.0;
}

double CRF_FeatureMap::computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab)
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

QNUInt32 CRF_FeatureMap::getStateFeatureIdx(QNUInt32 clab, QNUInt32 fno)
{
	return 0;
}

QNUInt32 CRF_FeatureMap::getTransFeatureIdx(QNUInt32 plab, QNUInt32 clab, QNUInt32 fno)
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

string CRF_FeatureMap::getMapDescriptor(QNUInt32 lambdaNum)
{
	return string("");
}
