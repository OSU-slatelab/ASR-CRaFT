#include "CRF_FeatureMap.h"
#include "CRF_StdFeatureMap.h"
#include "CRF_StdSparseFeatureMap.h"


CRF_FeatureMap::CRF_FeatureMap(QNUInt32 nlabs, QNUInt32 nfeas)
{
}

CRF_FeatureMap::CRF_FeatureMap(CRF_FeatureMap_config* cnf)
  :config(cnf)
{
	this->numFtrFuncs=0;
}


CRF_FeatureMap::~CRF_FeatureMap()
{
}

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

/*double CRF_FeatureMap::computeRi(float* ftr_buf, double* lambda, QNUInt32 &lc, QNUInt32 clab)
{
	return 0.0;
}

double CRF_FeatureMap::computeMij(float* ftr_buf, double* lambda, QNUInt32 &lc, QNUInt32 plab, QNUInt32 clab)
{
	return 0.0;
}*/

double CRF_FeatureMap::computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab)
{
	return 0.0;
}

double CRF_FeatureMap::computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab)
{
	return 0.0;
}


/*double CRF_FeatureMap::computeExpFState(float* ftr_buf, double* lambda, QNUInt32 &lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab)
{
	return 0.0;
}

double CRF_FeatureMap::computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32 &lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab)
{
	return 0.0;
}*/

double CRF_FeatureMap::computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab,bool compute_grad)
{
	return 0.0;
}

double CRF_FeatureMap::computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad)
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
	return QN_UINT32_MAX;
}

QNUInt32 CRF_FeatureMap::getStateBiasIdx(QNUInt32 clab)
{
	return QN_UINT32_MAX;
}

QNUInt32 CRF_FeatureMap::getTransFeatureIdx(QNUInt32 plab, QNUInt32 clab, QNUInt32 fno)
{
	return QN_UINT32_MAX;
}

QNUInt32 CRF_FeatureMap::getTransBiasIdx(QNUInt32 clab, QNUInt32 plab)
{
	return QN_UINT32_MAX;
}

QNUInt32 CRF_FeatureMap::getNumFtrFuncs()
{
	return this->numFtrFuncs;
}

QNUInt32 CRF_FeatureMap::getNumStates()
{
	return this->config->numStates;
}

void CRF_FeatureMap::setNumStates(QNUInt32 ns)
{
	this->config->numStates=ns;
}

QNUInt32 CRF_FeatureMap::recalc()
{
	return 0;
}

string CRF_FeatureMap::getMapDescriptor(QNUInt32 lambdaNum)
{
	return string("");
}

void CRF_FeatureMap::accumulateFeatures(float *ftr_buf,double *accumulator,QNUInt32 lab) {
	return;
}
