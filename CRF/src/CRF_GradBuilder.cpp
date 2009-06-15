#include "CRF_GradBuilder.h"

CRF_GradBuilder::CRF_GradBuilder(CRF_Model* crf_in)
	: crf(crf_in)
{
	this->num_labs=crf_in->getNLabs();
	this->alpha_base=new double[this->num_labs];
	this->tmp_beta=new double[this->num_labs];
	this->ftr_buf=NULL;
	this->lab_buf=NULL;
	this->nodeList=NULL;
}

CRF_GradBuilder::~CRF_GradBuilder()
{
	delete[] alpha_base;
	delete[] tmp_beta;
	if (this->ftr_buf != NULL) { delete[] ftr_buf;}
	if (this->lab_buf != NULL) { delete[] lab_buf;}
	if (this->nodeList != NULL) { delete nodeList;}
}

double CRF_GradBuilder::buildGradient(CRF_FeatureStream* ftr_strm,double* grad, double* Zx_out)
{
	return 0.0;
}

double CRF_GradBuilder::computeTransMatrix(CRF_Seq* seq)
{
	return 0.0;
}

double CRF_GradBuilder::computeTransMatrixLog(CRF_Seq* seq)
{
	return 1.0;
}

void CRF_GradBuilder::setNodeList(CRF_StateVector* nl)
{
	if (this->nodeList != NULL) {
		delete nodeList;
	}
	this->nodeList = nl;
}
