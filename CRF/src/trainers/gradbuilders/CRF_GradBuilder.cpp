/*
 * CRF_GradBuilder.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_GradBuilder.h"
#include "CRF_NewGradBuilder.h"
#include "CRF_NewGradBuilderSoft.h"
//#include "CRF_FerrGradBuilder.h"

/*
 * CRF_GradBuilder constructor
 *
 * Inputs: crf_in - CRF model used for gradient construction
 */

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

/*
 * CRF_GradBuilder destructor
 */
CRF_GradBuilder::~CRF_GradBuilder()
{
	cerr << "GradBuilder destructor" << endl;
	delete[] alpha_base;
	delete[] tmp_beta;
	if (this->ftr_buf != NULL) { delete[] ftr_buf;}
	if (this->lab_buf != NULL) { delete[] lab_buf;}
	if (this->nodeList != NULL) {
		//this->nodeList->deleteAll();
		delete nodeList;
	}
	//cerr << "exiting GradBuilder destructor" << endl;
}

/*
 * CRF_GradBuilder::buildGradient
 *
 * Input: *ftr_stream - input feature stream
 *        *grad - gradient vector return value
 *        *Zx_out - normalization constant return value
 *
 * Stub function.
 * Computes the gradient given the current CRF model and the features in ftr_strm and returns it in
 * grad.  Zx_out contains the normalization constant for the current sequence.
 */
double CRF_GradBuilder::buildGradient(CRF_FeatureStream* ftr_strm,double* grad, double* Zx_out)
{
	return 0.0;
}

/*
 * CRF_GradBuilder::setNodeList
 *
 * Mutator function to set the nodelist
 */

void CRF_GradBuilder::setNodeList(CRF_StateVector* nl)
{
	if (this->nodeList != NULL) {
		delete nodeList;
	}
	this->nodeList = nl;
}

/*
CRF_GradBuilder *CRF_GradBuilder::create(CRF_Model *crf_ptr,bool useLogspace, int nStates) {
	CRF_GradBuilder *gbuild;

	return new CRF_NewGradBuilder(crf_ptr);

}
*/

/*
 * CRF_GradBuilder::create
 *
 * Input: crf_ptr - pointer to CRF model
 *        ofunc - objective function type
 *
 * Factory function to create a gradient builder object automatically.
 */
CRF_GradBuilder* CRF_GradBuilder::create(CRF_Model *crf_ptr, objfunctype ofunc)
{
	CRF_GradBuilder* gbuild;
	switch (ofunc) {
	case EXPF :
		gbuild = new CRF_NewGradBuilder(crf_ptr);
		break;
	//case EXPFSOFT :
	//	gbuild = new CRF_NewGradBuilderSoft(crf_ptr);
	//	break;
	//case FERR :
	//	gbuild = new CRF_FerrGradBuilder(crf_ptr);
	//	break;
	default :
		gbuild = new CRF_NewGradBuilder(crf_ptr);
		break;
	}
	return gbuild;
}
