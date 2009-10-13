#include "CRF_GradBuilder.h"
#include "CRF_NewGradBuilder.h"

/*#include "CRF_NewGradBuilderLog.h"
#include "CRF_StdStateVector.h"
#include "CRF_StdNStateVector.h"
#include "CRF_StdStateVectorLog.h"
#include "CRF_StdNStateVectorLog.h"*/

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
	cerr << "called GradBuilder destructor" << endl;
	delete[] alpha_base;
	delete[] tmp_beta;
	if (this->ftr_buf != NULL) { delete[] ftr_buf;}
	if (this->lab_buf != NULL) { delete[] lab_buf;}
	if (this->nodeList != NULL) {
		//this->nodeList->deleteAll();
		delete nodeList;
	}
}

double CRF_GradBuilder::buildGradient(CRF_FeatureStream* ftr_strm,double* grad, double* Zx_out)
{
	return 0.0;
}

/*double CRF_GradBuilder::computeTransMatrix(CRF_Seq* seq)
{
	return 0.0;
}

double CRF_GradBuilder::computeTransMatrixLog(CRF_Seq* seq)
{
	return 1.0;
}*/

void CRF_GradBuilder::setNodeList(CRF_StateVector* nl)
{
	if (this->nodeList != NULL) {
		delete nodeList;
	}
	this->nodeList = nl;
}

CRF_GradBuilder *CRF_GradBuilder::create(CRF_Model *crf_ptr,bool useLogspace, int nStates) {
	CRF_GradBuilder *gbuild;

/*	if (useLogspace) {
		if (nStates==1) {
			gbuild=new CRF_NewGradBuilderLog(crf_ptr);
			gbuild->setNodeList(new CRF_StdStateVectorLog());
			//cout << "Using Logspace training..." << endl;
		}
		else	 {
			gbuild=new CRF_NewGradBuilderLog(crf_ptr);
			gbuild->setNodeList(new CRF_StdNStateVectorLog());
		}
	}
	else{
		if (nStates==1) {
			gbuild=new CRF_NewGradBuilder(crf_ptr);
			gbuild->setNodeList(new CRF_StdStateVector());
			//cout << "Using normal space training..." << endl;
		}
		else {
			gbuild=new CRF_NewGradBuilder(crf_ptr);
			gbuild->setNodeList(new CRF_StdNStateVector());
		}
	}*/

	return new CRF_NewGradBuilder(crf_ptr);

}
