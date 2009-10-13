#ifndef CRF_GRADBUILDER_H_
#define CRF_GRADBUILDER_H_

#include "../../CRF.h"
#include "../../CRF_Model.h"
#include "../../io/CRF_FeatureStream.h"
#include "../../nodes/CRF_StateVector.h"
//#include "CRF_Seq.h"

class CRF_GradBuilder
{
protected:
	CRF_Model* crf;
	double* alpha_base;
	double* tmp_beta;
	double* ExpF;
	float* ftr_buf;
	QNUInt32* lab_buf;
	QNUInt32 num_labs;
	QNUInt32 lambda_len;
	CRF_StateVector* nodeList;
public:
	CRF_GradBuilder(CRF_Model* crf_in);
	virtual ~CRF_GradBuilder();
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
/*	virtual double computeTransMatrix(CRF_Seq*);
	virtual double computeTransMatrixLog(CRF_Seq*);*/
	virtual void setNodeList(CRF_StateVector* nl);
	// factory method
	static CRF_GradBuilder *create(CRF_Model *crf_ptr, bool useLogspace,int nStates);
};

#endif /*CRF_GRADBUILDER_H_*/
