#ifndef CRF_GRADBUILDER_H_
#define CRF_GRADBUILDER_H_

#include "CRF.h"
#include "CRF_Model.h"
#include "CRF_FeatureStream.h"
#include "CRF_Seq.h"
#include "CRF_StateVector.h"

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
	virtual double computeTransMatrix(CRF_Seq*);
	virtual double computeTransMatrixLog(CRF_Seq*);
	virtual void setNodeList(CRF_StateVector* nl);
};

#endif /*CRF_GRADBUILDER_H_*/
