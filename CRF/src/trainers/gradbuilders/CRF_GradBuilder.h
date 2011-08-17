#ifndef CRF_GRADBUILDER_H_
#define CRF_GRADBUILDER_H_
/*
 * CRF_GradBuilder.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Implements the CRF_GradBuilder class
 */

#include "../../CRF.h"
#include "../../CRF_Model.h"
#include "../../io/CRF_FeatureStream.h"
#include "../../nodes/CRF_StateVector.h"

/*
 * class CRF_GradBuilder
 *
 * Used to construct a gradient for gradient-based training methods.  This class is an interface class
 * and should be used to construct subclasses to build the gradient as needed by a model topology.  See
 * NewGradBuilder for one example.
 */

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
	virtual void setNodeList(CRF_StateVector* nl);
	// factory method
//	static CRF_GradBuilder *create(CRF_Model *crf_ptr, bool useLogspace,int nStates);
	static CRF_GradBuilder* create(CRF_Model *crf_ptr, objfunctype ofunc);
};

#endif /*CRF_GRADBUILDER_H_*/
