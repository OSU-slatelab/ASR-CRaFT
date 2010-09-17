#ifndef CRF_NEWGRADBUILDER_H_
#define CRF_NEWGRADBUILDER_H_
/*
 * CRF_NewGradBuilder.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Implements the CRF_NewGradBuilder class
 */
#include <vector>

#include "../../CRF.h"
#include "CRF_GradBuilder.h"
#include "../../nodes/CRF_StdStateNode.h"
#include "../../nodes/CRF_StateVector.h"

/*
 * class CRF_NewGradBuilder
 *
 * Used to construct a gradient for gradient-based training methods.  This is a basic gradient builder
 * that is suitable for a linear chain model with either a single state or multi-state toplogy.
 */

class CRF_NewGradBuilder : public CRF_GradBuilder
{

public:
	CRF_NewGradBuilder(CRF_Model* crf_in);
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
	virtual ~CRF_NewGradBuilder();
};

#endif /*CRF_NEWGRADBUILDER_H_*/
