/*
 * CRF_FerrGradBuilder.h
 *
 *  Created on: Oct 2, 2009
 *      Author: morrijer
 */

#ifndef CRF_FERRGRADBUILDER_H_
#define CRF_FERRGRADBUILDER_H_

#include "CRF_GradBuilder.h"
#include "../../decoders/CRF_LatticeBuilder.h"

class CRF_FerrGradBuilder: public CRF_GradBuilder {
private:
	CRF_LatticeBuilder* lb;
public:
	CRF_FerrGradBuilder(CRF_Model* crf_in);
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
	virtual ~CRF_FerrGradBuilder();
};

#endif /* CRF_FERRGRADBUILDER_H_ */
