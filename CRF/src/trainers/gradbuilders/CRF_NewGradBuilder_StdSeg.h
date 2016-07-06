/*
 * CRF_NewGradBuilder_StdSeg.h
 *
 *  Created on: Nov 11, 2011
 *      Author: Yanzhang (Ryan) He
 */

#ifndef CRF_NEWGRADBUILDER_STDSEG_H_
#define CRF_NEWGRADBUILDER_STDSEG_H_

#include <vector>

#include "../../CRF.h"
#include "../../nodes/CRF_StdStateNode.h"
#include "../../nodes/CRF_StateVector.h"
#include "CRF_NewGradBuilder.h"

class CRF_NewGradBuilder_StdSeg : public CRF_NewGradBuilder {
public:
	CRF_NewGradBuilder_StdSeg(CRF_Model* crf_in);
	virtual ~CRF_NewGradBuilder_StdSeg();
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
};

#endif /* CRF_NEWGRADBUILDER_STDSEG_H_ */
