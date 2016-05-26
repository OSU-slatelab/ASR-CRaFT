/*
 * CRF_NewGradBuilder_StdSeg_BrokenClass.h
 *
 *  Created on: Nov 11, 2011
 *      Author: Yanzhang (Ryan) He
 */

#ifndef CRF_NEWGRADBUILDER_STDSEG_BROKENCLASS_H_
#define CRF_NEWGRADBUILDER_STDSEG_BROKENCLASS_H_

#include <vector>

#include "../../CRF.h"
#include "../../nodes/CRF_StdStateNode.h"
#include "../../nodes/CRF_StateVector.h"
#include "CRF_NewGradBuilder_StdSeg.h"

class CRF_NewGradBuilder_StdSeg_BrokenClass : public CRF_NewGradBuilder_StdSeg {
public:
	CRF_NewGradBuilder_StdSeg_BrokenClass(CRF_Model* crf_in);
	virtual ~CRF_NewGradBuilder_StdSeg_BrokenClass();
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
};

#endif /* CRF_NEWGRADBUILDER_STDSEG_BROKENCLASS_H_ */
