/*
 * CRF_NewGradBuilder_StdSeg_NoDur_NoTrans.h
 *
 *  Created on: Dec 1, 2011
 *      Author: Yanzhang (Ryan) He
 */

#ifndef CRF_NEWGRADBUILDER_STDSEG_NODUR_NOTRANS_H_
#define CRF_NEWGRADBUILDER_STDSEG_NODUR_NOTRANS_H_

#include <vector>

#include "../../CRF.h"
#include "../../nodes/CRF_StdStateNode.h"
#include "../../nodes/CRF_StateVector.h"
#include "CRF_NewGradBuilder_StdSeg.h"

class CRF_NewGradBuilder_StdSeg_NoDur_NoTrans : public CRF_NewGradBuilder_StdSeg {
public:
	CRF_NewGradBuilder_StdSeg_NoDur_NoTrans(CRF_Model* crf_in);
	virtual ~CRF_NewGradBuilder_StdSeg_NoDur_NoTrans();
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
};

#endif /* CRF_NEWGRADBUILDER_STDSEG_NODUR_NOTRANS_H_ */
