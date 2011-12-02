/*
 * CRF_NewGradBuilder_StdSeg_WithoutDurLab_WithoutTransFtr.h
 *
 *  Created on: Dec 1, 2011
 *      Author: hey
 */

#ifndef CRF_NEWGRADBUILDER_STDSEG_WITHOUTDURLAB_WITHOUTTRANSFTR_H_
#define CRF_NEWGRADBUILDER_STDSEG_WITHOUTDURLAB_WITHOUTTRANSFTR_H_

#include <vector>

#include "../../CRF.h"
#include "../../nodes/CRF_StdStateNode.h"
#include "../../nodes/CRF_StateVector.h"
#include "CRF_NewGradBuilder_StdSeg.h"

class CRF_NewGradBuilder_StdSeg_WithoutDurLab_WithoutTransFtr : public CRF_NewGradBuilder_StdSeg {
public:
	CRF_NewGradBuilder_StdSeg_WithoutDurLab_WithoutTransFtr(CRF_Model* crf_in);
	virtual ~CRF_NewGradBuilder_StdSeg_WithoutDurLab_WithoutTransFtr();
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
};

#endif /* CRF_NEWGRADBUILDER_STDSEG_WITHOUTDURLAB_WITHOUTTRANSFTR_H_ */
