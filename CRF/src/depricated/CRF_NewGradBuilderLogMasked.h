/*
 * CRF_NewGradBuilderLogMasked.h
 *
 *  Created on: Jul 27, 2009
 *      Author: morrijer
 */

#ifndef CRF_NEWGRADBUILDERLOGMASKED_H_
#define CRF_NEWGRADBUILDERLOGMASKED_H_

#include "CRF_NewGradBuilderLog.h"

class CRF_NewGradBuilderLogMasked: public CRF_NewGradBuilderLog {
public:
	CRF_NewGradBuilderLogMasked(CRF_Model* crf_in);
	virtual ~CRF_NewGradBuilderLogMasked();
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
};

#endif /* CRF_NEWGRADBUILDERLOGMASKED_H_ */
