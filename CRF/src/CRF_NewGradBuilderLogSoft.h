/*
 * CRF_NewGradBuilderLogSoft.h
 *
 *  Created on: Oct 6, 2009
 *      Author: morrijer
 */

#ifndef CRF_NEWGRADBUILDERLOGSOFT_H_
#define CRF_NEWGRADBUILDERLOGSOFT_H_

#include "CRF_GradBuilder.h"
#include "CRF_LatticeBuilder.h"
#include <vector>

class CRF_NewGradBuilderLogSoft: public CRF_GradBuilder {
private:
	CRF_LatticeBuilder* lb;
	vector<double> alphaAlignBase;
public:
	CRF_NewGradBuilderLogSoft(CRF_Model* crf_in);
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
	virtual ~CRF_NewGradBuilderLogSoft();
};

#endif /* CRF_NEWGRADBUILDERLOGSOFT_H_ */
