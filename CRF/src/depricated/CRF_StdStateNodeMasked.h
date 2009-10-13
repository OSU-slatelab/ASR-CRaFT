/*
 * CRF_StdStateNodeMasked.h
 *
 *  Created on: Jul 29, 2009
 *      Author: morrijer
 */

#ifndef CRF_STDSTATENODEMASKED_H_
#define CRF_STDSTATENODEMASKED_H_

#include "CRF_StdStateNode.h"

class CRF_StdStateNodeMasked: public CRF_StdStateNode {
public:
	CRF_StdStateNodeMasked(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdStateNodeMasked();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
};

#endif /* CRF_STDSTATENODEMASKED_H_ */
