/*
 * CRF_StdNStateNodeMasked.h
 *
 *  Created on: Jul 29, 2009
 *      Author: morrijer
 */

#ifndef CRF_STDNSTATENODEMASKED_H_
#define CRF_STDNSTATENODEMASKED_H_

#include "CRF_StdNStateNode.h"

class CRF_StdNStateNodeMasked: public CRF_StdNStateNode {
public:
	CRF_StdNStateNodeMasked(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdNStateNodeMasked();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
};

#endif /* CRF_STDNSTATENODEMASKED_H_ */
