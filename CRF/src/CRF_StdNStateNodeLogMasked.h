/*
 * CRF_StdNStateNodeLogMasked.h
 *
 *  Created on: Jul 29, 2009
 *      Author: morrijer
 */

#ifndef CRF_STDNSTATENODELOGMASKED_H_
#define CRF_STDNSTATENODELOGMASKED_H_

#include "CRF_StdNStateNodeLog.h"

class CRF_StdNStateNodeLogMasked: public CRF_StdNStateNodeLog {
public:
	CRF_StdNStateNodeLogMasked(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdNStateNodeLogMasked();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
};

#endif /* CRF_STDNSTATENODELOGMASKED_H_ */
