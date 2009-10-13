/*
 * CRF_StdStateNodeLogMasked.h
 *
 *  Created on: Jul 27, 2009
 *      Author: morrijer
 */

#ifndef CRF_STDSTATENODELOGMASKED_H_
#define CRF_STDSTATENODELOGMASKED_H_

#include "CRF_StdStateNodeLog.h"

class CRF_StdStateNodeLogMasked: public CRF_StdStateNodeLog {
public:
	CRF_StdStateNodeLogMasked(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StdStateNodeLogMasked();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
};

#endif /* CRF_STDSTATENODELOGMASKED_H_ */
