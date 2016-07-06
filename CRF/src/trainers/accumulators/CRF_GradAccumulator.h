/*
 * CRF_GradAccumulator.h
 *
 *  Created on: Jul 7, 2009
 *      Author: fosler
 */

#ifndef CRF_GRADACCUMULATOR_H_
#define CRF_GRADACCUMULATOR_H_

#include "../../CRF_Model.h"
#include "../../io/CRF_FeatureStreamManager.h"

class CRF_GradAccumulator {
protected:
	CRF_Model* crf;
	bool useLogspace;
	int nStates;
	int uttReport;
	objfunctype objective;

	// Added by Ryan
	QNUInt32 minibatch;

public:
	CRF_GradAccumulator(CRF_Model *myCrf,bool myLogspace,int myNStates);
	virtual ~CRF_GradAccumulator();
	virtual double accumulateGradient(CRF_FeatureStreamManager* ftr_strm, int nStreams, double* grad,QNUInt32 *uttCount);

	// Added by Ryan
	virtual double accumulateGradientMinibatch(CRF_FeatureStreamManager* ftr_strm, int nStreams, double* grad, double* Zx_out, QNUInt32 *uttCount, bool *endOfIter);

	inline void setUttReport(int u) { uttReport=u; }
	inline void setObjectiveFunction(objfunctype ofunc) {objective=ofunc;}

	// Added by Ryan
	inline void setMinibatch(QNUInt32 minibatch) { this->minibatch = minibatch; }
};

#endif /* CRF_GRADACCUMULATOR_H_ */
