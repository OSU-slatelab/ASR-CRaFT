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
public:
	CRF_GradAccumulator(CRF_Model *myCrf,bool myLogspace,int myNStates);
	virtual ~CRF_GradAccumulator();
	virtual double accumulateGradient(CRF_FeatureStreamManager* ftr_strm, int nStreams, double* grad,QNUInt32 *uttCount);
	inline void setUttReport(int u) { uttReport=u; }
	inline void setObjectiveFunction(objfunctype ofunc) {objective=ofunc;}
};

#endif /* CRF_GRADACCUMULATOR_H_ */
