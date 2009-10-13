/*
 * CRF_Pthread_GradAccumulator.h
 *
 *  Created on: Jul 7, 2009
 *      Author: fosler
 */

#ifndef CRF_PTHREAD_GRADACCUMULATOR_H_
#define CRF_PTHREAD_GRADACCUMULATOR_H_

#include "CRF_GradAccumulator.h"
#include "../gradbuilders/CRF_GradBuilder.h"
#include <pthread>

class CRF_Pthread_GradAccumulator_Thread
{
public:
	CRF_Pthread_GradAccumulator_Thread(CRF_GradBuilder *builder,double *grad);
	~CRF_Pthread_GradAccumulator_Thread();
	int start(CRF_FeatureStream *ftr_str);
	inline double getTmpZx() { return tmpZx; }
	inline double getLogLi() { return logLi; }
	inline CRF_GradBuilder* getBuilder() { return builder; }
	inline QNUInt32 getUttCount() { return uttCount; }
	int join();
protected:
	int run();
	static void * threadEntry(void*);
private:
	pthread_t threadId;
	CRF_GradBuilder *builder;
	CRF_FeatureStream *ftr_str;
	double *grad;
	double tmpZx;
	double logLi;
	QNUInt32 uttCount;
};


class CRF_Pthread_GradAccumulator: public CRF_GradAccumulator {
public:
	CRF_Pthread_GradAccumulator(CRF_Model *myCrf,bool myLogspace,int myNStates);
	virtual ~CRF_Pthread_GradAccumulator();
	double accumulateGradient(CRF_FeatureStreamManager* ftr_str_mgr, int nStreams, double* grad, QNUInt32 *uttCount);
	inline void setUttReport(int u) { uttReport=u; }
};

#endif /* CRF_PTHREAD_GRADACCUMULATOR_H_ */
