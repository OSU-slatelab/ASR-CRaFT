/*
 * CRF_Minibatch_GradAccumulator.h
 *
 *  Created on: Nov 22, 2014
 *      Author: Ryan He
 */

#ifndef CRF_MINIBATCH_GRADACCUMULATOR_H_
#define CRF_MINIBATCH_GRADACCUMULATOR_H_

#include "../../CRF_Model.h"
#include "../../io/CRF_FeatureStreamManager.h"
#include "../gradbuilders/CRF_GradBuilder.h"
#include "../../ftrmaps/CRF_FeatureMap.h" // for QN_UINT32_MAX -- move elsewhere?
#include <pthread.h>

class CRF_Minibatch_GradAccumulator_Thread
{
public:
	CRF_Minibatch_GradAccumulator_Thread(CRF_GradBuilder *builder, double *grad);
	~CRF_Minibatch_GradAccumulator_Thread();
	int start(CRF_FeatureStream *ftr_str);
	int join();
	inline double getBunchLogLiDenom() { return this->bunchLogLiDenom; }
	inline double getBunchLogLiNumer() { return this->bunchLogLiNumer; }
	inline double getBunchLogLi() { return this->bunchLogLi; }
	inline CRF_GradBuilder* getBuilder() { return this->gBuilder; }
	inline QNUInt32 getUttCount() { return this->uttCount; }
	inline void setSegID(QN_SegID id) { this->segid = id; }
	inline QN_SegID getSegID() { return this->segid; }
	inline void setUttReport(int u) { this->uttReport = u; }
	inline void setMinibatchPerThread(QNUInt32 minibatchPerThread) { this->minibatchPerThread = minibatchPerThread; }

protected:
	int run();
	static void * threadEntry(void*);
private:
	pthread_t threadId;
	CRF_GradBuilder *gBuilder;
	CRF_FeatureStream *ftr_str;
	double *grad;
	double bunchLogLiDenom;
	double bunchLogLiNumer;
	double bunchLogLi;
	QNUInt32 uttCount;
	QN_SegID segid;
	int uttReport;
	QNUInt32 minibatchPerThread;
};

class CRF_Minibatch_GradAccumulator {
protected:
	CRF_Model* crf;
	CRF_FeatureStream **ftrStrms;
	QN_SegID *strmsSegids;
	int uttReport;
	objfunctype objective;
	QNUInt32 minibatch;
	CRF_GradBuilder **gBuilders;
	const QNUInt32 nStreams;

public:
	CRF_Minibatch_GradAccumulator(CRF_Model *myCrf, CRF_FeatureStreamManager* myFtrStrmMgr, QNUInt32 myNStreams);
	virtual ~CRF_Minibatch_GradAccumulator();
	//virtual double accumulateGradient(CRF_FeatureStreamManager* ftr_strm, int nStreams, double* grad,QNUInt32 *uttCount);
	virtual double accumulateGradient(
			double* grad, double* Zx_out, QNUInt32 *uttCount, bool *isEndOfIter);
	void setUttReport(int u) { this->uttReport = u; }
	void setObjectiveFunction(objfunctype ofunc) { this->objective = ofunc; }
	void setMinibatch(QNUInt32 minibatch);
	QNUInt32 getNStreams() { return this->nStreams; }
	void rewindAllAndNextSegs();
};

#endif /* CRF_MINIBATCH_GRADACCUMULATOR_H_ */
