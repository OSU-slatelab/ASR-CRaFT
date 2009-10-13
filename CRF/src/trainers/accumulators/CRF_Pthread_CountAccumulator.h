/*
 * CRF_Pthread_CountAccumulator.h
 *
 *  Created on: Aug 24, 2009
 *      Author: fosler
 */


#ifndef CRF_PTHREAD_COUNTACCUMULATOR_H_
#define CRF_PTHREAD_COUNTACCUMULATOR_H_
#include "CRF_CountAccumulator.h"
#include "../../io/CRF_FeatureStreamManager.h"
#include <pthread>

class CRF_Pthread_CountAccumulator_Thread : public CRF_CountAccumulator {
public:
	CRF_Pthread_CountAccumulator_Thread(CRF_Model *crf_ptr_in,CRF_FeatureStream *ftr_strm_in, int uttRpt_in);
	virtual ~CRF_Pthread_CountAccumulator_Thread();
	int join();
	int start();
protected:

	int run();
	static void * threadEntry(void*);
private:
	pthread_t threadId;
};

class CRF_Pthread_CountAccumulator : public CRF_CountAccumulator {
protected:
	int nthreads;
	CRF_Pthread_CountAccumulator_Thread **children;
	CRF_FeatureStreamManager *ftr_str_mgr;
public:
	CRF_Pthread_CountAccumulator(CRF_Model *crf_ptr_in,CRF_FeatureStreamManager *ftr_strm_in, int nthreads, int uttRpt_in);
	virtual ~CRF_Pthread_CountAccumulator();
	virtual void accumulate();
	virtual void reset();
};

#endif /* CRF_PTHREAD_COUNTACCUMULATOR_H_ */
