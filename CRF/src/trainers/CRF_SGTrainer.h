#ifndef CRF_SGTRAINER_H_
#define CRF_SGTRAINER_H_
/*
 * CRF_SGTrainer.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Implements the CRF_SGTrainer class
 */
#include "../CRF.h"
#include "CRF_Trainer.h"
#include "gradbuilders/CRF_GradBuilder.h"
#include "../nodes/CRF_StateVector.h"
//#include "accumulators/CRF_GradAccumulator.h"
//#include "accumulators/CRF_Pthread_GradAccumulator.h"
#include "accumulators/CRF_Minibatch_GradAccumulator.h"

/*
 * class CRF_SGTrainer
 *
 * Used in training.  This class implements stochastic gradient descent training for CRFs.
 */

class CRF_SGTrainer: public CRF_Trainer
{
protected:
	int logspaceTrain;

	// Added by Ryan, for AdaGrad
	double eta; // the scaling factor for AdaGrad
	bool useAdagrad;
	// very small epsilon added to the dividend to avoid zero dividend for AdaGrad,
	// when gradient square sum = 0
	double eps;

	// Added by Ryan, for minibatch training
	int nThreads;
	int minibatch;

public:
	CRF_SGTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str_mgr, char* wt_fname);
	void train();

	// Added by Ryan
	void setNThreads(int n);
	void setMinibatch(int m);
	void setEta(double eta);
	void setUseAdagrad(double useAdagrad);

private:
	void sgtrain();
	void sgtrainMinibatch();
};

#endif /*CRF_SGTRAINER_H_*/
