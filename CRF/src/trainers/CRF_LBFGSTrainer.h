#ifndef CRF_LBFGSTRAINER_H_
#define CRF_LBFGSTRAINER_H_
/*
 * CRF_LBFGSTrainer.h
 *
 * Copyright (c) 2010
 * Author: Eric Fosler-Lussier
 *
 * Implements the CRF_LBFGSTrainer class
 */

#include "../CRF.h"
#include "CRF_Trainer.h"
#include "gradbuilders/CRF_GradBuilder.h"
#include "accumulators/CRF_Pthread_GradAccumulator.h"
#include "accumulators/CRF_GradAccumulator.h"
#include "../utils/lbfgs.h"

/*
 * class CRF_LBFGSTrainer
 *
 * Used in training.  This class implements L-BFGS gradient descent training for CRFs.
 */
class CRF_LBFGSTrainer: public CRF_Trainer
{
protected:
	int logspaceTrain;
	bool start;
	double invSquareVar;
	int iCounter;
	CRF_GradBuilder* gbuild;
	CRF_GradAccumulator* gaccum;
	double *grad;

public:
	CRF_LBFGSTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str, char* wt_fname);
	void train();
	lbfgsfloatval_t evaluateGradient(const lbfgsfloatval_t *myLambda,lbfgsfloatval_t *gradient, const int lambdaLen, const lbfgsfloatval_t step);
	static lbfgsfloatval_t _evaluateGradient(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,const lbfgsfloatval_t step);
	int progress(const lbfgsfloatval_t *lambda, const lbfgsfloatval_t *gradient, const lbfgsfloatval_t logLikelihood, const lbfgsfloatval_t lambdaNorm, const lbfgsfloatval_t gradientNorm, const lbfgsfloatval_t step, int lambdaLen, int lbfgsIter, int numberOfEvals);
	static int _progress(void *instance,const lbfgsfloatval_t *x,const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);
};

#endif /*CRF_LBFGSTRAINER_H_*/
