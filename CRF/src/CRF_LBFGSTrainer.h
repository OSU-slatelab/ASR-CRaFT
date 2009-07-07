#ifndef CRF_LBFGSTRAINER_H_
#define CRF_LBFGSTRAINER_H_

#include "CRF.h"
#include "CRF_Trainer.h"

#include "CRF_NewGradBuilder.h"
#include "CRF_NewGradBuilderLog.h"
#include "CRF_StdStateVector.h"
#include "CRF_StdStateVectorLog.h"
#include "CRF_StdNStateVector.h"
#include "CRF_StdNStateVectorLog.h"
//#include "CRF_StdGradBuilder.h"
//#include "CRF_NstateGradBuilder.h"
//#include "CRF_StdGradBuilderLog.h"

#include "lbfgs.h"

class CRF_LBFGSTrainer: public CRF_Trainer
{
protected:
	int logspaceTrain;
	bool start;
	double invSquareVar;
	int iCounter;
	CRF_GradBuilder* gbuild;
	double *grad;

public:
	CRF_LBFGSTrainer(CRF_Model* crf_in, CRF_FeatureStream* ftr_str, char* wt_fname);
	void train();
	lbfgsfloatval_t evaluateGradient(const lbfgsfloatval_t *myLambda,lbfgsfloatval_t *gradient, const int lambdaLen, const lbfgsfloatval_t step);
	static lbfgsfloatval_t _evaluateGradient(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,const lbfgsfloatval_t step);
	int progress(const lbfgsfloatval_t *lambda, const lbfgsfloatval_t *gradient, const lbfgsfloatval_t logLikelihood, const lbfgsfloatval_t lambdaNorm, const lbfgsfloatval_t gradientNorm, const lbfgsfloatval_t step, int lambdaLen, int lbfgsIter, int numberOfEvals);
	static int _progress(void *instance,const lbfgsfloatval_t *x,const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);
};

#endif /*CRF_LBFGSTRAINER_H_*/
