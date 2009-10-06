#ifndef CRF_AISTRAINER_H_
#define CRF_AISTRAINER_H_

#include "CRF.h"
#include "CRF_Trainer.h"
#include "CRF_FeatureStream.h"
#include "CRF_FeatureStreamManager.h"
#include "CRF_StateVector.h"
#include "CRF_LatticeBuilder.h"
#include "CRF_FeatureMap.h" // for QN_UINT32_MAX -- move elsewhere?
#include "CRF_NewLocalPosteriorBuilder.h"
#include "CRF_CountAccumulator.h"
#include "CRF_Pthread_CountAccumulator.h"

class CRF_AISTrainer: public CRF_Trainer
{
protected:
	int logspaceTrain;
	bool start;
	double invSquareVar;
	int iCounter;
	double *grad;
	float l1alpha;

public:
	CRF_AISTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str, char* wt_fname);
	void train();
	inline void setl1alpha(float a) {l1alpha=a;}
	inline float getl1apha() {return l1alpha;}
};

#endif /*CRF_AISTRAINER_H_*/
