#ifndef CRF_AISTRAINER_H_
#define CRF_AISTRAINER_H_
/*
 * CRF_SGTrainer.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Implements the CRF_AISTrainer class
 */

#include "../CRF.h"
#include "CRF_Trainer.h"
#include "../io/CRF_FeatureStream.h"
#include "../io/CRF_FeatureStreamManager.h"
#include "../nodes/CRF_StateVector.h"
#include "../decoders/CRF_LatticeBuilder.h"
#include "../ftrmaps/CRF_FeatureMap.h" // for QN_UINT32_MAX -- move elsewhere?
#include "../decoders/CRF_NewLocalPosteriorBuilder.h"
#include "accumulators/CRF_CountAccumulator.h"
#include "accumulators/CRF_Pthread_CountAccumulator.h"

/*
 * class CRF_AISTrainer
 *
 * Used in training.  This class implements AIS training for CRFs.
 */
class CRF_AISTrainer: public CRF_Trainer
{
protected:
	int logspaceTrain;
	bool start;
	double invSquareVar;

	// commented by Ryan. TODO: When is iCounter initialized?
	// How can we set the initial starting training iteration?
	// In other trainers, they are set by this->iCounter=this->crf_ptr->getInitIter();
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
