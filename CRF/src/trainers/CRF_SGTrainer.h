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

/*
 * class CRF_SGTrainer
 *
 * Used in training.  This class implements stochastic gradient descent training for CRFs.
 */

class CRF_SGTrainer: public CRF_Trainer
{
protected:
	int logspaceTrain;
public:
	CRF_SGTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str_mgr, char* wt_fname);
	void train();
};

#endif /*CRF_SGTRAINER_H_*/
