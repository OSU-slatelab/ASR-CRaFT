/*
 * CRF_Trainer.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_Trainer.h"

/*
 * CRF_Trainer constructor
 *
 * Inputs: crf_in - CRF model
 *         ftr_str - feature stream manager for featuers and labels
 *         wt_fname - file name of file for storing weights
 */
CRF_Trainer::CRF_Trainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str, char* wt_fname)
	: crf_ptr(crf_in),
	  ftr_strm_mgr(ftr_str),
	  weight_fname(wt_fname)
{
	this->uttRpt = 100;
	this->useLogspace=0;
	this->gvar=0.0;
	this->useGvar=false;
	this->useLabelMask=false;
	this->objective=EXPF;
}

/*
 * CRF_Trainer destructor
 */
CRF_Trainer::~CRF_Trainer()
{
}

/*
 * CRF_Trainer::train
 *
 * Stub function.
 * Performs the actual training algorithm.
 */
void CRF_Trainer::train()
{
}

/*
 * CRF_Trainer::display_ftrstrm
 *
 * Displays the feature stream.
 * Used in debugging.
 */
void CRF_Trainer::display_ftrstrm()
{
	this->ftr_strm_mgr->trn_stream->display();
}

/*
 * CRF_Trainer:setMaxIters
 *
 * Mutator function to set the maximum number of training iterations
 */
void CRF_Trainer::setMaxIters(int mi) {
	this->maxIters=mi;
}

/*
 * CRF_Trainer::setLR
 *
 * Mutator function to set the learning rate (for training algorithms that use a learning rate)
 */
void CRF_Trainer::setLR(float lr_in) {
	this->lr = lr_in;
}

/*
 * CRF_Trainer::setUttRpt
 *
 * Mutator function to set the number of segments to processes before issuing a logging report
 */
void CRF_Trainer::setUttRpt(QNUInt32 rpt_in) {
	this->uttRpt = rpt_in;
}

/*
 * CRF_Trainer::setLogSpace
 *
 * Muatator function to toggle the computation space from log to real
 *
 * No longer used - computation in logspace by default
 */
void CRF_Trainer::setLogSpace(int lsf) {
	this->useLogspace=lsf;
}

void CRF_Trainer::setGaussVar(float gvar_in) {
	this->gvar=gvar_in;
	if (this->gvar != 0.0) {
		this->useGvar=true;
	}
}

void CRF_Trainer::setLabelMask(bool useMask) {
	this->useLabelMask=useMask;
}

/*
 * CRF_Trainer::setObjectiveFunction
 *
 * Mutator function to set the objective function for training
 */
void CRF_Trainer::setObjectiveFunction(objfunctype ofunc) {
	this->objective=ofunc;
}

