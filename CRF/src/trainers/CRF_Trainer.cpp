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
	// Added by Ryan, supporting auto-resume
	// get the directory of the weight files
	string weight_fname_str = weight_fname;
	unsigned found = weight_fname_str.find_last_of('/');  // assuming Unix-like file system
	if (found == string::npos) {
		this->weight_dir = ".";
	} else {
		this->weight_dir = weight_fname_str.substr(0, found);
	}

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

// Added by Ryan
/*
 * CRF_Trainer::getWeightDir()
 *
 * Accessor function to get the weight directory
 *
 */
string CRF_Trainer::getWeightDir() {
	return this->weight_dir;
}

// Added by Ryan
/*
 * CRF_Trainer::touchDoneFileIter
 *
 * Touch the .done.train.i* file after the ith iteration finishes.
 *
 */
bool CRF_Trainer::touchDoneFileIter(int iter) {
	string donefile;
	stringstream ss;
	ss << this->weight_dir << "/.done.train.i" << iter;
	ss >> donefile;
	std::ofstream outfile (donefile.c_str());
	if (outfile.is_open()) {
		outfile.close();
		return true;
	} else {
		cerr << "ERROR: cannot touch the done file " << donefile << endl;
		return false;
	}
}

// Added by Ryan
/*
 * CRF_Trainer::touchDoneFileAll
 *
 * Touch the .done.train file after the training is done.
 *
 */
bool CRF_Trainer::touchDoneFileFinal() {
	string donefile;
	stringstream ss;
	ss << this->weight_dir << "/.done.train";
	ss >> donefile;
	std::ofstream outfile (donefile.c_str());
	if (outfile.is_open()) {
		outfile.close();
		return true;
	} else {
		cerr << "ERROR: cannot touch the done file " << donefile << endl;
		return false;
	}
}

// Added by Ryan
/*
 * CRF_Trainer::setLRDecayRate
 *
 * Mutator function to set the learning rate decaying factor (for training algorithms that use a learning rate)
 */
void CRF_Trainer::setLRDecayRate(float lr_decay_rate_in) {
	if (lr_decay_rate_in < -1e-6 || lr_decay_rate_in - 1.0 > 1e-6)  {
		string errstr="CRF_Trainer::setLRDecayRate() threw exception: "
				"lr_decay_rate has to be between 0.0 and 1.0. "
				"Current value: " + stringify(lr_decay_rate_in);
		throw runtime_error(errstr);
	}
	this->lr_decay_rate = lr_decay_rate_in;
}
