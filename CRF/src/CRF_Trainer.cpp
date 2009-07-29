#include "CRF_Trainer.h"

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
}

CRF_Trainer::~CRF_Trainer()
{
}

void CRF_Trainer::train()
{
}

void CRF_Trainer::display_ftrstrm()
{
	this->ftr_strm_mgr->trn_stream->display();
}

void CRF_Trainer::setMaxIters(int mi) {
	this->maxIters=mi;
}

void CRF_Trainer::setLR(float lr_in) {
	this->lr = lr_in;
}

void CRF_Trainer::setUttRpt(QNUInt32 rpt_in) {
	this->uttRpt = rpt_in;
}

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
