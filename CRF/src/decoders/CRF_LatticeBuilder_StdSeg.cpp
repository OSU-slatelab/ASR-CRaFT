/*
 * CRF_LatticeBuilder_StdSeg.cpp
 *
 *  Created on: Oct 16, 2011
 *      Author: hey
 */

#include "CRF_LatticeBuilder_StdSeg.h"
CRF_LatticeBuilder_StdSeg::CRF_LatticeBuilder_StdSeg(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in)
	: CRF_LatticeBuilder(ftr_strm_in, crf_in)
{

	// just for debugging
	//cout << "beginning of CRF_LatticeBuilder_StdSeg constructor." << endl;

//	this->nodeList= new CRF_StateVector();

	// bunch_size can be only equal to 1 for
	// CRF_InFtrStream_SeqMultiWindow and CRF_InLabStream_SeqMultiWindow.
	this->bunch_size=1;

//	this->num_ftrs=this->ftr_strm->num_ftrs();
//	this->num_labs=this->crf->getNLabs();

	// Changed by Ryan
//#ifndef SEGMENTAL_CRF
//	this->ftr_buf = new float[num_ftrs*bunch_size];
//	this->lab_buf = new QNUInt32[bunch_size];
//#else
	delete [] this->ftr_buf;
	delete [] this->lab_buf;

	this->lab_max_dur = this->crf->getLabMaxDur();
	this->ftr_buf_size = this->num_ftrs * this->lab_max_dur;
	this->ftr_buf = new float[ftr_buf_size];
	this->labs_width = this->ftr_strm->num_labs();
	//this->lab_buf_size = labs_width * lab_max_dur;
	this->lab_buf_size = labs_width;
	if (this->labs_width == 0)
	{
		this->lab_buf = NULL;
	} else {
		this->lab_buf = new QNUInt32[lab_buf_size];
	}
	this->nActualLabs = this->crf->getNActualLabs();
	this->nodeStartStates = new vector<QNUInt32>();
//#endif

//	this->alpha_base = new double[this->num_labs];
//	for (QNUInt32 i=0; i<this->num_labs; i++) {
//		this->alpha_base[i]=0.0;
//	}

	// just for debugging
	//cout << "end of CRF_LatticeBuilder_StdSeg constructor." << endl;

}

CRF_LatticeBuilder_StdSeg::~CRF_LatticeBuilder_StdSeg() {
	delete this->nodeStartStates;
}

void CRF_LatticeBuilder_StdSeg::setNodeStartState(QNUInt32 nodeCnt, QNUInt32 nodeStartState)
{
	// to make sure nodeCnt is in sequential order
	assert(nodeCnt <= this->nodeStartStates->size());

	if (nodeCnt < this->nodeStartStates->size()) {
		this->nodeStartStates->at(nodeCnt) = nodeStartState;
	} else {
		this->nodeStartStates->push_back(nodeStartState);
	}
}

