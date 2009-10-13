#include "CRF_Seq.h"

using namespace std;

CRF_Seq::CRF_Seq(float* fb, QNUInt32 sizeof_fb, QNUInt32 l, QNUInt32 sizeof_mats, CRF_Seq* prev_node)
	: ftr_buf(fb),
	  ftr_len(sizeof_fb),
	  lab(l),
	  mat_len(sizeof_mats),
	  prev(prev_node)
{
	this->alpha_buf=new double[sizeof_mats];
	this->beta_buf=new double[sizeof_mats];
	this->alpha_beta_buf=new double[sizeof_mats];
	this->M_trans_buf=new double[sizeof_mats*sizeof_mats];
	this->M_state_buf=new double[sizeof_mats];
	for (QNUInt32 i=0; i<sizeof_mats; i++) {
		this->alpha_buf[i]=0;
		this->beta_buf[i]=0;
	}
	this->mat_len=sizeof_mats;
	this->scale=1.0;
	this->next=NULL;
}

CRF_Seq::~CRF_Seq()
{
	delete [] this->ftr_buf;
	delete [] this->alpha_buf;
	delete [] this->beta_buf;
	delete [] this->alpha_beta_buf;
	delete [] this->M_trans_buf;
	delete [] this->M_state_buf;
}

void CRF_Seq::reset(float* fb, QNUInt32 sizeof_fb, QNUInt32 l, QNUInt32 sizeof_mats)
{
	memcpy(this->ftr_buf,fb,sizeof_fb*sizeof(double));
	this->ftr_len=sizeof_fb;
	this->lab=l;
	this->mat_len=sizeof_mats;
	this->scale=1.0;
}


void CRF_Seq::setNext(CRF_Seq* next_node)
{
	this->next=next_node;
}

float* CRF_Seq::getFtr()
{
	return this->ftr_buf;
}

QNUInt32 CRF_Seq::getFtrLen()
{
	return this->ftr_len;
}

void CRF_Seq::setMtrans(double* M_in) {
	this->M_trans_buf=M_in;
}

void CRF_Seq::setMstate(double* R_in) {
	this->M_state_buf=R_in;
}

void CRF_Seq::setAlpha(double* alpha_in) {
	this->alpha_buf = alpha_in;
}

void CRF_Seq::setBeta(double* beta_in) {
	this->beta_buf = beta_in;
}

void CRF_Seq::setScale(double scale_in) {
	this->scale=scale_in;
}

double* CRF_Seq::getMtrans() {
	return this->M_trans_buf;
}

double* CRF_Seq::getMstate() {
	return this->M_state_buf;
}

double* CRF_Seq::getAlpha() {
	return this->alpha_buf;
}

double* CRF_Seq::getBeta() {
	return this->beta_buf;
}

double* CRF_Seq::getAlphaBeta() {
	return this->alpha_beta_buf;
}

double CRF_Seq::getScale() {
	return this->scale;
}

QNUInt32 CRF_Seq::getLab()
{
	return this->lab;
}

CRF_Seq* CRF_Seq::getNext()
{
	return this->next;
}

CRF_Seq* CRF_Seq::getPrev()
{
	return this->prev;
}
