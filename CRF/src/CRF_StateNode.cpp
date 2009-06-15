#include "CRF_StateNode.h"

CRF_StateNode::CRF_StateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf_in)
	: ftrBuf(fb),
	  ftrBuf_size(sizeof_fb),
	  label(lab),
	  crf_ptr(crf_in)
{
}

CRF_StateNode::~CRF_StateNode()
{
	delete [] this->ftrBuf;
	if (this->alphaArray != NULL) { delete this->alphaArray; }
	if (this->betaArray != NULL) { delete this->betaArray; }
	if (this->alphaBetaArray != NULL) { delete this->betaArray; }
}

double CRF_StateNode::computeTransMatrix()
{
	return 0;
}


double CRF_StateNode::computeTransMatrixLog()
{
	return 0;
}

double CRF_StateNode::computeAlpha(double* prev_alpha)
{
	return 0;
}

double CRF_StateNode::computeFirstAlpha(double* prev_alpha)
{
	return this->computeAlpha(prev_alpha);
}

double CRF_StateNode::computeBeta(double* result_beta, double scale)
{
	return 0;
}

double* CRF_StateNode::computeAlphaBeta(double Zx)
{
	return NULL;
}

void CRF_StateNode::setTailBeta()
{
}

double CRF_StateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	return 0;
}

double CRF_StateNode::computeAlphaSum()
{
	return 0;
}

void CRF_StateNode::reset(float *fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf_in)
{
	//memcpy(fb,this->ftrBuf,sizeof_fb);
	if (this->ftrBuf != NULL) { delete[] this->ftrBuf; }
	this->ftrBuf=fb;
	this->ftrBuf_size=sizeof_fb;
	this->label=lab;
	this->crf_ptr=crf_in;
	this->alphaScale=0.0;
}

double* CRF_StateNode::getAlpha()
{
	return this->alphaArray;
}

double* CRF_StateNode::getBeta()
{
	return this->betaArray;
}

double* CRF_StateNode::getAlphaBeta()
{
	return this->alphaBetaArray;
}

QNUInt32 CRF_StateNode::getLabel()
{
	return this->label;
}

double CRF_StateNode::getAlphaScale()
{
	return this->alphaScale;
}

double CRF_StateNode::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	return 0.0;
}

double CRF_StateNode::getStateValue(QNUInt32 cur_lab)
{
	return 0.0;
}

double CRF_StateNode::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	return 0.0;
}
