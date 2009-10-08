#ifndef CRF_STATENODE_H_
#define CRF_STATENODE_H_

#include <string>

#include "CRF.h"
#include "CRF_Model.h"
#include <vector>

class CRF_StateNode
{
protected:
	float* ftrBuf;
	QNUInt32 ftrBuf_size;
	QNUInt32 label;
	CRF_Model* crf_ptr;
	double* alphaArray;
	double* betaArray;
	double* alphaBetaArray;
	vector<double> alphaArrayAligned;
	vector<double> betaArrayAligned;
	QNUInt32 alphaSize;
	double alphaScale;
public:
	CRF_StateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual ~CRF_StateNode();
	virtual double computeTransMatrix();
	virtual double computeTransMatrixLog();
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual double* computeAlphaBeta(double Zx);
	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
	virtual double computeAlphaSum();
	virtual double computeAlphaAlignedSum();
	virtual void reset(float *fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf_in);
	virtual double* getAlpha();
	virtual double* getBeta();
	virtual double* getAlphaBeta();
	virtual vector<double>* getAlphaAligned();
	virtual vector<double>* getBetaAligned();
	virtual QNUInt32 getLabel();
	virtual double getAlphaScale();
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	virtual double getStateValue(QNUInt32 cur_lab);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	static CRF_StateNode* createStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf);
	virtual float *getFtrBuffer();
	virtual QNUInt32 getFtrBufferSize();
};

#endif /*CRF_STATENODE_H_*/
