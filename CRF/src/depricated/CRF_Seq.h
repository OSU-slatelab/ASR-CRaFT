#ifndef CRF_SEQ_H_
#define CRF_SEQ_H_
#include "CRF.h"

class CRF_Seq
{
protected:
	float* ftr_buf;
	QNUInt32 ftr_len;
	double* alpha_buf;
	double* beta_buf;
	double* alpha_beta_buf;
	double* M_trans_buf;
	double* M_state_buf;
	QNUInt32 lab;
	QNUInt32 mat_len;
	CRF_Seq* next;
	CRF_Seq* prev;
	float scale;
public:
	CRF_Seq(float* fb, QNUInt32 sizeof_fb, QNUInt32 l, QNUInt32 sizeof_mats, CRF_Seq* prev_node);
	virtual ~CRF_Seq();
	void reset(float* fb, QNUInt32 sizeof_fb, QNUInt32 l, QNUInt32 sizeof_mats);
	void setNext(CRF_Seq*);
	float* getFtr();
	QNUInt32 getFtrLen();
	void setMtrans(double*);
	void setMstate(double*);
	void setAlpha(double*);
	void setBeta(double*);
	void setScale(double);
	double* getMtrans();
	double* getMstate();
	double* getAlpha();
	double* getBeta();
	double* getAlphaBeta();
	double getScale();
	QNUInt32 getLab();
	CRF_Seq* getNext();
	CRF_Seq* getPrev();
};

#endif /*CRF_SEQ_H_*/
