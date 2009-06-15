#ifndef CRF_TRAINER_H_
#define CRF_TRAINER_H_

#include "CRF.h"
#include "CRF_Model.h"
#include "CRF_FeatureStream.h"

class CRF_Trainer
{
protected:
	CRF_Model* crf_ptr;
	CRF_FeatureStream* ftr_strm;
	int maxIters;
	char* weight_fname;
	float lr;
	QNUInt32 uttRpt;
	int useLogspace;
	float gvar;
	bool useGvar;
public:
	CRF_Trainer(CRF_Model* crf_in, CRF_FeatureStream* ftr_str, char* wt_fname);
	virtual ~CRF_Trainer();
	void display_ftrstrm();
	virtual void train();
	virtual void setMaxIters(int);
	virtual void setLR(float lr_in);
	virtual void setUttRpt(QNUInt32 rpt_in);
	virtual void setLogSpace(int);
	virtual void setGaussVar(float gvar_in);
};

#endif /*CRF_TRAINER_H_*/
