#ifndef CRF_NEWGRADBUILDER_H_
#define CRF_NEWGRADBUILDER_H_
#include <vector>

#include "../../CRF.h"
#include "CRF_GradBuilder.h"
#include "../../nodes/CRF_StdStateNode.h"
#include "../../nodes/CRF_StateVector.h"
//#include "CRF_StdStateVector.h"

class CRF_NewGradBuilder : public CRF_GradBuilder
{
protected:
	//vector <CRF_StateNode*> nodeList;

public:
	CRF_NewGradBuilder(CRF_Model* crf_in);
	virtual double buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out);
	virtual ~CRF_NewGradBuilder();
};

#endif /*CRF_NEWGRADBUILDER_H_*/
