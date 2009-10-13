#ifndef CRF_SGTRAINER_H_
#define CRF_SGTRAINER_H_

#include "../CRF.h"
#include "CRF_Trainer.h"

#include "gradbuilders/CRF_NewGradBuilder.h"
#include "gradbuilders/CRF_FerrGradBuilder.h"
#include "../nodes/CRF_StateVector.h"
//#include "CRF_NewGradBuilderLog.h"
//#include "CRF_StdStateVector.h"
//#include "CRF_StdStateVectorLog.h"
//#include "CRF_StdNStateVector.h"
//#include "CRF_StdNStateVectorLog.h"
//#include "CRF_StdGradBuilder.h"
//#include "CRF_NstateGradBuilder.h"
//#include "CRF_StdGradBuilderLog.h"

class CRF_SGTrainer: public CRF_Trainer
{
protected:
	int logspaceTrain;
public:
	CRF_SGTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str_mgr, char* wt_fname);
	void train();
};

#endif /*CRF_SGTRAINER_H_*/
