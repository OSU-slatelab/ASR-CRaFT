#ifndef CRF_LATTICEBUILDER_H_
#define CRF_LATTICEBUILDER_H_

#include "fst/lib/fstlib.h"
#include "CRF.h"
#include "CRF_Model.h"
#include "CRF_FeatureStream.h"
#include "CRF_StateVector.h"
#include "CRF_StdStateVectorLog.h"
#include "CRF_StdNStateVectorLog.h"

using namespace fst;

class CRF_LatticeBuilder
{
protected:
	CRF_StateVector* nodeList;
	CRF_Model* crf;
	CRF_FeatureStream* ftr_strm;
	float* ftr_buf;
	QNUInt32* lab_buf;
	QNUInt32 bunch_size;
	QNUInt32 num_ftrs;
	QNUInt32 num_labs;
public:
	CRF_LatticeBuilder(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_LatticeBuilder();
	virtual StdVectorFst* testBuild();
	virtual StdVectorFst* buildLattice();
	virtual StdVectorFst* bestPath(bool align);
	virtual StdVectorFst* LMBestPath(bool align, StdFst* lmFst);
	virtual StdVectorFst* nStateBuildLattice();
	virtual StdVectorFst* nStateBuildLattice(StdVectorFst* labFst);
	virtual StdVectorFst* nStateBestPath(bool align);
	virtual StdVectorFst* nStateLMBestPath(bool align, StdFst* lmFst);
	virtual StdVectorFst* nStateBestPath_old(bool align);
};

#endif /*CRF_LATTICEBUILDER_H_*/
