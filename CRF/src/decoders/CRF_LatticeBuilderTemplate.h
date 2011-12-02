/*
 * CRF_LatticeBuilderTemplate.h
 *
 *  Created on: Oct 26, 2011
 *      Author: hey
 */

#ifndef CRF_LATTICEBUILDERTEMPLATE_H_
#define CRF_LATTICEBUILDERTEMPLATE_H_

#include "fst/fstlib.h"
#include "../CRF.h"
#include "../CRF_Model.h"
#include "../io/CRF_FeatureStream.h"
#include "../nodes/CRF_StateVector.h"

template <class LatticeBuilder>
class CRF_LatticeBuilderTemplate
{
public:
	const LatticeBuilder lb;
	CRF_LatticeBuilderTemplate(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_LatticeBuilderTemplate();

//public:
//	CRF_LatticeBuilderTemplate(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
//	virtual ~CRF_LatticeBuilderTemplate();
//	template <class Arc> int buildLattice(VectorFst<Arc>*fst, bool align=false, VectorFst<Arc>*alignFst=NULL,bool norm=true);
//	template <class Arc> int nStateBuildLattice(VectorFst<Arc>*fst, bool align=false, VectorFst<Arc>*alignFst=NULL, bool norm=true);
//	virtual StdVectorFst* buildLattice();
//	virtual int getAlignmentGammas(vector<double> *denominatorStateGamma,
//									vector<double> *numeratorStateGamma,
//									vector<double> *denominatorTransGamma,
//									vector<double> *numeratorTransGamma);
//	virtual int getAlignmentGammasNState(vector<double> *denominatorStateGamma,
//									vector<double> *numeratorStateGamma,
//									vector<double> *denominatorTransGamma,
//									vector<double> *numeratorTransGamma);
//	virtual CRF_StateVector* getNodeList();
//	virtual void computeAlignedAlphaBeta(Fst<LogArc>& fst, int nstates);
};

template <class LatticeBuilder>
CRF_LatticeBuilderTemplate<LatticeBuilder>::CRF_LatticeBuilderTemplate(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in)
	: lb(ftr_strm_in, crf_in) {
}

template <class LatticeBuilder>
CRF_LatticeBuilderTemplate<LatticeBuilder>::~CRF_LatticeBuilderTemplate() {
}

#endif /* CRF_LATTICEBUILDERTEMPLATE_H_ */
