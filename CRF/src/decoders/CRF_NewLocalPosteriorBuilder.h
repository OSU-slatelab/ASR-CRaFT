#ifndef CRF_NEWLOCALPOSTERIORBUILDER_H_
#define CRF_NEWLOCALPOSTERIORBUILDER_H_

/*
 * CRF_NewLocalPosteriorBuilder.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */

#include "../CRF_Model.h"
#include "../nodes/CRF_StateVector.h"
#include "../io/CRF_FeatureStream.h"

/*
 * class CRF_NewLocalPosteriorBuilder
 *
 * Used to build local posterior probability values per frame using the
 * forward-backward algorithm.  Takes as input a CRF_Model and a boolean
 * flag to control whether the posteriors returned are normalized using the Z(x)
 * normalization constant or not.
 *
 * References
 * -----------
 * - E. Fosler-Lussier and J. Morris.
 *       "Crandem Systems: Conditional Random Field Acoustic Models for
 *          Hidden Markov Models",
 *       IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP),
 *        Las Vegas, Nevada, 2008.
 * - J. Morris and E. Fosler-Lussier.
 *       "CRANDEM: Conditional Random Fields for Word Recognition",
 *       10th Conference of the International Speech Communication Association (Interspeech/ISCA),
 *        Brighton, UK, 2009.
 */

class CRF_NewLocalPosteriorBuilder
{
protected:
	CRF_Model* crf;
	CRF_StateVector* nodeList;
	bool ownsNodeList;
	float* ftr_buf;
	QNUInt32* lab_buf;
	double* alpha_base;
	bool normalize;
public:
	CRF_NewLocalPosteriorBuilder(CRF_Model* crf_in, bool norm=true);
	virtual ~CRF_NewLocalPosteriorBuilder();
	virtual CRF_StateVector* buildFtrSeq(CRF_FeatureStream* ftr_strm);
	//virtual void computeAlphaBeta(CRF_StateVector *nodeList);
	virtual CRF_StateVector* buildFtrSeqNState(CRF_FeatureStream* ftr_strm);
};

#endif /*CRF_NEWLOCALPOSTERIORBUILDER_H_*/
