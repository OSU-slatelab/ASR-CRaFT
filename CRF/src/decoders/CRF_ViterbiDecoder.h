#ifndef CRF_VITERBIDECODER_H_
#define CRF_VITERBIDECODER_H_

/*
 * CRF_ViterbiDecoder.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_ViterbiState and CRF_ViterbiDecoder
 */

#include "fst/fstlib.h"
#include "../CRF.h"
#include "../CRF_Model.h"
#include "../io/CRF_FeatureStream.h"
#include "../nodes/CRF_StateVector.h"
#include <vector>
#include <map>

using namespace fst;
using namespace std;
typedef StdArc::StateId StateId;
typedef StdArc::Weight Weight;

/*
 * class CRF_ViterbiState
 *
 * Used by the CRF_ViterbiDecoder object to perform a best pass Viterbi decoding
 * for word recognition.
 *
 */
class CRF_ViterbiState {
public:
	uint label; // This is the phone label into the state
	uint wrd_label; // This is the word label into the state (if it exists)
	uint state; // This is the state_id number from the fst
	float weight; // This is the pruning weight for the state (state-internal weights stored separately)

	// Changed by Ryan. change uint to int, so that end_idx can be a null pointer -1
//	uint end_idx;
	int end_idx;

	// Changed by Ryan
	CRF_ViterbiState(uint st, uint lab, float wt, int idx) :
		label(lab), state(st), weight(wt), end_idx(idx) { wrd_label=0;}

	bool operator == (const CRF_ViterbiState p) const {
		return (state==p.state && label==p.label);
	}

	bool operator < (const CRF_ViterbiState p) const {
		// definite order - order by state, then by label
		// order is independent of weight and two ViterbiStates with different weights
		// can be considered the same state
		if (state == p.state) {
			return (label < p.label);
		}
		else {
			return (state<p.state);
		}
	}
};

// Added by Ryan
class CRF_ComposedLatState {
public:
	uint state; // the state id used in the lm/dict fst
	uint phn_id; // the phone id used in phone fst

	CRF_ComposedLatState(uint st, uint phn) : state(st), phn_id(phn) {}
	bool operator == (const CRF_ComposedLatState p) const {
		return (state==p.state && phn_id==p.phn_id);
	}

	bool operator < (const CRF_ComposedLatState p) const {
		// definite order - order by state, then by label
		if (state == p.state) {
			return (phn_id < p.phn_id);
		}
		else {
			return (state < p.state);
		}
	}
};

/*
 * class CRF_ViterbiDecoder
 *
 * Performs a best-pass Viterbi decoding over the possible CRF outputs, constrained
 * by a language model in finite-state transducer format to provide word-level
 * transcriptions.
 *
 */

class CRF_ViterbiDecoder
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
	double* alpha_base;
	vector<uint> tmpViterbiStateIds_new;
	vector<uint> tmpViterbiPhnIds;
	vector<uint> tmpViterbiWrdIds;
	vector<uint> prevViterbiWrdIds;

	// Added by Ryan
	// if the first state of the current hypothesized phone
	// is a starting boundary of this phone
	vector<bool> tmpIsPhoneStartBoundary;

	// Changed by Ryan
//	map<uint,uint> tmpViterbiStateIdMap_new;
	map<CRF_ComposedLatState,uint> tmpViterbiStateIdMap_new;


	vector<float> tmpViterbiWts;
	vector<int> tmpViterbiPtrs;

	vector<float> tmpPruneWts;
	vector<float>* curViterbiWts;
	vector<float>* prevViterbiWts;
	vector<uint>* curViterbiStateIds_new;
	vector<uint>* prevViterbiStateIds_new;
	uint updateCounter;
	uint epsilonCounter;
	uint addedCounter;
	uint checkedCounter;
	uint updateCheckCounter;
	uint dupCounter;

	// Added by Ryan
	uint crossUpdateSamePhoneCounter;
	uint crossUpdateDiffPhoneCounter;
	uint crossUpdateSameWordCounter;
	uint crossUpdateDiffWordCounter;
	uint crossUpdateSameStateCounter;
	uint crossUpdateDiffStateCounter;

public:
	CRF_ViterbiDecoder(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_ViterbiDecoder();
	virtual float internalStateUpdate(int nodeCnt, uint phn_id, int prevIdx, vector<float>* prevWts,vector<float>* curWts, vector<int>* curPtrs);
	void crossStateUpdate_new(uint check_state, uint end_idx, float transw, uint phn_lab, uint wrd_lab);

	virtual int nStateDecode(VectorFst<StdArc>* fst, VectorFst<StdArc>* lm_fst, double beam=0.0, uint min_hyps=0, uint max_hyps=0, float beam_inc=0.05);
	virtual CRF_StateVector* getNodeList();

	// added by Ryan
	virtual void createFreePhoneLmFst(VectorFst<StdArc>* new_lm_fst);
};

#endif /*CRF_VITERBIDECODER_H_*/
