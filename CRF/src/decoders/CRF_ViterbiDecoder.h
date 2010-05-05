#ifndef CRF_VITERBIDECODER_H_
#define CRF_VITERBIDECODER_H_

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


class CRF_ViterbiState {
public:
	uint label; // This is the phone label into the state
	uint wrd_label; // This is the word label into the state (if it exists)
	uint state; // This is the state_id number from the fst
	float weight; // This is the pruning weight for the state (state-internal weights stored separately)
	uint end_idx;
	CRF_ViterbiState(uint st, uint lab, float wt, uint idx) :
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
	//vector<CRF_ViterbiState> tmpViterbiStateIds;
	vector<uint> tmpViterbiStateIds_new;
	vector<uint> tmpViterbiPhnIds;
	vector<uint> tmpViterbiWrdIds;
	vector<uint> prevViterbiWrdIds;

	map<uint,uint> tmpViterbiStateIdMap_new;
	vector<float> tmpViterbiWts;
	vector<int> tmpViterbiPtrs;

	vector<float> tmpPruneWts;
	vector<float>* curViterbiWts;
	vector<float>* prevViterbiWts;
	//vector<CRF_ViterbiState>* curViterbiStateIds;
	//vector<CRF_ViterbiState>* prevViterbiStateIds;
	vector<uint>* curViterbiStateIds_new;
	vector<uint>* prevViterbiStateIds_new;
	uint updateCounter;
	uint epsilonCounter;
	uint addedCounter;
	uint checkedCounter;
	uint updateCheckCounter;
	uint dupCounter;

public:
	CRF_ViterbiDecoder(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_ViterbiDecoder();
	//virtual float internalStateUpdate(int nodeCnt, StdArc prevState, int prevIdx, vector<float>* prevWts,vector<float>* curWts, vector<int>* curPtrs);
	virtual float internalStateUpdate(int nodeCnt, uint phn_id, int prevIdx, vector<float>* prevWts,vector<float>* curWts, vector<int>* curPtrs);
	//float crossStateUpdate(int nodeCnt, const StdArc* prev_arc, int prev_idx,float base_weight, VectorFst<StdArc>* lm_fst, float min_weight, float beam);
	//float crossStateUpdate(int nodeCnt, CRF_ViterbiState prev_state, uint prev_phn, uint prev_idx,float base_weight, VectorFst<StdArc>* lm_fst, float min_weight, float beam);
	//float crossStateUpdate(int nodeCnt, uint prev_state, uint prev_phn, uint prev_idx,float base_weight, VectorFst<StdArc>* lm_fst, float min_weight, float beam);
	void crossStateUpdate_new(uint check_state, uint end_idx, float transw, uint phn_lab, uint wrd_lab);

	virtual int nStateDecode(VectorFst<StdArc>* fst, VectorFst<StdArc>* lm_fst, double beam=0.0, uint min_hyps=0, uint max_hyps=0, float beam_inc=0.05);
	virtual CRF_StateVector* getNodeList();
};

#endif /*CRF_VITERBIDECODER_H_*/
