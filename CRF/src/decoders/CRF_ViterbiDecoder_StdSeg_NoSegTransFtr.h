/*
 * CRF_ViterbiDecoder_StdSeg_NoSegTransFtr.h
 *
 *  Created on: Feb 6, 2012
 *      Author: hey
 */

#ifndef CRF_VITERBIDECODER_STDSEG_NOSEGTRANSFTR_H_
#define CRF_VITERBIDECODER_STDSEG_NOSEGTRANSFTR_H_

#include "fst/fstlib.h"
#include "../CRF.h"
#include "../CRF_Model.h"
#include "../io/CRF_FeatureStream.h"
#include "../nodes/CRF_StateVector.h"
#include "CRF_ViterbiDecoder.h"  // for the definition of class CRF_ViterbiState
#include <vector>
#include <map>
#include <deque>
//#include <hash_map>
//#include <unordered_map>

using namespace fst;
using namespace std;
typedef StdArc::StateId StateId;
typedef StdArc::Weight Weight;

const int NOT_VALID_FST_STATE_ID = -1;  // Not a valid state ID
const int BEFORE_START_TIME_STAMP = -1;  // the starting point before the first frame. This MUST be -1.

template <class CRF_VtbNode> class CRF_ViterbiDecoder_StdSeg_NoSegTransFtr;  // forward declaration

// Added by Ryan
class CRF_TimedState {
public:
	int time_stamp;  // time stamp of the utterance. -1 means starting point before the first frame.
	int state_id;    // state id of the lm lattice

	CRF_TimedState(int ts, int sid) : time_stamp(ts), state_id(sid) {}
	bool operator == (const CRF_TimedState& p) const {
		return (time_stamp == p.time_stamp && state_id == p.state_id);
	}
	bool operator < (const CRF_TimedState& p) const {
		// definite order - order by time stamp, then by state id
		if (time_stamp == p.time_stamp) {
			return (state_id < p.state_id);
		}
		else {
			return (time_stamp < p.time_stamp);
		}
	}
};

// Added by Ryan
/*
 * Segmental lattice state composed of lm/dict fst state and phone fst state
 */
class CRF_ComposedLat_Seg_State {
public:
	uint state; // the state id used in lm/dict fst
	uint phn_id; // the phone id used in phone fst
	uint dur; // duration of the segment

	CRF_ComposedLat_Seg_State(uint st, uint phn, uint d) : state(st), phn_id(phn), dur(d) {}
	bool operator == (const CRF_ComposedLat_Seg_State p) const {
		return (state==p.state && phn_id==p.phn_id && dur==p.dur);
	}

	bool operator < (const CRF_ComposedLat_Seg_State p) const {
		// definite order - order first by state, second by label, last by dur
		if (state == p.state) {
			if (phn_id == p.phn_id)
			{
				return (dur < p.dur);
			}
			else {
				return (phn_id < p.phn_id);
			}
		}
		else {
			return (state < p.state);
		}
	}
};

/*
 * class CRF_ViterbiDecoder_StdSeg_NoSegTransFtr
 *
 * Performs a best-pass Viterbi decoding over the possible CRF outputs, constrained
 * by a language model in finite-state transducer format to provide word-level
 * transcriptions.
 *
 * This decoder is specifically for segmental CRF without segmental transition
 * features (but transition bias and frame-level transition features are allowed).
 * When the maximum segment length is 1, it reduces to a frame level CRF decoder.
 *
 */

/*
 * class CRF_ViterbiNode
 *
 * The list of all viterbi states at the node of current time step. It also
 * maintains an ID map for mapping states to their index in the list, and
 * the minimum weight (the weight of the best path) at the current time step.
 *
 */
class CRF_ViterbiNode
{
protected:

	// Each element in these lists corresponds to a hypothesized phone/word state
	// in the lm/dict lattice.
	// Every state must have a non-epsilon input phone.
	// The StateIdMap is a map from state ID to the index of it in these lists.
	// The min_weight is the minimum weight of the best path to any of these hypothesized states.
	vector<uint> viterbiStateIds;
	vector<uint> viterbiPhnIds;
	vector<uint> viterbiWrdIds;
	vector<float> viterbiWts; // each weight in this list is the minimum weight of the weights of n states for the same phone/word in viterbiWts_nState.
	map<CRF_ComposedLat_Seg_State, uint> viterbiStateIdMap;
	float min_weight;		// minimum weight of any path up to the current state
	vector<bool> isPhoneStartBoundary; // if the first state of the current hypothesized phone is a starting boundary of this phone

	// Each element in these lists corresponds to one of the internal multi-states of a phone
	// for a hypothesized state above.
	// For each hypothesized phone state stored in the lists above, there must be n individual
	// ordered internal phone states in the lists below, according to the n-state phone model.
	// So, the size of the lists below are always n times of that of the lists above.
	vector<float> viterbiWts_nStates;
	vector<int> viterbiPtrs_nStates;
	vector<uint> viterbiDurs_nStates;
	uint nStates;
	vector<float> viterbiAcouWts_nStates;
	vector<float> viterbiLmWts_nStates;

	vector<uint> bestSegPointers_nStates;
	map<CRF_ComposedLatState, uint> startState_pointer_of_bestSegPointers_nStates_map;
	vector<float> bestSegWts;

	uint trans_addedCounter;
	uint trans_updateCheckCounter;
	uint trans_updateCounter;
	uint addedCounter;
	uint updateCheckCounter;
	uint updateCounter;

public:
	friend class CRF_ViterbiDecoder_StdSeg_NoSegTransFtr<CRF_ViterbiNode>;

	CRF_ViterbiNode(uint nStates_in) : nStates(nStates_in)
	{
		// infinity for our initial minimum weight
		min_weight = 99999.0;

		trans_addedCounter = 0;
		trans_updateCheckCounter = 0;
		trans_updateCounter = 0;
		addedCounter = 0;
		updateCheckCounter = 0;
		updateCounter = 0;
	}
	virtual ~CRF_ViterbiNode(){}

	void clear(){
		viterbiStateIds.clear();
		viterbiPhnIds.clear();
		viterbiWrdIds.clear();
		viterbiWts.clear();
		viterbiStateIdMap.clear();
		min_weight = 99999.0;
		isPhoneStartBoundary.clear();

		viterbiAcouWts_nStates.clear();
		viterbiLmWts_nStates.clear();

		viterbiWts_nStates.clear();
		viterbiPtrs_nStates.clear();
		viterbiDurs_nStates.clear();

		bestSegPointers_nStates.clear();
		startState_pointer_of_bestSegPointers_nStates_map.clear();
		bestSegWts.clear();

		trans_addedCounter = 0;
		trans_updateCheckCounter = 0;
		trans_updateCounter = 0;
		addedCounter = 0;
		updateCheckCounter = 0;
		updateCounter = 0;
	}

	void addNonEpsVtbState(uint stateId, uint phnId, uint wrdId, uint dur,
			float min_wt, float* wts_nStates, int* ptrs_nStates, bool isStartBound,
			float* acou_wts_nStates, float* lm_wts_nStates, double beam=0.0)
	{

		// just for debugging
//		cout << "addNonEpsVtbState(): stateId=" << stateId << " phnId=" << phnId
//				<< " wrdId=" << wrdId << " dur=" << dur << " min_wt=" << min_wt;
//		for (uint st = 0; st < nStates; st++)
//		{
//			cout << " wts_nStates[" << st << "]=" << wts_nStates[st];
//		}
//		cout << " isStartBound=" << isStartBound << endl;


//		bool prune = true;
//		if (beam <= 0.0) { prune = false; }
//
//		if (min_wt < min_weight + beam || !prune) {
//
//			if (min_wt < min_weight) {
//				min_weight = min_wt;
//			}

		// just for debugging
//		if (stateId == 2 || phnId == 2)
//		{
//			assert(stateId == 2 || phnId == 2);
//			cout << " stateId=" << stateId << ", phnId=" << phnId << ", wrdId=" << wrdId
//					<< ", dur=" << dur << ", min_wt=" << min_wt << ", wts_nStates[0]=" << wts_nStates[0]
//					<< ", ptrs_nStates[0]=" << ptrs_nStates[0] << ", isStartBound=" << isStartBound
//					<< ", acou_wts_nStates[0]=" << acou_wts_nStates[0]
//				    << ", lm_wts_nStCRF_ViterbiNodeates[0]=" << lm_wts_nStates[0]
//					<< ", beam=" << beam << endl;
//		}

			CRF_ComposedLat_Seg_State composed_lat_seg_state(stateId, phnId, dur);
			map<CRF_ComposedLat_Seg_State,uint>::iterator vtbStateIt;
			vtbStateIt = viterbiStateIdMap.find(composed_lat_seg_state);
			if (vtbStateIt == viterbiStateIdMap.end()) {
				trans_addedCounter++;

				// just for debugging
//				if (stateId == 2 || phnId == 2)
//				{
//				cout << "  addNonEpsVtbState(): addedCounter++, addedCounter=" << addedCounter << endl;
//				}

				// the above means that it is not in our expansion list
				// so we expand it as new and add it
				viterbiStateIds.push_back(stateId);
				viterbiPhnIds.push_back(phnId);
				viterbiWrdIds.push_back(wrdId);
				viterbiWts.push_back(min_wt);
				viterbiStateIdMap[composed_lat_seg_state] = viterbiStateIds.size() - 1;
				for (uint st = 0; st < nStates; st++) {
					viterbiWts_nStates.push_back(wts_nStates[st]);
					viterbiPtrs_nStates.push_back(ptrs_nStates[st]);
					viterbiDurs_nStates.push_back(dur);

					viterbiAcouWts_nStates.push_back(acou_wts_nStates[st]);
					viterbiLmWts_nStates.push_back(lm_wts_nStates[st]);

					// just for debugging
//					if (stateId == 2 || phnId == 2)
//					{
//					cout << "-->[Inserted!] nState=" << st << ", viterbiWts_nStates["
//							<< viterbiWts_nStates.size()-1 << "]=" << viterbiWts_nStates[viterbiWts_nStates.size()-1] <<
//							", viterbiDurs_nStates[" << viterbiDurs_nStates.size()-1 << "]=" <<
//							viterbiDurs_nStates[viterbiDurs_nStates.size()-1] << endl;
//					}
				}
				isPhoneStartBoundary.push_back(isStartBound);
			}
			else {
				trans_updateCheckCounter++;

				// just for debugging
//				if (stateId == 2 || phnId == 2)
//				{
//				cout << "  addNonEpsVtbState(): updateCheckCounter++, updateCheckCounter=" << updateCheckCounter << endl;
//				}

				uint stateIdxInVtbList = vtbStateIt->second;
				if (min_wt < viterbiWts[stateIdxInVtbList]) {
					trans_updateCounter++;

					// just for debugging
//					if (stateId == 2 || phnId == 2)
//					{
//					cout << "  addNonEpsVtbState(): updateCounter++, updateCounter=" << updateCounter << endl;
//					}

					// this transition is better than the one we've previously
					// expanded on, so replace the one we've got with it
//					viterbiPhnIds[stateIdxInVtbList] = phnId;  // viterbiStateIds and viterbiPhnIds remains unchanged
					viterbiWrdIds[stateIdxInVtbList] = wrdId;
					viterbiWts[stateIdxInVtbList] = min_wt;
				}
				for (uint st = 0; st < nStates; st++) {
					int nStateIdxInVtbList = stateIdxInVtbList * nStates + st;

					// just for debugging
//					if (stateId == 2 || phnId == 2)
//					{
//					cout << "  wts_nStates[" << st << "]=" << wts_nStates[st] <<
//							", viterbiWts_nStates[" << nStateIdxInVtbList << "]=" <<
//							viterbiWts_nStates[nStateIdxInVtbList] <<
//							", viterbiDurs_nStates[" << nStateIdxInVtbList << "]=" <<
//							viterbiDurs_nStates[nStateIdxInVtbList] << endl;
//					}

					if (wts_nStates[st] < viterbiWts_nStates[nStateIdxInVtbList]) {

						viterbiWts_nStates[nStateIdxInVtbList] = wts_nStates[st];
						viterbiPtrs_nStates[nStateIdxInVtbList] = ptrs_nStates[st];
						viterbiDurs_nStates[nStateIdxInVtbList] = dur;

						viterbiAcouWts_nStates[nStateIdxInVtbList] = acou_wts_nStates[st];
						viterbiLmWts_nStates[nStateIdxInVtbList] = lm_wts_nStates[st];

						if (st == 0)
						{
							isPhoneStartBoundary[stateIdxInVtbList] = isStartBound;
						}

						// just for debugging
//						if (stateId == 2 || phnId == 2)
//						{
//						cout << "-->[Updated!] nState=" << st << ", viterbiWts_nStates["
//								<< nStateIdxInVtbList << "]=" << viterbiWts_nStates[nStateIdxInVtbList] <<
//								", viterbiDurs_nStates[" << nStateIdxInVtbList << "]=" <<
//								viterbiDurs_nStates[nStateIdxInVtbList] << endl;
//						}
					}
				}
			}
//		}


			// just for debugging
//			if (stateId == 2 || phnId == 2)
//			{
//				assert(stateId == 2 || phnId == 2);
//			}
	}

	// Choose the best segment (i.e. the best duration) for each nState of each phone state in the node.
	// This function has to be run after both transition and state values are updated for all segments.
	void choose_nState_Best_Seg()
	{
		// just for debugging
//		cout << "viterbiStateIds.size()=" << viterbiStateIds.size() << endl;

		bestSegPointers_nStates.clear();
		startState_pointer_of_bestSegPointers_nStates_map.clear();
		bestSegWts.clear();

		for (uint stateIdxInVtbList = 0; stateIdxInVtbList < viterbiStateIds.size(); stateIdxInVtbList++)
		{
			uint state_id = viterbiStateIds[stateIdxInVtbList];
			uint phn_id = viterbiPhnIds[stateIdxInVtbList];

			CRF_ComposedLatState composed_lat_state(state_id, phn_id);
			map<CRF_ComposedLatState,uint>::iterator vtbStateIt;
			vtbStateIt = startState_pointer_of_bestSegPointers_nStates_map.find(composed_lat_state);

			if (vtbStateIt == startState_pointer_of_bestSegPointers_nStates_map.end()) {
				addedCounter++;

				// this means this state is not yet in the best segment list.
				// so we add it to the list as the best segment we've seen so far.
				startState_pointer_of_bestSegPointers_nStates_map[composed_lat_state] =
						bestSegPointers_nStates.size();
				uint stateIdxInVtbList_nState_startState = stateIdxInVtbList * nStates;
				for (uint st = 0; st < nStates; st++) {
					bestSegPointers_nStates.push_back(stateIdxInVtbList_nState_startState + st);
				}
			}
			else {
				updateCheckCounter++;

				// this means this state has a certain segment already saved in our best segment list.
				// we have to check if the current segment is even better than the one in the best list.

				uint startState_pointer_of_bestSegPointers_nStates = vtbStateIt->second;
				for (uint st = 0; st < nStates; st++) {
					uint nStateCurSegPointer = stateIdxInVtbList * nStates + st;
					uint nStateBestSegPointer = bestSegPointers_nStates[startState_pointer_of_bestSegPointers_nStates + st];

					if (viterbiWts_nStates[nStateCurSegPointer] < viterbiWts_nStates[nStateBestSegPointer])
					{
						// this segment is better than the one we've previously
						// saved for the best segment, so replace the one we've got with it
						updateCounter++;
						bestSegPointers_nStates[startState_pointer_of_bestSegPointers_nStates + st] = nStateCurSegPointer;
					}
				}
			}
		}

		// just for debugging
//		cout << "bestSegPointers_nStates.size()=" << bestSegPointers_nStates.size() <<
//				", startState_pointer_of_bestSegPointers_nStates_map.size()=" << startState_pointer_of_bestSegPointers_nStates_map.size() <<
//				", bestSegWts.size()=" << bestSegWts.size() << endl;
//		uint bestSegWtsCount = 0;

		for (uint startState_pointer_of_bestSegPointers_nStates = 0;
				startState_pointer_of_bestSegPointers_nStates < bestSegPointers_nStates.size();
				startState_pointer_of_bestSegPointers_nStates += nStates)
		{
			// just for debugging
//			cout << "startState_pointer_of_bestSegPointers_nStates=" << startState_pointer_of_bestSegPointers_nStates << ", ";

			uint nStateBestSegPointer = bestSegPointers_nStates[startState_pointer_of_bestSegPointers_nStates];
			float bestSegWeight = viterbiWts_nStates[nStateBestSegPointer];
			for (uint st = 1; st < nStates; st++) {
				uint nStateBestSegPointer = bestSegPointers_nStates[startState_pointer_of_bestSegPointers_nStates + st];
				if (viterbiWts_nStates[nStateBestSegPointer] < bestSegWeight)
				{
					bestSegWeight = viterbiWts_nStates[nStateBestSegPointer];
				}
			}
			bestSegWts.push_back(bestSegWeight);

			// just for debugging
//			cout << "bestSegWts[" << bestSegWtsCount++ << "]=" << bestSegWeight << endl;
		}

		// just for debugging
//		cout << "bestSegWts.size()=" << bestSegWts.size() << endl;
	}

	// Find the minimum weight
	// This function has to be run after both transition and state values are updated for all segments.
	void findMinWeight()
	{
		for (uint idx = 0; idx < viterbiStateIds.size(); idx++)
		{
			if (viterbiWts[idx] < min_weight)
			{
				min_weight = viterbiWts[idx];
			}
		}
	}
};

template <class CRF_VtbNode> class CRF_ViterbiDecoder_StdSeg_NoSegTransFtr
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
	QNUInt32 nStates;

	uint epsilonCounter;
	uint checkedCounter;
//	uint addedCounter;
//	uint updateCheckCounter;
//	uint updateCounter;
	uint epsAdded;
	uint epsModified;
	uint epsMapSize;
	uint num_pruned;
	uint int_time;
	uint cross_time;
	uint merge_time;

	// the list of next lab_max_dur nodes after the current node before pruning.
	// reason to use deque: nextViterbiNodes_unPruned would be a node window with
	// size lab_max_dur to succeeding nodes following the current. The window shifts
	// to the right by 1 node every time, which pops the head node to curViterbiNode_unPruned
	// and pushes a new node to the back. A deque would be better than a vector for
	// this purpose since it can efficiently pop its head.
//	deque<CRF_ViterbiNode*> nextViterbiNodes_unPruned;
	deque<CRF_VtbNode*> nextViterbiNodes_unPruned;  // its size should be lab_max_dur-1 instead of lab_max_dur.

	// the current node before pruning
//	CRF_ViterbiNode* curViterbiNode_unPruned;
	CRF_VtbNode* curViterbiNode_unPruned;

	// the current node after pruning
	// a simplified version of CRF_ViterbiNode.
	// each vector definition as in CRF_ViterbiNode.
	vector<uint>* curViterbiStateIds;
	vector<uint>* curViterbiWrdIds;
	vector<float>* curViterbiWts_nStates;
	vector<float>* curViterbiAcouWts_nStates;
	vector<float>* curViterbiLmWts_nStates;

	// the previous node after pruning
	// a simplified version of CRF_ViterbiNode.
	// each vector definition as in CRF_ViterbiNode.
	vector<uint>* prevViterbiStateIds;
	vector<uint>* prevViterbiWrdIds;
	vector<float>* prevViterbiWts_nStates;
	double prev_min_weight;
	vector<float>* prevViterbiAcouWts_nStates;
	vector<float>* prevViterbiLmWts_nStates;

	size_t labs_width;
	QNUInt32 lab_max_dur;
	QNUInt32 ftr_buf_size;
	QNUInt32 lab_buf_size;
	QNUInt32 nActualLabs;

	bool if_output_full_fst;
	VectorFst<StdArc>* output_full_fst;
//	__gnu_cxx::hash_map<TimedState,int>* fullfst_timedState_id_map;
//	unordered_map<TimedState,int> fullfst_timedState_id_map;
	map<CRF_TimedState,int> timedLmState_to_fullFstStateId_map;
	vector<CRF_TimedState> fullFstStateId_to_timedLmState_vector;

public:
	CRF_ViterbiDecoder_StdSeg_NoSegTransFtr(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_ViterbiDecoder_StdSeg_NoSegTransFtr();
	virtual void setIfOutputFullFst(bool ifFull);
	virtual void stateValueUpdate(uint nodeCnt);
	virtual void stateValueUpdate_onOutputFullFst(uint nodeCnt);
	virtual float internalStateTransUpdate(uint nodeCnt, uint stateId, uint phn_id, uint wrd_id, uint prevStateIdxInPrunedVtbList, double beam=0.0);
	virtual float crossStateTransUpdate(uint nodeCnt, uint state_id, uint phn_id, uint wrd_id, uint prevPhn_EndStateIdx_InPrunedVtbList, uint prev_phn_id, float expand_wt, double beam=0.0);
	virtual void expandCrossStateFromPrevNode(uint nodeCnt, VectorFst<StdArc>* lm_fst, double beam=0.0);
	virtual void pruning(uint nodeCnt, double beam=0.0);
	virtual int nStateDecode(VectorFst<StdArc>* fst, VectorFst<StdArc>* lm_fst, VectorFst<StdArc>* out_full_fst, double beam=0.0, uint min_hyps=0, uint max_hyps=0, float beam_inc=0.05);
	virtual CRF_StateVector* getNodeList();
	virtual void createFreePhoneLmFst(VectorFst<StdArc>* new_lm_fst);
	virtual StateId getOutputFullFstStateIdByTimedLmStateId(int nodeCnt, int lm_stateId);
	virtual StateId findOrInsertLmStateToOutputFullFst(int nodeCnt, int lm_stateId);
	virtual bool getTimedPhnNStateIdByOutputFullFstStateId(int fst_state_id, int& ret_nodeCnt, int& ret_lm_stateId, int& ret_phn_nState);
	virtual void insertArcToOutputFullFst(
			uint nodeCnt, uint stateId, uint phnId, uint wrdId, uint phn_nState, uint dur,
			uint prev_stateId, uint prev_phn_nState, float trans_wt, bool isCrossWordBoundary);
};

#endif /* CRF_VITERBIDECODER_STDSEG_NOSEGTRANSFTR_H_ */
