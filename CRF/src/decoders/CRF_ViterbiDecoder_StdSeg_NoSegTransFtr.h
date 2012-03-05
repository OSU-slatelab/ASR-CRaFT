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

using namespace fst;
using namespace std;
typedef StdArc::StateId StateId;
typedef StdArc::Weight Weight;

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
public:

	// Each element in these lists corresponds to a hypothesized phone/word state
	// in the lm/dict lattice.
	// Every state must have a non-epsilon input phone.
	// The StateIdMap is a map from state ID to the index of it in these lists.
	// The min_weight is the minimum weight of the best path to any of these hypothesized states.
	vector<uint> viterbiStateIds;
	vector<uint> viterbiPhnIds;
	vector<uint> viterbiWrdIds;
	vector<float> viterbiWts; // each weight in this list is the minimum weight of the weights of n states for the same phone/word in viterbiWts_nState.
	map<CRF_ComposedLatState, uint> viterbiStateIdMap;
	float min_weight;
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

	uint addedCounter;
	uint updateCheckCounter;
	uint updateCounter;

	CRF_ViterbiNode(uint nStates_in) : nStates(nStates_in)
	{
		// infinity for our initial minimum weight
		min_weight = 99999.0;

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

			CRF_ComposedLatState composed_lat_state(stateId, phnId);
			map<CRF_ComposedLatState,uint>::iterator vtbStateIt;
			vtbStateIt = viterbiStateIdMap.find(composed_lat_state);
			if (vtbStateIt == viterbiStateIdMap.end()) {
				addedCounter++;

				// just for debugging
//				cout << "addNonEpsVtbState(): addedCounter++, addedCounter=" << addedCounter << endl;

				// the above means that it is not in our expansion list
				// so we expand it as new and add it
				viterbiStateIds.push_back(stateId);
				viterbiPhnIds.push_back(phnId);
				viterbiWrdIds.push_back(wrdId);
				viterbiWts.push_back(min_wt);
				viterbiStateIdMap[composed_lat_state] = viterbiStateIds.size() - 1;
				for (uint st = 0; st < nStates; st++) {
					viterbiWts_nStates.push_back(wts_nStates[st]);
					viterbiPtrs_nStates.push_back(ptrs_nStates[st]);
					viterbiDurs_nStates.push_back(dur);

					viterbiAcouWts_nStates.push_back(acou_wts_nStates[st]);
					viterbiLmWts_nStates.push_back(lm_wts_nStates[st]);
				}
				isPhoneStartBoundary.push_back(isStartBound);
			}
			else {
				updateCheckCounter++;

				// just for debugging
//				cout << "addNonEpsVtbState(): updateCheckCounter++, updateCheckCounter=" << updateCheckCounter << endl;

				uint stateIdxInVtbList = vtbStateIt->second;
				if (min_wt < viterbiWts[stateIdxInVtbList]) {
					updateCounter++;

					// just for debugging
//					cout << "addNonEpsVtbState(): updateCounter++, updateCounter=" << updateCounter << endl;

					// this transition is better than the one we've previously
					// expanded on, so replace the one we've got with it
//					viterbiPhnIds[stateIdxInVtbList] = phnId;  // viterbiStateIds and viterbiPhnIds remains unchanged
					viterbiWrdIds[stateIdxInVtbList] = wrdId;
					viterbiWts[stateIdxInVtbList] = min_wt;
				}
				for (uint st = 0; st < nStates; st++) {
					int nStateIdxInVtbList = stateIdxInVtbList * nStates + st;
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
					}
				}
			}
//		}
	}
};

class CRF_ViterbiDecoder_StdSeg_NoSegTransFtr
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
	deque<CRF_ViterbiNode*> nextViterbiNodes_unPruned;

	// the current node before pruning
	CRF_ViterbiNode* curViterbiNode_unPruned;

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

public:
	CRF_ViterbiDecoder_StdSeg_NoSegTransFtr(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_ViterbiDecoder_StdSeg_NoSegTransFtr();
	virtual void stateFeatureUpdate(uint nodeCnt);
	virtual float internalStateTransUpdate(uint nodeCnt, uint stateId, uint phn_id, uint wrd_id, uint prevStateIdxInPrunedVtbList, double beam=0.0);
	virtual float crossStateTransUpdate(uint nodeCnt, uint state_id, uint phn_id, uint wrd_id, uint prevPhn_EndStateIdx_InPrunedVtbList, uint prev_phn_id, float expand_wt, double beam=0.0);
	virtual void expandCrossStateFromPrevNode(uint nodeCnt, VectorFst<StdArc>* lm_fst, double beam=0.0);
	virtual void pruning(uint nodeCnt, double beam=0.0);
	virtual int nStateDecode(VectorFst<StdArc>* fst, VectorFst<StdArc>* lm_fst, double beam=0.0, uint min_hyps=0, uint max_hyps=0, float beam_inc=0.05);
	virtual CRF_StateVector* getNodeList();
};

#endif /* CRF_VITERBIDECODER_STDSEG_NOSEGTRANSFTR_H_ */
