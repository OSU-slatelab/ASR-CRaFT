/*
 * CRF_ViterbiNode_PruneTrans.h
 *
 *  Created on: Mar 6, 2012
 *      Author: hey
 */

#ifndef CRF_VITERBINODE_PRUNETRANS_H_
#define CRF_VITERBINODE_PRUNETRANS_H_

#include "fst/fstlib.h"
#include "../CRF.h"
#include "../CRF_Model.h"
#include "../io/CRF_FeatureStream.h"
#include "../nodes/CRF_StateVector.h"
#include "CRF_ViterbiDecoder.h"  // for the definition of class CRF_ViterbiState
#include "CRF_ViterbiDecoder_StdSeg_NoSegTransFtr.h" // for the definition of class CRF_ViterbiDecoder_StdSeg_NoSegTransFtr
#include <vector>
#include <map>
#include <deque>

using namespace fst;
using namespace std;
typedef StdArc::StateId StateId;
typedef StdArc::Weight Weight;

/*
 * class CRF_ViterbiNode_PruneTrans
 *
 * The list of all viterbi states at the node of current time step. It also
 * maintains an ID map for mapping states to their index in the list, and
 * the minimum weight (the weight of the best path) at the current time step.
 *
 */
class CRF_ViterbiNode_PruneTrans
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
	map<CRF_ComposedLatState, uint> viterbiStateIdMap;
	float min_weight;		// minimum weight of any path up to the current state
	float min_trans_weight;	// minimum weight of any path up to the current state excluding current state feature values
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

public:
	friend class CRF_ViterbiDecoder_StdSeg_NoSegTransFtr<CRF_ViterbiNode_PruneTrans>;

	CRF_ViterbiNode_PruneTrans(uint nStates_in) : nStates(nStates_in)
	{
		// infinity for our initial minimum weight
		min_weight = 99999.0;
		min_trans_weight = 99999.0;

		addedCounter = 0;
		updateCheckCounter = 0;
		updateCounter = 0;
	}
	virtual ~CRF_ViterbiNode_PruneTrans(){}

	void clear(){
		viterbiStateIds.clear();
		viterbiPhnIds.clear();
		viterbiWrdIds.clear();
		viterbiWts.clear();
		viterbiStateIdMap.clear();
		min_weight = 99999.0;
		min_trans_weight = 99999.0;
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


		bool prune = true;
		if (beam <= 0.0) { prune = false; }

//		if (min_wt < min_weight + beam || !prune) {
		if (min_wt < min_trans_weight + beam || !prune) {

			// just for debugging
//			cout << "incoming min_wt = " << min_wt << ", min_trans_weight = "
//					<< min_trans_weight << ", beam size = "
//					<< beam << ", keep the hypothesis." << endl;

//			if (min_wt < min_weight) {
//				min_weight = min_wt;
//			}

			if (min_wt < min_trans_weight) {
				min_trans_weight = min_wt;
			}

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
		}
		else {
			// just for debugging
//			cout << "incoming min_wt = " << min_wt << ", min_trans_weight = "
//					<< min_trans_weight << ", beam size = "
//					<< beam << ", prune away the hypothesis." << endl;
		}
	}
};

#endif /* CRF_VITERBINODE_PRUNETRANS_H_ */
