/*
 * CRF_ViterbiDecoder_StdSeg_NoSegTransFtr.cpp
 *
 *  Created on: Feb 6, 2012
 *      Author: hey
 */

#include "CRF_ViterbiDecoder_StdSeg_NoSegTransFtr.h"
#include <cmath>
#include <vector>
#include <deque>
#include <map>

CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::CRF_ViterbiDecoder_StdSeg_NoSegTransFtr(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in)
	: crf(crf_in),
	  ftr_strm(ftr_strm_in)
{
	this->nodeList= new CRF_StateVector();

	// bunch_size (number of windows ending at next frame) for CRF_InFtrStream_SeqMultiWindow,
	// starts from 1, and is added by 1 in each iteration, until being equal to lab_max_dur.
	this->bunch_size = 1;

	this->num_ftrs = this->ftr_strm->num_ftrs();
	this->num_labs = this->crf->getNLabs();

	this->lab_max_dur = this->crf->getLabMaxDur();
	this->ftr_buf_size = this->num_ftrs * this->lab_max_dur;
	this->ftr_buf = new float[ftr_buf_size];
	this->labs_width = this->ftr_strm->num_labs();
	//this->lab_buf_size = labs_width * lab_max_dur;
	this->lab_buf_size = this->labs_width;
	if (this->labs_width == 0)
	{
		this->lab_buf = NULL;
	} else {
		this->lab_buf = new QNUInt32[lab_buf_size];
	}
	this->nActualLabs = this->crf->getNActualLabs();

	this->alpha_base = new double[this->num_labs];
	for (QNUInt32 i=0; i<this->num_labs; i++) {
		this->alpha_base[i]=0.0;
	}

	this->nStates = this->crf->getFeatureMap()->getNumStates();

	for (QNUInt32 i = 0; i < this->lab_max_dur; i++) {
		this->nextViterbiNodes_unPruned.push_back(new CRF_ViterbiNode(this->nStates));
	}
	this->curViterbiNode_unPruned = new CRF_ViterbiNode(this->nStates);

	this->curViterbiStateIds = new vector<uint>();
	this->curViterbiWrdIds = new vector<uint>();
	this->curViterbiWts_nStates = new vector<float>();
	this->curViterbiAcouWts_nStates = new vector<float>();
	this->curViterbiLmWts_nStates = new vector<float>();

	this->prevViterbiStateIds = new vector<uint>();
	this->prevViterbiWrdIds = new vector<uint>();
	this->prevViterbiWts_nStates = new vector<float>();
	this->prevViterbiAcouWts_nStates = new vector<float>();
	this->prevViterbiLmWts_nStates = new vector<float>();

}

CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::~CRF_ViterbiDecoder_StdSeg_NoSegTransFtr()
{
	delete [] this->ftr_buf;
	delete [] this->lab_buf;
	delete [] this->alpha_base;
	delete this->nodeList;
	for (QNUInt32 i = 0; i < this->lab_max_dur; i++) {
		delete this->nextViterbiNodes_unPruned[i];
	}
	delete this->curViterbiNode_unPruned;

	delete this->curViterbiStateIds;
	delete this->curViterbiWrdIds;
	delete this->curViterbiWts_nStates;
	delete this->curViterbiAcouWts_nStates;
	delete this->curViterbiLmWts_nStates;

	delete this->prevViterbiStateIds;
	delete this->prevViterbiWrdIds;
	delete this->prevViterbiWts_nStates;
	delete this->prevViterbiAcouWts_nStates;
	delete this->prevViterbiLmWts_nStates;
}

/*
 * CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::getNodeList
 *
 * Accessor function for nodeList created during decode, if further processing
 * is desired.
 *
 */
CRF_StateVector * CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::getNodeList() {
	return this->nodeList;
}

void CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::stateFeatureUpdate(uint nodeCnt)
{
	// just for debugging
//	cout << "State feature update" << endl;

	for (uint lat_state_vtb_idx = 0;
			lat_state_vtb_idx < this->curViterbiNode_unPruned->viterbiStateIds.size();
			lat_state_vtb_idx++)
	{
		int phn_lab = this->curViterbiNode_unPruned->viterbiPhnIds[lat_state_vtb_idx] - 1;
		float min_wt = 99999.0;
		for (uint st = 0; st < this->nStates; st++)
		{
			int phn_nState_vtb_idx = lat_state_vtb_idx * this->nStates + st;
			if (this->curViterbiNode_unPruned->viterbiWts_nStates[phn_nState_vtb_idx] < 99999.0)
			{
				uint cur_phn_nState_lab = phn_lab * this->nStates + st;
				uint dur = this->curViterbiNode_unPruned->viterbiDurs_nStates[phn_nState_vtb_idx];

				// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
				// implemented yet, getStateValue() function in the multi-state frame node class should be
				// called instead until the multi-state segmental node class is implemented.
				float phnStateVal;
				if (this->crf->getModelType() == STDFRAME)
				{
					phnStateVal = -1 * this->nodeList->at(nodeCnt)->getStateValue(cur_phn_nState_lab);
				}
				else
				{
					phnStateVal = -1 * this->nodeList->at(nodeCnt)->getStateValue(cur_phn_nState_lab, dur);
				}


				float nState_wt = this->curViterbiNode_unPruned->viterbiWts_nStates[phn_nState_vtb_idx] + phnStateVal;
				this->curViterbiNode_unPruned->viterbiWts_nStates[phn_nState_vtb_idx] = nState_wt;
				this->curViterbiNode_unPruned->viterbiAcouWts_nStates[phn_nState_vtb_idx] += phnStateVal;
				if (nState_wt < min_wt)
				{
					min_wt = nState_wt;
					this->curViterbiNode_unPruned->viterbiWts[lat_state_vtb_idx] = min_wt;
				}
				// just for debugging
//				cout << " state " << st << ", state value = " << phnStateVal << endl;
			}
		}
		// just for debugging
//		cout << endl;

		// just for debugging
//        cout << "stateFeatureUpdate(): stateId=" << this->curViterbiNode_unPruned->viterbiStateIds[lat_state_vtb_idx]
//			<< " phnId=" << phn_lab + 1 << " wrdId=" << this->curViterbiNode_unPruned->viterbiWrdIds[lat_state_vtb_idx]
//		 	<< " min_wt=" << min_wt;
//		for (uint st = 0; st < this->nStates; st++)
//		{
//			int phn_nState_vtb_idx = lat_state_vtb_idx * this->nStates + st;
//			cout << " wts_nStates[" << st << "]=" << this->curViterbiNode_unPruned->viterbiWts_nStates[phn_nState_vtb_idx];
//			cout << " dur[" << st << "]=" << this->curViterbiNode_unPruned->viterbiDurs_nStates[phn_nState_vtb_idx];
//		}
//		cout << " isStartBound=" << this->curViterbiNode_unPruned->isPhoneStartBoundary[lat_state_vtb_idx] << endl;
	}

	// just for debugging
//	cout << "State feature update ended." << endl;
}

/*
 * CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::internalStateUpdate
 *
 * Input: nodeCnt - index into nodeList for current node updating next lab_max_dur nodes
 * 		  state_id - the state id in OpenFst lattice
 *        phn_id - phone label used in OpenFst lattice
 *        wrd_id - word label used in OpenFst lattice
 *        prevStateIdxInVtbList - previous state id in prevViterbiStateIds
 *        beam - beam width for pruning
 *
 * Performs an intra-state update of weights and backward pointers for the
 * Viterbi best-path decode.
 *
 * For a given phone label (phn_id), a given current node (nodeCnt) and a given
 * previous state(prevIdx), updates the possible paths that result in that label
 * into the current node and the next lab_max_dur-1 nodes. This update assumes a
 * multi-state model in the underlying phone structure.
 *
 */
float CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::internalStateTransUpdate(uint nodeCnt,
		uint state_id, uint phn_id, uint wrd_id, uint prevStateIdxInPrunedVtbList, double beam)
{
	time_t time1 = time(NULL);

	// just for debugging
//	cout << "Internal state trans update" << endl;

	// the start node should not invoke this function, so nodeCnt can't be 0.
	assert(nodeCnt > 0);

	int phn_lab = phn_id - 1;
	// should never be less than 0 - if there is this is a problem
	assert(phn_lab >= 0);
	float* wts_nStates = new float[this->nStates];
	int* ptrs_nStates = new int[this->nStates];
	float min_wt = 99999.0;
	bool isPhoneStartBoundary = false;

	float* acou_wts_nStates = new float[this->nStates];
	float* lm_wts_nStates = new float[this->nStates];

	// Computing the transition weights between multi-states for the same phone label (phn_lab).
	// Only two types of transitions are allowed:
	// 1. self-transition (i-th -> i-th states), and
	// 2. transitions to the next adjacent state (i-th -> (i+1)-th states).
	for (uint st = 0; st < this->nStates; st++)
	{
		int state_idx = prevStateIdxInPrunedVtbList * this->nStates + st;
		int cur_phn_nState_lab = phn_lab * this->nStates + st;

		// just for debugging
//		cout << "State " << st << endl;

		if (st == 0)
		{
			// if we are in a start state, we could get here only by self-transition.
			float old_wt = this->prevViterbiWts_nStates->at(state_idx);
			float old_acou_wt = this->prevViterbiAcouWts_nStates->at(state_idx);
			float old_lm_wt = this->prevViterbiLmWts_nStates->at(state_idx);


			// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
			// implemented yet, getTransValue() function in the multi-state frame node class should be
			// called instead until the multi-state segmental node class is implemented.
			float trans_wt = -1 * this->nodeList->at(nodeCnt)->getTransValue(cur_phn_nState_lab, cur_phn_nState_lab); //for dur=1


			wts_nStates[st] = old_wt + trans_wt;
			ptrs_nStates[st] = state_idx;
			acou_wts_nStates[st] = old_acou_wt + trans_wt;
			lm_wts_nStates[st] = old_lm_wt;

			// just for debugging
//			cout << "State " << st << " -> state " << st
//					<< ", old weight = " << old_wt << ", trans weight = " << trans_wt << endl;
		}
		else
		{
			// if we are not in a start state, we could get here by either way.

			float old_wt1 = this->prevViterbiWts_nStates->at(state_idx);
			float old_acou_wt1 = this->prevViterbiAcouWts_nStates->at(state_idx);
			float old_lm_wt1 = this->prevViterbiLmWts_nStates->at(state_idx);

			// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
			// implemented yet, getTransValue() function in the multi-state frame node class should be
			// called instead until the multi-state segmental node class is implemented.
			float trans_wt1 = -1 * this->nodeList->at(nodeCnt)->getTransValue(cur_phn_nState_lab, cur_phn_nState_lab); //for dur=1


			float new_wt1 = old_wt1 + trans_wt1;
			float new_acou_wt1 = old_acou_wt1 + trans_wt1;
			float new_lm_wt1 = old_lm_wt1;


			float old_wt2 = this->prevViterbiWts_nStates->at(state_idx - 1);
			float old_acou_wt2 = this->prevViterbiAcouWts_nStates->at(state_idx - 1);
			float old_lm_wt2 = this->prevViterbiLmWts_nStates->at(state_idx - 1);

			// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
			// implemented yet, getTransValue() function in the multi-state frame node class should be
			// called instead until the multi-state segmental node class is implemented.
			float trans_wt2 = -1 * this->nodeList->at(nodeCnt)->getTransValue(cur_phn_nState_lab - 1, cur_phn_nState_lab); //for dur=1


			float new_wt2 = old_wt2 + trans_wt2;
			float new_acou_wt2 = old_acou_wt2 + trans_wt2;
			float new_lm_wt2 = old_lm_wt2;

			// just for debugging
//			cout << "State " << st << " -> state " << st
//					<< ", old weight = " << old_wt1 << ", trans weight = " << trans_wt1 << endl;
//			cout << "State " << st - 1 << " -> state " << st
//					<< ", old weight = " << old_wt2 << ", trans weight = " << trans_wt2 << endl;

			if (new_wt1 < new_wt2)
			{
				wts_nStates[st] = new_wt1;
				ptrs_nStates[st] = state_idx;
				acou_wts_nStates[st] = new_acou_wt1;
				lm_wts_nStates[st] = new_lm_wt1;
			}
			else
			{
				wts_nStates[st] = new_wt2;
				ptrs_nStates[st] = state_idx - 1;
				acou_wts_nStates[st] = new_acou_wt2;
				lm_wts_nStates[st] = new_lm_wt2;
			}
		}

		if (wts_nStates[st] < min_wt)
		{
			min_wt = wts_nStates[st];
		}

		// just for debugging
//		cout << "State " << st << ", weight = " << wts_nStates[st] << ", min_weight = " << min_wt << endl;
	}

	int dur = 1;
	this->curViterbiNode_unPruned->addNonEpsVtbState(state_id, phn_id, wrd_id, dur, min_wt, wts_nStates, ptrs_nStates, isPhoneStartBoundary, acou_wts_nStates, lm_wts_nStates, beam);
	// for a segmental model that does not use segmental transition features, transition value
	// would be the same for the segment states in succeeding nodes if those segments share the
	// same start frame and share the same phone/word the transition is coming from and share
	// the same phone/word the transition is going to.
	for (int i = 0; i < this->lab_max_dur - 1; i++)
	{
		dur++; // this guarantees all the states being added share the same start frame
		this->nextViterbiNodes_unPruned[i]->addNonEpsVtbState(state_id, phn_id, wrd_id, dur, min_wt, wts_nStates, ptrs_nStates, isPhoneStartBoundary, acou_wts_nStates, lm_wts_nStates, beam);
	}

	delete [] wts_nStates;
	delete [] ptrs_nStates;
	delete [] acou_wts_nStates;
	delete [] lm_wts_nStates;

	time_t time2 = time(NULL);
	int_time += time2 - time1;

	// just for debugging
//	cout << "Internal state trans update ended." << endl;

	return min_wt;
}

/*
 * CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::crossStateTransUpdate
 *
 * Input: check_state - index into nodeList for current state to be updated
 *        end_idx - phone label used in OpenFst lattice
 *        transw - previous state Id
 *        phn_lab - pointer to previous state weights
 *        wrd_lab - pointer to current state weights to be updated
 *
 *
 * Performs an inter-state update of weights and backward pointers for the
 * Viterbi best-path decode.
 *
 * For a given phone label (phn_id), and a given previous state (end_idx), updates
 * the possible paths into the state that result in that label from the previous
 * state.
 *
 */
float CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::crossStateTransUpdate(uint nodeCnt,
		uint state_id, uint phn_id, uint wrd_id, uint prevPhn_EndStateIdx_InPrunedVtbList,
		uint prev_phn_id, float expand_wt, double beam)
{
	// just for debugging
//	cout << "Cross state trans update" << endl;
	//cout << "Expansion" << endl;
//	cout << "Old expansion weight = " << expand_wt;

	// non-epsilon transition
	float trans_wt;
	if (nodeCnt == 0)
	{
		// for the start node, there is not any phone transition
		// but only lattice state transition with weight exp_wt
		// which is added below.
		trans_wt = 0.0;

		// just for debugging
//		cout << " trans weight (nodeCnt=0) = " << 0.0;
	}
	else
	{
		// phone transition weight
		int prev_lab = (prev_phn_id - 1) * this->nStates + this->nStates - 1;  // the end state of the previous phone
		int cur_lab = (phn_id - 1) * this->nStates;  // the start state of the current phone

		// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
		// implemented yet, getTransValue() function in the multi-state frame node class should be
		// called instead until the multi-state segmental node class is implemented.
		trans_wt = -1 * this->nodeList->at(nodeCnt)->getTransValue(prev_lab,cur_lab); //for dur=1


		// just for debugging
//		cout << " trans weight (prev_lab=" << prev_lab << ",cur_lab=" << cur_lab << ") = " << trans_wt;
	}

	float prev_acou_wt = this->prevViterbiAcouWts_nStates->at(prevPhn_EndStateIdx_InPrunedVtbList);
	float new_lm_wt = expand_wt - prev_acou_wt;

	// add in the value on this arc
	// and all arcs travelled to get to this point
	expand_wt += trans_wt;

	float new_acou_wt = prev_acou_wt + trans_wt;

	// just for debugging
//	cout << " new expansion weight = " << expand_wt << endl;

	float* wts_nStates = new float[this->nStates];
	int* ptrs_nStates = new int[this->nStates];
	float min_wt = expand_wt;
	bool isPhoneStartBoundary = true;

	float* acou_wts_nStates = new float[this->nStates];
	float* lm_wts_nStates = new float[this->nStates];

	wts_nStates[0] = expand_wt;
	ptrs_nStates[0] = prevPhn_EndStateIdx_InPrunedVtbList;
	acou_wts_nStates[0] = new_acou_wt;
	lm_wts_nStates[0] = new_lm_wt;
	for (int i = 1; i < this->nStates; i++)
	{
		wts_nStates[i] = 99999.0;
		ptrs_nStates[i] = -1;
		acou_wts_nStates[i] = 99999.0;
		lm_wts_nStates[i] = 99999.0;
	}

	int dur = 1;
	this->curViterbiNode_unPruned->addNonEpsVtbState(state_id, phn_id, wrd_id, dur, min_wt, wts_nStates, ptrs_nStates, isPhoneStartBoundary, acou_wts_nStates, lm_wts_nStates, beam);
	// for a segmental model that does not use segmental transition features, transition value
	// would be the same for the segment states in succeeding nodes if those segments share the
	// same start frame and share the same phone/word the transition is coming from and share
	// the same phone/word the transition is going to.
	for (int i = 0; i < this->lab_max_dur - 1; i++)
	{
		dur++; // this guarantees all the states being added share the same start frame
		this->nextViterbiNodes_unPruned[i]->addNonEpsVtbState(state_id, phn_id, wrd_id, dur, min_wt, wts_nStates, ptrs_nStates, isPhoneStartBoundary, acou_wts_nStates, lm_wts_nStates, beam);
	}

	delete [] wts_nStates;
	delete [] ptrs_nStates;
	delete [] acou_wts_nStates;
	delete [] lm_wts_nStates;

	//cout << "Expansion ended" << endl;
	//cout << "Ending crossStateTransUpdate" << endl;
	//cout << "Counter: " << counter << endl;

	// just for debugging
//	cout << "Cross state trans update ended." << endl;

	return min_wt;
}

void CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::expandCrossStateFromPrevNode(uint nodeCnt, VectorFst<StdArc>* lm_fst, double beam)
{
	bool prune = true;
	if (beam <= 0.0) { prune = false; }

	deque<CRF_ViterbiState> epsList;
	map<CRF_ViterbiState,uint> epsMap;
	map<CRF_ViterbiState,uint>::iterator epsIt;

	// just for debugging
//	cout << "Current node: " << nodeCnt << endl;
//	cout << "Expand cross state from previous node." << endl;
//	cout << "Previous minimum weight: " << this->prev_min_weight << endl;
//	cout << "Beam: " << beam << endl;

	time_t time1 = time(NULL);
	for (uint prevState_listIdx = 0;
			prevState_listIdx < this->prevViterbiStateIds->size(); prevState_listIdx++) {

		uint prev_state = this->prevViterbiStateIds->at(prevState_listIdx);

		//TODO: what if nodeCnt = 0?
		uint prev_phn_id;
		if (nodeCnt == 0)
		{
			prev_phn_id = 0;
		}
		else
		{
			prev_phn_id = this->nodeList->at(nodeCnt-1)->viterbiPhnIds[prevState_listIdx];
		}

		int end_idx = prevState_listIdx * this->nStates + this->nStates - 1;
		float prev_state_exp_wt = this->prevViterbiWts_nStates->at(end_idx);

		// just for debugging
//		cout << "Previous state expansion weight: " << prev_state_exp_wt << endl;

		if (nodeCnt == 0 || prev_state_exp_wt < this->prev_min_weight + beam || !prune) {
			//only perform this update if our exit state from the previous iteration
			// is under the pruning threshold - cut down on unneeded expansions that
			// will only be pruned anyway.

			// just for debugging
//			cout << "Weight checking: within the beam." << endl;

			//uint prev_lab = (prev_phn_id - 1) * nStates + nStates - 1;

			// just for debugging
//			cout << "Starting LM iteration" << endl;

			for (ArcIterator< VectorFst<StdArc> > aiter(*lm_fst,prev_state);
							!aiter.Done();
							aiter.Next()) {
				uint check_state = aiter.Value().nextstate;
				if (aiter.Value().ilabel == 0) {

					// just for debugging
//					cout << "Epsilon handling" << endl;

					// epsilon transition - put the arc at the end of the expansionList and
					// expand it later
					float new_wt = prev_state_exp_wt + aiter.Value().weight.Value();
					// first check the epsMap and see if we already have a representative
					// for this state
//					CRF_ViterbiState test_state(check_state,prev_lab,new_wt,end_idx);
					CRF_ViterbiState test_state(check_state, prev_phn_id, new_wt, end_idx);
					epsIt = epsMap.find(test_state);
					if (epsIt == epsMap.end()) {
						epsList.push_back(test_state);
						epsMap[test_state] = epsList.size() - 1;
						epsAdded++;
					}
					else {
						uint eps_idx = epsIt->second;
						if (new_wt < epsList[eps_idx].weight) {
							epsList[eps_idx].weight = new_wt;
							epsList[eps_idx].end_idx = end_idx;
						}
						epsModified++;
					}
					// just for debugging
//					cout << "Epsilon handling ended" << endl;
				}
				else {

					// just for debugging
//					cout << "Expansion" << endl;

					checkedCounter++;

//					// non-epsilon transition
//					// First check to see if we've already done this arc on this cycle
//					// First figure out our values
//					int cur_lab = (aiter.Value().ilabel-1) * nStates;
//					float transw = -1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
//					// add in the value on this arc
//					// and all arcs travelled to get to this point
//					transw += prev_state_wt + aiter.Value().weight.Value();
//					if (transw < min_weight) {min_weight = transw;}
//					// see if we need to add it to our list of
//					// currently active transitions or if we just need to check a previously
//					// added transition
//					if (transw<min_weight+beam || !prune) {
//						crossStateUpdate_new(check_state, end_idx, transw, aiter.Value().ilabel, aiter.Value().olabel);
//					}

					uint phn_id = aiter.Value().ilabel;
					uint wrd_id = aiter.Value().olabel;
					float exp_wt = prev_state_exp_wt + aiter.Value().weight.Value();
					crossStateTransUpdate(nodeCnt, check_state, phn_id,
							wrd_id, end_idx, prev_phn_id, exp_wt, beam);

					// just for debugging
//					cout << "Expansion ended" << endl;
				}
			}
			// just for debugging
//			cout << "LM ended" << endl;
		}
		// just for debugging
//		cout << "Weight check ended" << endl;
	}
	// just for debugging
//	cout << "Cross state update ends" << endl;

	// just for debugging
//	cout << "Epsilon update" << endl;

	uint eps_offset = 0;
	while (!epsList.empty()) {
		CRF_ViterbiState eps_state = epsList.front(); epsList.pop_front();
		eps_offset++;
		epsMap.erase(eps_state);
		uint prev_state = eps_state.state;

//		uint prev_lab = eps_state.label;
		uint prev_phn_id = eps_state.label;

		float prev_state_exp_wt = eps_state.weight;
		uint end_idx = eps_state.end_idx;
		for (ArcIterator< VectorFst<StdArc> > aiter(*lm_fst,prev_state);
													!aiter.Done();
													aiter.Next()) {

			uint check_state = aiter.Value().nextstate;
			if (aiter.Value().ilabel == 0) {

				// just for debugging
//				cout << "Epsilon update: meet another epsilon." << endl;

				// epsilon transition - put the arc at the end of the expansionList and
				// expand it later
				float new_wt = prev_state_exp_wt + aiter.Value().weight.Value();
				// first check the epsMap and see if we already have a representative
				// for this state
//				CRF_ViterbiState test_state(check_state,prev_lab,new_wt,end_idx);
				CRF_ViterbiState test_state(check_state, prev_phn_id, new_wt, end_idx);
				//epsList.push_back(test_state);

				epsIt = epsMap.find(test_state);
				if (epsIt == epsMap.end()) {
					epsList.push_back(test_state);
					epsMap[test_state] = epsList.size() + eps_offset - 1;
					//adding eps_offset ensures that our list remains in synch
					epsAdded++;
				}
				else {
					uint eps_idx = epsIt->second - eps_offset;
					// we have to correct for the fact that we've removed eps_offset
					// elements from our list
					if (new_wt < epsList[eps_idx].weight) {
						epsList[eps_idx].weight = new_wt;
						epsList[eps_idx].end_idx = end_idx;
					}
					epsModified++;
				}

				//expansionList.push_back(check_state);
				//expansionWt.push_back(new_wt);
				epsilonCounter++;

				// just for debugging
//				cout << "Epsilon update: ended handling the extra epsilon." << endl;
			}
			else {

				// just for debugging
//				cout << "Epsilon update: expanding a non-epsilon arc." << endl;

				checkedCounter++;

//				// non-epsilon transition
//				// First check to see if we've already done this arc on this cycle
//				// First figure out our values
//				int cur_lab=(aiter.Value().ilabel-1)*nStates;
//				float transw=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
//				// add in the value on this arc
//				// and all arcs travelled to get to this point
//				transw+=aiter.Value().weight.Value()+exp_wt;
//				if (transw<min_weight) {min_weight=transw;}
//				// see if we need to add it to our list of
//				// currently active transitions or if we just need to check a previously
//				// added transition
//				if (transw<min_weight+beam || !prune) {
//					crossStateUpdate_new(check_state, end_idx, transw, aiter.Value().ilabel, aiter.Value().olabel);
//				}

				uint phn_id = aiter.Value().ilabel;
				uint wrd_id = aiter.Value().olabel;
				//TODO: it is needed to decide whether to use prev_phn_id or prev_lab for crossStateTransUpdate
////				uint prev_lab = (prev_phn_id - 1) * nStates + nStates - 1;
//				uint prev_phn_id = (prev_lab - nStates + 1) / nStates + 1;
				float exp_wt = prev_state_exp_wt + aiter.Value().weight.Value();
				crossStateTransUpdate(nodeCnt, check_state, phn_id,
						wrd_id, end_idx, prev_phn_id, exp_wt, beam);

				// just for debugging
//				cout << "Epsilon update: finished expanding a non-epsilon arc." << endl;
			}
		}
	}
	// just for debugging
//	cout << "Epsilon update ends" << endl;

	time_t time2 = time(NULL);
	cross_time += time2 - time1;
	epsMapSize = epsMap.size();

	// just for debugging
//	cout << "Ended Expanding cross state from previous node." << endl;
}


void CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::pruning(uint nodeCnt, double beam)
{
	bool prune = true;
	if (beam <= 0.0) { prune = false; }

	num_pruned = 0;
	// Now that we have our minimum weights at this time frame, we need to go through
	// our list again to see what we will keep and what we will prune
	this->nodeList->at(nodeCnt)->viterbiPhnIds.clear();
	this->nodeList->at(nodeCnt)->viterbiPointers.clear();
	this->nodeList->at(nodeCnt)->viterbiDurs.clear();
	this->nodeList->at(nodeCnt)->isPhoneStartBoundary.clear();
	this->nodeList->at(nodeCnt)->viterbiPhnIds.reserve(this->curViterbiNode_unPruned->viterbiPhnIds.size());
	this->nodeList->at(nodeCnt)->viterbiPointers.reserve(this->curViterbiNode_unPruned->viterbiPtrs_nStates.size());
	this->nodeList->at(nodeCnt)->viterbiDurs.reserve(this->curViterbiNode_unPruned->viterbiDurs_nStates.size());
	this->nodeList->at(nodeCnt)->isPhoneStartBoundary.reserve(this->curViterbiNode_unPruned->isPhoneStartBoundary.size());
	this->curViterbiStateIds->clear();
	this->curViterbiWrdIds->clear();
	this->curViterbiWts_nStates->clear();
	this->curViterbiAcouWts_nStates->clear();
	this->curViterbiLmWts_nStates->clear();
	this->curViterbiStateIds->reserve(this->curViterbiNode_unPruned->viterbiStateIds.size());
	this->curViterbiWrdIds->reserve(this->curViterbiNode_unPruned->viterbiWrdIds.size());
	this->curViterbiWts_nStates->reserve(this->curViterbiNode_unPruned->viterbiWts_nStates.size());
	this->curViterbiAcouWts_nStates->reserve(this->curViterbiNode_unPruned->viterbiAcouWts_nStates.size());
	this->curViterbiLmWts_nStates->reserve(this->curViterbiNode_unPruned->viterbiLmWts_nStates.size());

	//cout << "final move" << endl;
	time_t time1 = time(NULL);
	// find the minimum weight for the current node.
	for (uint idx = 0; idx < this->curViterbiNode_unPruned->viterbiStateIds.size(); idx++)
	{
		if (this->curViterbiNode_unPruned->viterbiWts[idx] < this->curViterbiNode_unPruned->min_weight)
		{
			this->curViterbiNode_unPruned->min_weight = this->curViterbiNode_unPruned->viterbiWts[idx];
		}
	}
	float min_weight = this->curViterbiNode_unPruned->min_weight;
	for (uint idx = 0; idx < this->curViterbiNode_unPruned->viterbiStateIds.size(); idx++)
	{
		// First, check to see if we want to keep this state at all
		if (this->curViterbiNode_unPruned->viterbiWts[idx] < min_weight + beam || !prune)
		{
			// We'll keep this for the next iteration

			uint state = this->curViterbiNode_unPruned->viterbiStateIds[idx];
			uint phn = this->curViterbiNode_unPruned->viterbiPhnIds[idx];
			uint wrd = this->curViterbiNode_unPruned->viterbiWrdIds[idx];
			bool isStartBound = this->curViterbiNode_unPruned->isPhoneStartBoundary[idx];

			this->curViterbiStateIds->push_back(state);
			this->curViterbiWrdIds->push_back(wrd);
			this->nodeList->at(nodeCnt)->viterbiPhnIds.push_back(phn);
			this->nodeList->at(nodeCnt)->isPhoneStartBoundary.push_back(isStartBound);

			for (uint st = 0; st < this->nStates; st++)
			{
				int nState_idx = idx * this->nStates + st;
				float nState_wt = this->curViterbiNode_unPruned->viterbiWts_nStates[nState_idx];
				int nState_ptr = this->curViterbiNode_unPruned->viterbiPtrs_nStates[nState_idx];
				uint nState_dur = this->curViterbiNode_unPruned->viterbiDurs_nStates[nState_idx];
				float nState_acou_wt = this->curViterbiNode_unPruned->viterbiAcouWts_nStates[nState_idx];
				float nState_lm_wt = this->curViterbiNode_unPruned->viterbiLmWts_nStates[nState_idx];

				// TODO: change the name of viterbiPointers to viterbiPointers_nStates and change
				// 		 the name of viterbiDurs to viterbiDurs_nStates
				this->curViterbiWts_nStates->push_back(nState_wt);
				this->nodeList->at(nodeCnt)->viterbiPointers.push_back(nState_ptr);
				this->nodeList->at(nodeCnt)->viterbiDurs.push_back(nState_dur);
				this->curViterbiAcouWts_nStates->push_back(nState_acou_wt);
				this->curViterbiLmWts_nStates->push_back(nState_lm_wt);
			}
		}
		else {
			num_pruned++;
		}
	}
	time_t time2 = time(NULL);
	merge_time += time2 - time1;
	//cout << "Final move ends" << endl;
}


/*
 * CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::nStateDecode
 *
 * Input: *result_fst - empty fst to store final result
 *        *lm_fst - fst containing the language model to decode against
 *        input_beam - beam width to decode against
 *
 *  Performs Viterbi decoding using time-synchronous pruning controlled by the
 *  input_beam parameter and constrained by the language model in lm_fst.
 *  result_fst contains a finite-state transducer with the
 *  single best word sequence as determined by the Viterbi decode.
 *
 *  NOTE:  Currently this function is hard-coded to assume that the end of
 *  sentence token (e.g. SENT_END) is mapped to symbol 1 in the language model.
 *  This should be updated to take a parameter that holds the value of the end
 *  of sentence token.
 *
 *
 */
int CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::nStateDecode(VectorFst<StdArc>* result_fst, VectorFst<StdArc>* lm_fst, double input_beam, uint min_hyps, uint max_hyps, float beam_inc)
{
	// Takes in an empty fst and an fst containing a dictionary/language model network
	// both of these should be on the tropical semiring
	// Returns a weighted fst containing the best path through the current utterance

	QNUInt32 ftr_count;
	VectorFst<StdArc>* fst = new VectorFst<StdArc>();
	StateId startState=fst->AddState();   // 1st state will be state 0 (returned by AddState)
	fst->SetStart(startState);  // arg is state ID
	bool prune = true;
	if (input_beam<=0.0) { prune=false;}

	// first we need to set up our hypothesis list
	// The "word" id vector curViterbiWrdIds gives us the list of states that we are
	// considering at this time step and is used to index the curViterbiWts and
	// curViterbiStateIds vectors.

	for (QNUInt32 i = 0; i < this->lab_max_dur; i++) {
		this->nextViterbiNodes_unPruned[i]->clear();
	}
	this->curViterbiNode_unPruned->clear();

	this->curViterbiStateIds->clear();
	this->curViterbiWrdIds->clear();
	this->curViterbiWts_nStates->clear();
	this->curViterbiAcouWts_nStates->clear();
	this->curViterbiLmWts_nStates->clear();

	this->prevViterbiStateIds->clear();
	this->prevViterbiWrdIds->clear();
	this->prevViterbiWts_nStates->clear();
	this->prevViterbiAcouWts_nStates->clear();
	this->prevViterbiLmWts_nStates->clear();

	StateIterator< VectorFst<StdArc> > lmIter(*lm_fst);
	StateId lm_start=lmIter.Value();
	this->prevViterbiStateIds->push_back(lm_start);

	// only the ending state of nStates can have cross transitions,
	// so the put a 0.0 weight for that state and infinity(99999.0)
	// for the others
	for (uint st = 0; st < this->nStates - 1; st++)
	{
		this->prevViterbiWts_nStates->push_back(99999.0);
	}
	this->prevViterbiWts_nStates->push_back(0.0);
	for (uint st = 0; st < this->nStates - 1; st++)
	{
		this->prevViterbiAcouWts_nStates->push_back(99999.0);
	}
	this->prevViterbiAcouWts_nStates->push_back(0.0);
	for (uint st = 0; st < this->nStates - 1; st++)
	{
		this->prevViterbiLmWts_nStates->push_back(99999.0);
	}
	this->prevViterbiLmWts_nStates->push_back(0.0);

	this->prev_min_weight = 0.0;

//	// The candidate phones, words, states and corresponding weights for the initial node.
//	// (only the first lab_max_dur nodes of the utterance could be an initial node)
//	vector<uint> viterbiStartPhns;
//	vector<uint> viterbiStartWrds;
//	vector<uint> viterbiStartStateIds;
//	vector<float> viterbiStartWts;
//
//	// we look for arcs off the start state and make them hypotheses
//	StateIterator< VectorFst<StdArc> > lmIter(*lm_fst);
//	StateId lm_start=lmIter.Value();
//	deque<uint> expansionList;
//	deque<float> expansionWt;
//
//	expansionList.push_back(lm_start);
//	expansionWt.push_back(0.0);
//
//	// Expand the states with epsilon input to the ones without, treated as the starting
//	// state for the initial node.
//	// (only the first lab_max_dur nodes of the utterance could be an initial node)
//	while (!expansionList.empty()) {
//		uint exp_state=expansionList.front();expansionList.pop_front();
//		float exp_wt=expansionWt.front(); expansionWt.pop_front();
//		for (ArcIterator< VectorFst<StdArc> > aiter(*lm_fst,exp_state); !aiter.Done(); aiter.Next())
//		{
//			if (aiter.Value().ilabel==0) {
//				expansionList.push_back(aiter.Value().nextstate);
//				expansionWt.push_back(aiter.Value().weight.Value()+exp_wt);
//			}
//			else {
//				viterbiStartStateIds.push_back(aiter.Value().nextstate);
//				viterbiStartPhns.push_back(aiter.Value().ilabel);
//				viterbiStartWrds.push_back(aiter.Value().olabel);
//
//				for (uint state=0; state<nStates; state++) {
//					if (state==0) {
//						// changed by Ryan, the accumulated weight should be included.
//						//viterbiStartWts.push_back(aiter.Value().weight.Value());
//						viterbiStartWts.push_back(aiter.Value().weight.Value() + exp_wt);
//					}
//					else {
//						viterbiStartWts.push_back(99999.0);
//					}
//				}
//			}
//		}
//	}

	int seq_len = 0;
	int nodeCnt = 0;
	double beam = input_beam;
	int totalStates = 0;
	num_pruned = 0;
	int_time = 0;
	cross_time = 0;
	merge_time = 0;
	time_t loopstart = time(NULL);

	// this is important! Every time bunch_size has to be reset to 1.
	// since it must start from 1 and be added by 1 every iteration up to lab_max_dur.
	this->bunch_size = 1;
	do {
		ftr_count=this->ftr_strm->read(this->bunch_size,this->ftr_buf,this->lab_buf);

		if (ftr_count > 0) {

			// just for debugging
//			cout << "before creating new_buf." << endl;

			QNUInt32 cur_ftr_buf_size = this->num_ftrs * ftr_count;

			// Just for debugging
//			cout << "QNUInt32 cur_ftr_buf_size = num_ftrs * ftr_count = " << num_ftrs << " * " << ftr_count << " = " << cur_ftr_buf_size << endl;

			float* new_buf = new float[cur_ftr_buf_size];
			for (QNUInt32 j = 0; j < cur_ftr_buf_size; j++) {
				new_buf[j] = this->ftr_buf[j];
				//cout << " " << new_buf[j];
			}

			QNUInt32 label = CRF_LAB_BAD;
			//TODO: use this->labs_width>0 or this->lab_buf != NULL ?
			if (this->labs_width > 0)
			{
				//TODO: design a labmap class to do the label mapping.
				QNUInt32 actualLab = CRF_LAB_BAD;
//				QNUInt32 ifBrokenLab = CRF_LAB_BAD;
//				for (QNUInt32 dur = 1; dur <= ftr_count; dur++)
//				{
//					actualLab = this->lab_buf[this->labs_width * (dur - 1)];
//					if (actualLab == CRF_LAB_BAD)
//						continue;
//					QNUInt32 begin = this->lab_buf[this->labs_width * (dur - 1) + 1];
//					QNUInt32 end = this->lab_buf[this->labs_width * (dur - 1) + 2];
//					assert(end - begin + 1 == dur);
//					//ifBrokenLab = this->lab_buf[this->labs_width * (dur - 1) + 3];
//					if (label != CRF_LAB_BAD)
//					{
//						string errstr="CRF_LatticeBuilder::buildLattice() caught exception: two valid labels found for the same frame.";
//						throw runtime_error(errstr);
//					}
//					label = this->nActualLabs * (dur - 1) + actualLab;
//					//label = this->nActualLabs * 2 * (dur - 1) + actualLab * 2 + ifBrokenLab;
//				}
					actualLab = this->lab_buf[0];

					// just for debugging
//					cout << "actualLab=" << actualLab << endl;

					if (actualLab != CRF_LAB_BAD)
					{
						QNUInt32 begin = this->lab_buf[1];
						QNUInt32 end = this->lab_buf[2];

						QNUInt32 dur = end - begin + 1;

						// just for debugging
//						cout << "begin=" << begin << " end=" << end;

						// for phn_dur label
						label = this->nActualLabs * (dur - 1) + actualLab;

						// for phn_dur_brokenClasses2 label
//						ifBrokenLab = this->lab_buf[3];
//						label = nActualLabs * (dur - 1) + actualLab * 2 + ifBrokenLab;
					}

					// just for debugging
//					cout << " label=" << label << endl;
			}

			QNUInt32 nodeMaxDur;
			if (nodeCnt + 1 <= this->lab_max_dur)
			{
				nodeMaxDur = nodeCnt + 1;
			}
			else
			{
				nodeMaxDur = this->lab_max_dur;
			}
			// TODO: For segmental CRFs which have different structures for different nodes, these parameters need to be changed.
			QNUInt32 prevNode_nLabs = this->crf->getNLabs();
			QNUInt32 nextNode_nActualLabs = this->nActualLabs;


			// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
			// implemented yet, for the multi-state case, we have to initiate the multi-state frame node
			// class instead of multi-state, for single-state case, we could initiate the frame node
			// class or segmental node class.
			if (this->crf->getModelType() == STDFRAME)
			{
				this->nodeList->set(nodeCnt,new_buf,cur_ftr_buf_size,label,this->crf);
			}
			else
			{
				this->nodeList->set(nodeCnt,new_buf,cur_ftr_buf_size,label,this->crf,nodeMaxDur,prevNode_nLabs,nextNode_nActualLabs);
			}


			// just for debugging
//			cout << "Label: " << label << endl;

			QNUInt32 numPrevNodes;
			if (nodeCnt + 1 <= this->lab_max_dur)
			{
				numPrevNodes = nodeCnt;
			}
			else
			{
				numPrevNodes = this->lab_max_dur;
			}
			assert(numPrevNodes + 1 == nodeMaxDur || numPrevNodes == nodeMaxDur);
			CRF_StateNode** prevNodes = NULL;
			if (numPrevNodes > 0)
			{

				// just for debugging
//				cout << "before creating prevNodes." << endl;

				prevNodes = new CRF_StateNode*[numPrevNodes];
				for (QNUInt32 i = 0; i < numPrevNodes; i++)
				{
					QNUInt32 ni = nodeCnt - numPrevNodes + i;
					prevNodes[i] = this->nodeList->at(ni);
				}
			}
			this->nodeList->at(nodeCnt)->setPrevNodes(prevNodes, numPrevNodes);

			// just for debugging
			//cout << "Before computing state and transition matrix for node " << nodeCnt << ": ";
			//int pauseTemp;
			//cin >> pauseTemp;

			float value = this->nodeList->at(nodeCnt)->computeTransMatrix();


			// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
			// implemented yet, computeAlpha() function in the multi-state frame node class should be
			// called instead until the multi-state segmental node class is implemented.
			double scale;
			if (this->crf->getModelType() == STDFRAME)
			{
				double* prev_alpha;
				if (nodeCnt == 0) {
					prev_alpha=this->alpha_base;
					scale=this->nodeList->at(nodeCnt)->computeFirstAlpha(prev_alpha);
				}
				else {
					prev_alpha=this->nodeList->at(nodeCnt-1)->getAlpha();
					scale=this->nodeList->at(nodeCnt)->computeAlpha(prev_alpha);
				}
			}
			else
			{
				if (nodeCnt == 0) {
					scale = this->nodeList->at(nodeCnt)->computeFirstAlpha();
				}
				else {
					scale = this->nodeList->at(nodeCnt)->computeAlpha();
				}
			}


			seq_len++;

			// Now we do our forward viterbi processing
			// First make sure the various bookkeeping vectors are empty before we begin.
			//cout << "Clearing viterbi bookkeeping structs" << endl;
//			this->nodeList->at(nodeCnt)->viterbiPhnIds.clear();
//			this->nodeList->at(nodeCnt)->viterbiPointers.clear();
//			this->curViterbiStateIds->clear();
//			this->curViterbiWrdIds->clear();
//			this->curViterbiWts_nStates->clear();
			epsAdded=0;
			epsModified=0;
			epsMapSize=0;

			epsilonCounter=0;
			checkedCounter=0;
//			addedCounter=0;
//			updateCheckCounter=0;
//			updateCounter=0;
//			dupCounter=0;

			// Cross State Transition Weight Update
			expandCrossStateFromPrevNode(nodeCnt, lm_fst, beam);

			// Internal State Transition Feature Update
			// The start node does not have internal state transition. So nodeCnt starts from 1.
			if (nodeCnt > 0)
			{
				for (uint idx = 0; idx < this->prevViterbiStateIds->size(); idx++)
				{
					uint state_id = this->prevViterbiStateIds->at(idx);
					uint phn_id = this->nodeList->at(nodeCnt - 1)->viterbiPhnIds[idx];
					uint wrd_id = this->prevViterbiWrdIds->at(idx);
					uint prevStateIdxInPrunedVtbList = idx;
					internalStateTransUpdate(nodeCnt, state_id, phn_id, wrd_id, prevStateIdxInPrunedVtbList, beam);
				}
			}

			// State Feature Update
			stateFeatureUpdate(nodeCnt);

			// Pruning
			pruning(nodeCnt, beam);



			// Printing

			totalStates += this->curViterbiStateIds->size();

			// As a check on where we are, let's dump the current ViterbiStates and Weights at every
			// cycle.  Just to see if what we're doing looks sane
			double average = (totalStates / (nodeCnt + 1.0));

			// changed by Ryan
			//if (nodeCnt%50 == 0) {
			if (nodeCnt%100 == 0) {
//			if (nodeCnt%1 == 0) {
				cout << "time: " << nodeCnt << " " << this->curViterbiStateIds->size() << " hyps ";
				//cout << "time: " << nodeCnt << " " << curViterbiStateIds->size() << " hyps ";
				cout << this->curViterbiWts_nStates->size() << " nStates ";
				cout << this->prevViterbiWts_nStates->size() << " prev_nStates";
				cout << " pruned away " << num_pruned << " this cycle" << endl;
				cout << "       checked: " << checkedCounter;
				cout << " added: " << this->curViterbiNode_unPruned->addedCounter;
				cout << " updated: " << this->curViterbiNode_unPruned->updateCounter << endl;
				cout << "       checkUpdate: " << this->curViterbiNode_unPruned->updateCheckCounter;
				cout << " epsilons: " << epsilonCounter << endl;
				cout << "      Average hyps per timestep: " << average;
				cout << "      Current pruning beam: " << beam << endl;
				cout << " epsAdded: " << epsAdded << "  epsModified: " << epsModified << endl;
				cout << " epsMapSize: " << epsMapSize << endl;
				cout << "Internal update time: " << int_time;
				cout << " Cross update time: " << cross_time;
				cout << " Merge time: " << merge_time << endl;

				// Added by Ryan
				cout << "minimum weight = " << this->curViterbiNode_unPruned->min_weight << endl;

				// Added by Ryan
//				cout << "Cross update: " << crossUpdateSamePhoneCounter << " same phones, "
//						<< crossUpdateDiffPhoneCounter << " different phones, "
//						<< crossUpdateSameWordCounter << " same words, "
//						<< crossUpdateDiffWordCounter << " different words, "
//						<< crossUpdateSameStateCounter << " same states, "
//						<< crossUpdateDiffStateCounter << " different states."
//						<< endl;
			}



			// swap the current and previous pruned vectors

			vector<uint>* holdStates = this->curViterbiStateIds;
			this->curViterbiStateIds = this->prevViterbiStateIds;
			this->prevViterbiStateIds = holdStates;

			vector<uint>* holdWrdIds = this->curViterbiWrdIds;
			this->curViterbiWrdIds = this->prevViterbiWrdIds;
			this->prevViterbiWrdIds = holdWrdIds;

			vector<float>* holdWts = this->curViterbiWts_nStates;
			this->curViterbiWts_nStates = this->prevViterbiWts_nStates;
			this->prevViterbiWts_nStates = holdWts;

			vector<float>* holdAcouWts = this->curViterbiAcouWts_nStates;
			this->curViterbiAcouWts_nStates = this->prevViterbiAcouWts_nStates;
			this->prevViterbiAcouWts_nStates = holdAcouWts;

			vector<float>* holdLmWts = this->curViterbiLmWts_nStates;
			this->curViterbiLmWts_nStates = this->prevViterbiLmWts_nStates;
			this->prevViterbiLmWts_nStates = holdLmWts;

			this->prev_min_weight = this->curViterbiNode_unPruned->min_weight;

			this->curViterbiStateIds->clear();
			this->curViterbiWrdIds->clear();
			this->curViterbiWts_nStates->clear();
			this->curViterbiAcouWts_nStates->clear();
			this->curViterbiLmWts_nStates->clear();

			// pop the first node of next unpruned nodes deque as the current unpruned node
			// push a new unpruned node to the back of the deque
			CRF_ViterbiNode* holdNode = this->curViterbiNode_unPruned;
			this->curViterbiNode_unPruned = this->nextViterbiNodes_unPruned.front();
			this->nextViterbiNodes_unPruned.pop_front();
			holdNode->clear();
			this->nextViterbiNodes_unPruned.push_back(holdNode);



			// just for debugging
//			cout << endl;

			if (prevNodes != NULL)
				delete [] prevNodes;
			nodeCnt++;
		}

		// bunch_size (number of windows ending at next frame) is added by 1 in each iteration, until being equal to lab_max_dur.
		if (this->bunch_size < this->lab_max_dur)
			this->bunch_size++;

	} while (ftr_count > 0);
	time_t loopend = time(NULL);
	cout << "Internal update time: " << int_time;
	cout << " Cross update time: " << cross_time;
	cout << " Merge time: " << merge_time << endl;
	cout << " Total loop time: " << (loopend - loopstart) << endl;
	this->nodeList->setNodeCount(nodeCnt);


	// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
	// implemented yet, computeAlphaSum() function in the multi-state frame node class should be
	// called instead until the multi-state segmental node class is implemented.
	double Zx = this->nodeList->at(nodeCnt - 1)->computeAlphaSum();


	int final_state = fst->AddState();
	fst->SetFinal(final_state, Zx);

	// Backtrack through the decoder to build the output lattice
	// First search through the final timestep and find the best result that ends in a final state
	float min_weight = 99999.0;
	int min_idx = -1;
	cout << "checking through " << this->prevViterbiStateIds->size() << " states for final state " << endl;
	int fstate_cnt = 0;
	for (uint idx = 0; idx < this->prevViterbiStateIds->size(); idx++) {
		uint tmp_wrd = this->prevViterbiWrdIds->at(idx);
		if (tmp_wrd == 1) {
			fstate_cnt++;
			int end_nState_idx = idx * this->nStates + this->nStates - 1;
			float weight = this->prevViterbiWts_nStates->at(end_nState_idx);
			if (weight < min_weight) {
				min_weight = weight;
				min_idx = end_nState_idx;
			}
		}
	}

	// changed by Ryan
	//cout << "Found " << fstate_cnt << " final states " << min_idx << endl;
	cout << "Found " << fstate_cnt << " final states " << min_idx << " with weight " << min_weight << endl;

	// Build the phone lattice, then compose it with the lm lattice.
	int cur_state=fst->AddState();
	int next_state=final_state;
	if (min_idx < 0) {
		// Failsafe - create a null transition
		cout << "ERROR: Could not reach end of utterance" << endl;
		fst->AddArc(startState,StdArc(0,0,8,final_state));

		// output the acoustic model weight and the language model weight
		cout << "Acoustic model weight (negative log potential) = " << 0
				<< ", -Z(X) = " << 0
				<< ", language model weight (negative log probability) = " << 0 << endl;

	}
	else {

		// output the acoustic model weight and the language model weight
		cout << "Acoustic model weight (negative log potential) = " << this->prevViterbiAcouWts_nStates->at(min_idx)
				<< ", -Z(X) = " << -1 * Zx
				<< ", language model weight (negative log probability) = " << this->prevViterbiLmWts_nStates->at(min_idx) << endl;

		// added by Ryan, just for debugging
		cout << "Reverse label sequence: ";

		int seg_end_idx = nodeCnt - 1;
		while(seg_end_idx >= 0)
//		for (int idx = nodeCnt - 1; idx >= 0; idx--)
		{
			int cur_min_idx_ptr = min_idx / this->nStates;
			int cur_min_idx_offset = min_idx % this->nStates;
			uint cur_lab = this->nodeList->at(seg_end_idx)->viterbiPhnIds[cur_min_idx_ptr] - 1;
			uint cur_dur = this->nodeList->at(seg_end_idx)->viterbiDurs[cur_min_idx_ptr];
			uint cur_lab_nState = cur_lab * this->nStates + cur_min_idx_offset;

			assert(cur_dur > 0 && cur_dur <= this->lab_max_dur);

			int seg_start_idx = seg_end_idx + 1 - cur_dur;

			if (seg_start_idx < 0)
			{
				string errstr="CRF_ViterbiDecoder_StdSeg_NoSegTransFtr::nStateDecode() caught exception: "
						"current segment end index = " + stringify(seg_end_idx) + ", start index = " + stringify(seg_start_idx) +
						", duration = " + stringify(cur_dur);
				throw runtime_error(errstr);
			}
			else if (seg_start_idx == 0)
			{
				// the current node is the hypothesized start segment of the current utterance

				// we're in our first timestep, so we add an arc that just has the state value here
				//cout << idx << " cur_lab: " << cur_lab << endl;

				// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
				// implemented yet, getStateValue() function in the multi-state frame node class should be
				// called instead until the multi-state segmental node class is implemented.
				float trans_w;
				if (this->crf->getModelType() == STDFRAME)
				{
					trans_w = -1 * this->nodeList->at(seg_end_idx)->getStateValue(cur_lab_nState);
				}
				else
				{
					trans_w = -1 * this->nodeList->at(seg_end_idx)->getStateValue(cur_lab_nState, cur_dur);
				}

				fst->AddArc(startState, StdArc(cur_lab_nState + 1, cur_lab + 1, trans_w, next_state));

				// added by Ryan, just for debugging
				cout << " " << cur_lab;
			}
			else
			{
				// seg_start_idx > 0
				// which means the current segment is not the start segment of the utterance

				int prev_seg_end_idx = seg_start_idx - 1;

				// just for debugging
				cout << "seg_end_idx = " << seg_end_idx << ", cur_dur = " << cur_dur
						<< ", cur_lab = " << cur_lab << ", cur_lab_nState = " << cur_lab_nState
						<< ", prev_seg_end_idx = " << prev_seg_end_idx << endl;

				assert(prev_seg_end_idx >= 0);

				int prev_min_idx = this->nodeList->at(seg_end_idx)->viterbiPointers[min_idx];
				int prev_min_idx_ptr = prev_min_idx / this->nStates;
				int prev_min_idx_offset = prev_min_idx % this->nStates;
				uint prev_lab = this->nodeList->at(prev_seg_end_idx)->viterbiPhnIds[prev_min_idx_ptr] - 1;
				uint prev_lab_nState = prev_lab * this->nStates + prev_min_idx_offset;

				// TODO: By Ryan, this is tricky now. Since the multi-state segmental node class has not been
				// implemented yet, getFullTransValue() function in the multi-state frame node class should be
				// called instead until the multi-state segmental node class is implemented.
				float trans_w;
				if (this->crf->getModelType() == STDFRAME)
				{
					trans_w = -1 * this->nodeList->at(seg_end_idx)->getFullTransValue(prev_lab_nState, cur_lab_nState);
				}
				else
				{
					trans_w = -1 * this->nodeList->at(seg_end_idx)->getFullTransValue(prev_lab_nState, cur_lab_nState, cur_dur);
				}

				// Changed by Ryan
				//
				// The original way does not apply to the case when nStates == 1
				// So we have to use the new way to include that case as well.
				//
				// ********* This is the old way. *********
//				// Now we add the arc.  We need to be careful about how we label it.  If we're making
//				// a transition from an end state to a start state, then we need to put the label for
//				// the new phone on this arc
//				if (cur_min_idx_offset==0 && prev_min_idx_offset==nStates-1) {
				// ****************************************
				//
				// ********* This is the new way. *********
				// Now we add the arc. We put the label for the new phone on this arc only when the
				// hypothesized phone state of the current node is the start boundary of the current
				// phone, which also means the hypothesized phone state of the previous node is the
				// end boundary of the previous phone.
				if (cur_min_idx_offset == 0 && this->nodeList->at(seg_end_idx)->isPhoneStartBoundary[cur_min_idx_ptr]) {
				// ****************************************
					//cout << idx << " cur_lab: " << cur_lab << endl;
					fst->AddArc(cur_state, StdArc(cur_lab_nState + 1, cur_lab + 1, trans_w, next_state));

					// added by Ryan, just for debugging
					cout << cur_lab << endl;
				}
				else {
					//cout << idx << " cur_lab: " << cur_lab << endl;
					fst->AddArc(cur_state, StdArc(cur_lab_nState + 1, 0, trans_w, next_state));
				}
				min_idx = prev_min_idx;
				next_state = cur_state;
				cur_state = fst->AddState();
			}
			for (int in_frame_id = seg_end_idx; in_frame_id >= seg_start_idx; in_frame_id--)
			{
				this->nodeList->at(in_frame_id)->viterbiPhnIds.clear();
				this->nodeList->at(in_frame_id)->viterbiPointers.clear();
				this->nodeList->at(in_frame_id)->isPhoneStartBoundary.clear();
				this->nodeList->at(in_frame_id)->viterbiDurs.clear();
			}
			seg_end_idx = seg_start_idx - 1; // move to the previous segment
			seg_start_idx = -1;  // clear seg_start_idx
		}

		// added by Ryan, just for debugging
		cout << endl;
	}
	Connect(fst);
	Compose(*fst,*lm_fst,result_fst);

	// added by Ryan, just for debugging
	cout << "result_fst has " << result_fst->NumStates() << " states." << endl;


	for (QNUInt32 i = 0; i < this->lab_max_dur; i++) {
		this->nextViterbiNodes_unPruned[i]->clear();
	}
	this->curViterbiNode_unPruned->clear();

	this->curViterbiStateIds->clear();
	this->curViterbiWrdIds->clear();
	this->curViterbiWts_nStates->clear();
	this->curViterbiAcouWts_nStates->clear();
	this->curViterbiLmWts_nStates->clear();

	this->prevViterbiStateIds->clear();
	this->prevViterbiWrdIds->clear();
	this->prevViterbiWts_nStates->clear();
	this->prevViterbiAcouWts_nStates->clear();
	this->prevViterbiLmWts_nStates->clear();
	this->prev_min_weight = 99999.0;


	delete fst;

	return nodeCnt;

}
