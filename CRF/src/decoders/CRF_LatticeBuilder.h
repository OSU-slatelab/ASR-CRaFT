#ifndef CRF_LATTICEBUILDER_H_
#define CRF_LATTICEBUILDER_H_
/*
 * CRF_LatticeBuilder.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "fst/fstlib.h"
#include "../CRF.h"
#include "../CRF_Model.h"
#include "../io/CRF_FeatureStream.h"
#include "../nodes/CRF_StateVector.h"

using namespace fst;

/*
 * class CRF_LatticeBuilder
 *
 * Used to create an OpenFST lattice for decoding.  Instantiating a
 * CRF_LatticeBuilder object requires an input feature stream
 * (see CRF_FeatureStream) and a CRF model (see CRF_Model).
 * Lattices are built for a feature sequence using a call to one of the
 * buildLattice functions.
 */

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
	double* alpha_base;
	virtual void _computeGamma(vector<double> *stategamma,vector<double> *transgamma,Fst<LogArc>& fst,int nstates);

	// Added by Ryan
	size_t labs_width;
	QNUInt32 lab_max_dur;
	QNUInt32 ftr_buf_size;
	QNUInt32 lab_buf_size;
	QNUInt32 nActualLabs;
	vector<QNUInt32>* nodeStartStates;

public:
	CRF_LatticeBuilder(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_LatticeBuilder();
	template <class Arc> int buildLattice(VectorFst<Arc>*fst, bool align=false, VectorFst<Arc>*alignFst=NULL,bool norm=true);
	template <class Arc> int nStateBuildLattice(VectorFst<Arc>*fst, bool align=false, VectorFst<Arc>*alignFst=NULL, bool norm=true);
	virtual StdVectorFst* buildLattice();
	virtual int getAlignmentGammas(vector<double> *denominatorStateGamma,
									vector<double> *numeratorStateGamma,
									vector<double> *denominatorTransGamma,
									vector<double> *numeratorTransGamma);
	virtual int getAlignmentGammasNState(vector<double> *denominatorStateGamma,
									vector<double> *numeratorStateGamma,
									vector<double> *denominatorTransGamma,
									vector<double> *numeratorTransGamma);
	virtual CRF_StateVector* getNodeList();
	virtual void computeAlignedAlphaBeta(Fst<LogArc>& fst, int nstates);
};

// These template functions are included here for compile-time generation

/*
 * CRF_LatticeBuilder::buildLattice
 *
 * Input: *fst - pointer to VectorFst where result should be stored
 *        align - true if resulting fst should be aligned to input labels
 *        *labFst - pointer to a VectorFst where label lattice should be stored
 *                  (for alignment)
 *        norm - true if the normalization constant should be computed and
 *               placed on the final arc of the lattice
 *
 * Returns: number of observation nodes in the sequence from the input stream
 *
 * Reads the current sequence of feature observations from the input feature
 * stream and uses them to construct an OpenFst lattice for CRF best-path
 * processing.  If the align flag is set, the resulting lattice is composed with
 * an fst created from the label sequence and only the paths that match the
 * label sequence taken from the input feature stream is returned.
 *
 */
template <class Arc> int CRF_LatticeBuilder::buildLattice(VectorFst<Arc>* fst,
									  bool align,
									  VectorFst<Arc>*labFst,
									  bool norm) {

	// just for debugging
	//cout << "Beginning of CRF_LatticeBuilder::buildLattice";

	// Returns the best path through the current segment
	QNUInt32 ftr_count;

	int startState=0;
	int labStartState=0;
	QNUInt32 curLab=0;
	int curLabState=labStartState;
	bool firstLab=true;
	fst->AddState();   // 1st state will be state 0 (returned by AddState)
	fst->SetStart(startState);  // arg is state ID
	if (align) {
		labFst->AddState(); // 1st state will be state 0 (resturned by AddState);
		labFst->SetStart(labStartState); // arg is stateID
	}

	int seq_len=0;
	QNUInt32 nodeCnt=0;

	// Changed by Ryan
#ifndef SEGMENTAL_CRF
	do {
		ftr_count=ftr_strm->read(this->bunch_size,ftr_buf,lab_buf);

		for (QNUInt32 i=0; i<ftr_count; i++) {
			float* new_buf = new float[this->num_ftrs];
			for (QNUInt32 j=0; j<this->num_ftrs; j++) {
				int idx=i*this->num_ftrs+j;
				new_buf[j]=ftr_buf[idx];
			}
			this->nodeList->set(nodeCnt,new_buf,num_ftrs,this->lab_buf[i],this->crf);

			// Added by Ryan, just for debugging
			//cout << "Before computing state and transition matrix for node " << nodeCnt << ": ";
			//int pauseTemp;
			//cin >> pauseTemp;

			float value=this->nodeList->at(nodeCnt)->computeTransMatrix();
			double scale;
			double* prev_alpha;
			if (nodeCnt == 0) {
				prev_alpha=this->alpha_base;
				scale=this->nodeList->at(nodeCnt)->computeFirstAlpha(prev_alpha);
			}
			else {
				prev_alpha=this->nodeList->at(nodeCnt-1)->getAlpha();
				scale=this->nodeList->at(nodeCnt)->computeAlpha(prev_alpha);
			}
			seq_len++;
			if (nodeCnt==startState) {
				// Add arcs from the startState to each possible label
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					float value=-1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
					int cur_state=fst->AddState();
					fst->AddArc(startState,Arc(cur_lab+1,cur_lab+1,value,cur_state));

					// just for debugging
					//cout << "AddArc(" << startState << ",Arc(" << cur_lab+1 << "," << cur_lab+1 << "," << value << "," << cur_state << "));" << endl;
				}
			}
			else {
				int cur_time=nodeCnt+1;
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					int cur_state=fst->AddState();
					for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
						float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						int prev_state=(this->num_labs)*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,Arc(cur_lab+1,cur_lab+1,value,cur_state));

						// just for debugging
						//cout << "AddArc(" << prev_state << ",Arc(" << cur_lab+1 << "," << cur_lab+1 << "," << value << "," << cur_state << "));" << endl;
					}
				}
			}
			if (align) {
				QNUInt32 lab=this->nodeList->at(nodeCnt)->getLabel()+1;
				if (firstLab or (lab != curLab)) {
					int prevLabState=curLabState;
					curLabState=labFst->AddState();
					labFst->AddArc(prevLabState,Arc(lab,lab,0,curLabState));
					labFst->AddArc(curLabState,Arc(lab,lab,0,curLabState)); // Add self loop
					firstLab=false;
					curLab=lab;
				}
			}

			// just for debugging
			//cout << endl;

			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	this->nodeList->setNodeCount(nodeCnt);
	double Zx=0;
	//double Zx=this->nodeList->at(nodeCnt-1)->computeAlphaSum();
	if (norm) { Zx=-1*this->nodeList->at(nodeCnt-1)->computeAlphaSum(); }
	//cout << "Zx: " << Zx << endl;
	int final_state = fst->AddState();
	for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
		int prev_state=(this->num_labs)*(nodeCnt-1)+(prev_lab+1);
		//fst->AddArc(prev_state,Arc(0,0,0,final_state));
		fst->AddArc(prev_state,Arc(0,0,Zx,final_state));

		// just for debugging
		//cout << "AddArc(" << prev_state << ",Arc(" << 0 << "," << 0 << "," << Zx << "," << final_state << "));" << endl;
	}
	fst->SetFinal(final_state,0);
	if (align) {
		labFst->SetFinal(curLabState,0);
	}
#else
	do {
		ftr_count=ftr_strm->read(this->bunch_size,ftr_buf,lab_buf);

		if (ftr_count > 0) {

			// just for debugging
//			cout << "before creating new_buf." << endl;

			QNUInt32 cur_ftr_buf_size = num_ftrs * ftr_count;

			// Just for debugging
//			cout << "QNUInt32 cur_ftr_buf_size = num_ftrs * ftr_count = " << num_ftrs << " * " << ftr_count << " = " << cur_ftr_buf_size << endl;

			float* new_buf = new float[cur_ftr_buf_size];
			for (QNUInt32 j=0; j<cur_ftr_buf_size; j++) {
				new_buf[j]=this->ftr_buf[j];
				//cout << " " << new_buf[j];
			}

			QNUInt32 label = CRF_LAB_BAD;
			//TODO: use this->labs_width>0 or this->lab_buf != NULL ?
			if (this->labs_width > 0)
			{
				//TODO: design a labmap class to do the label mapping.
				QNUInt32 actualLab = CRF_LAB_BAD;
				//QNUInt32* ifBrokenLab = CRF_LAB_BAD;
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

						//ifBrokenLab = this->lab_buf[3];
						label = nActualLabs * (dur - 1) + actualLab;
						//label = nActualLabs * 2 * (dur - 1) + actualLab * 2 + ifBrokenLab;
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
			// TODO: For segmental CRFs which have the different structures for different nodes, these parameters need to be changed.
			QNUInt32 prevNode_nLabs = this->crf->getNLabs();
			QNUInt32 nextNode_nActualLabs = this->nActualLabs;
			this->nodeList->set(nodeCnt,new_buf,cur_ftr_buf_size,label,this->crf,nodeMaxDur,prevNode_nLabs,nextNode_nActualLabs);

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
			double scale;
			if (nodeCnt == 0) {
			// just for debugging
			//	scale = this->nodeList->at(nodeCnt)->computeFirstAlpha();
			}
			else {
			// just for debugging
			//	scale = this->nodeList->at(nodeCnt)->computeAlpha();
			}
			seq_len++;
			if (nodeCnt == startState) {

				// just for debugging
//				cout << "beginning of the first node: " << nodeCnt << endl;

				assert(numPrevNodes == 0 && nodeMaxDur == 1);
				// Add arcs from the startState to each possible label
				int cur_lab = 0;
				QNUInt32 num_new_states = 0;

				//***** CRF_StdSegStateNode_WithoutDurLab *****//
				bool createNewNodes = true;
				//***********************************//
				int cur_state;
				//***** for CRF_StdSegStateNode *****//
				for (int dur = numPrevNodes + 1; dur <= nodeMaxDur; dur++)
				{
				//***********************************//
					for (int lab = 0; lab < this->nActualLabs; lab++) {

						//***** for CRF_StdSegStateNode *****//
						//cur_state = fst->AddState();
						//num_new_states++;
						//***********************************//

						//***** CRF_StdSegStateNode_WithoutDurLab *****//
						if (createNewNodes)
						{
							cur_state = fst->AddState();
							num_new_states++;
						}
						//***********************************//

					//***** CRF_StdSegStateNode_WithoutDurLab *****//
					//for (int dur = numPrevNodes + 1; dur <= nodeMaxDur; dur++)
					//{
					//***********************************//

						// for CRF_StdSegStateNode
						//float value = -1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
						// for CRF_StdSegStateNode_WithoutDurLab
						float value = -1*this->nodeList->at(nodeCnt)->getStateValue(lab, dur);

						// just for debugging
//						cout << "nodeCnt=" << nodeCnt << ", dur=" << dur << ", lab=" << lab << endl;
//						cout << "Arc value = -1*this->nodeList->at(" << nodeCnt << ")->getStateValue(lab=" << lab << ", dur=" << dur << ")=" << value << endl;

						fst->AddArc(startState,Arc(cur_lab+1,cur_lab+1,value,cur_state));

						// just for debugging
//						cout << "AddArc(" << startState << ",Arc(" << cur_lab+1 << "," << cur_lab+1 << "," << value << "," << cur_state << "));" << endl;

						cur_lab++;

						//***** CRF_StdSegStateNode_WithoutDurLab *****//
						cur_state++;
						//***********************************//
					}
					//***** CRF_StdSegStateNode_WithoutDurLab *****//
					cur_state -= this->nActualLabs;
					createNewNodes = false;
					//***********************************//
				}

				// just for debugging
//				cout << "To add start state " << startState << " to node " << nodeCnt << endl;

				assert(nodeCnt <= this->nodeStartStates->size());
				if (nodeCnt < this->nodeStartStates->size()) {
					this->nodeStartStates->at(nodeCnt) = startState + 1;
				} else {
					this->nodeStartStates->push_back(startState + 1);
				}

				QNUInt32 nextNodeCnt = nodeCnt + 1;
				QNUInt32 nextNodeStartState = this->nodeStartStates->at(nodeCnt) + num_new_states;

				// just for debugging
//				cout << "To add start state " << nextNodeStartState << " to node " << nextNodeCnt << endl;

				assert(nextNodeCnt <= this->nodeStartStates->size());
				if (nextNodeCnt < this->nodeStartStates->size()) {
					this->nodeStartStates->at(nextNodeCnt) = nextNodeStartState;
				} else {
					this->nodeStartStates->push_back(nextNodeStartState);
				}

				// just for debugging
//				cout << "end of the first node." << nodeCnt << endl;
			}
			else if (numPrevNodes < nodeMaxDur) {

				// just for debugging
//				cout << "beginning of one of the starting nodes: " << nodeCnt << endl;

				int cur_time = nodeCnt + 1;
				int cur_lab = 0;
				QNUInt32 num_new_states = 0;

				//***** CRF_StdSegStateNode_WithoutDurLab *****//
				bool createNewNodes = true;
				//***********************************//
				int cur_state;
				//***** for CRF_StdSegStateNode *****//
				for (int dur = 1; dur <= numPrevNodes; dur++)
				{
				//***********************************//
					for (int lab = 0; lab < this->nActualLabs; lab++)
					{
						//***** for CRF_StdSegStateNode *****//
						//cur_state = fst->AddState();
						//num_new_states++;
						//***********************************//

						//***** CRF_StdSegStateNode_WithoutDurLab *****//
						if (createNewNodes)
						{
							cur_state = fst->AddState();
							num_new_states++;
						}
						//***********************************//

					//***** CRF_StdSegStateNode_WithoutDurLab *****//
					//for (int dur = 1; dur <= numPrevNodes; dur++)
					//{
					//*********************************************//

						int prev_numAvailLabs = this->nodeList->at(nodeCnt-dur)->getNumAvailLabs();
						for (int prev_lab = 0; prev_lab < prev_numAvailLabs; prev_lab++) {
							// for CRF_StdSegStateNode
							//float value = -1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
							// for CRF_StdSegStateNode_WithoutDurLab
							float value = -1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,lab,dur);

							// just for debugging
//							cout << "nodeCnt=" << nodeCnt << ", dur=" << dur << ", lab=" << lab << ", prev_lab=" << prev_lab << endl;
//							cout << "Arc value = -1*this->nodeList->at(" << nodeCnt << ")->getFullTransValue(prev_lab=" << prev_lab << ", lab=" << lab << ", dur=" << dur << ")=" << value << endl;

							int prev_state = this->nodeStartStates->at(nodeCnt-dur) + prev_lab;
							fst->AddArc(prev_state,Arc(cur_lab+1,cur_lab+1,value,cur_state));

							// just for debugging
//							cout << "AddArc(" << prev_state << ",Arc(" << cur_lab+1 << "," << cur_lab+1 << "," << value << "," << cur_state << "));" << endl;
						}
						cur_lab++;

						//***** CRF_StdSegStateNode_WithoutDurLab *****//
						cur_state++;
						//***********************************//
					}
					//***** CRF_StdSegStateNode_WithoutDurLab *****//
					cur_state -= this->nActualLabs;
					createNewNodes = false;
					//***********************************//
				}
				//***** for CRF_StdSegStateNode *****//
				for (int dur = numPrevNodes + 1; dur <= nodeMaxDur; dur++)
				{
				//***********************************//
					for (int lab = 0; lab < this->nActualLabs; lab++) {

						//***** for CRF_StdSegStateNode *****//
						//cur_state = fst->AddState();
						//num_new_states++;
						//***********************************//

						//***** CRF_StdSegStateNode_WithoutDurLab *****//
						if (createNewNodes)
						{
							cur_state = fst->AddState();
							num_new_states++;
						}
						//***********************************//

					//***** CRF_StdSegStateNode_WithoutDurLab *****//
					//for (int dur = numPrevNodes + 1; dur <= nodeMaxDur; dur++)
					//{
					//***********************************//

						// for CRF_StdSegStateNode
						//float value = -1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
						// for CRF_StdSegStateNode_WithoutDurLab
						float value = -1*this->nodeList->at(nodeCnt)->getStateValue(lab, dur);

						// just for debugging
//						cout << "nodeCnt=" << nodeCnt << ", dur=" << dur << ", lab=" << lab << endl;
//						cout << "Arc value = -1*this->nodeList->at(" << nodeCnt << ")->getStateValue(lab=" << lab << ", dur=" << dur << ")=" << value << endl;

						fst->AddArc(startState,Arc(cur_lab+1,cur_lab+1,value,cur_state));

						// just for debugging
//						cout << "AddArc(" << startState << ",Arc(" << cur_lab+1 << "," << cur_lab+1 << "," << value << "," << cur_state << "));" << endl;

						cur_lab++;

						//***** CRF_StdSegStateNode_WithoutDurLab *****//
						cur_state++;
						//***********************************//
					}
					//***** CRF_StdSegStateNode_WithoutDurLab *****//
					cur_state -= this->nActualLabs;
					createNewNodes = false;
					//***********************************//
				}

				QNUInt32 nextNodeCnt = nodeCnt + 1;
				QNUInt32 nextNodeStartState = this->nodeStartStates->at(nodeCnt) + num_new_states;

				// just for debugging
//				cout << "To add start state " << nextNodeStartState << " to node " << nextNodeCnt << endl;

				assert(nextNodeCnt <= this->nodeStartStates->size());
				if (nextNodeCnt < this->nodeStartStates->size()) {
					this->nodeStartStates->at(nextNodeCnt) = nextNodeStartState;
				} else {
					this->nodeStartStates->push_back(nextNodeStartState);
				}

				// just for debugging
//				cout << "end of one of the starting nodes: " << nodeCnt << endl;
			}
			else {

				// just for debugging
//				cout << "beginning of a regular node: " << nodeCnt << endl;

				int cur_time = nodeCnt + 1;
				int cur_lab = 0;
				QNUInt32 num_new_states = 0;

				//***** CRF_StdSegStateNode_WithoutDurLab *****//
				bool createNewNodes = true;
				//***********************************//
				int cur_state;
				//***** for CRF_StdSegStateNode *****//
				for (int dur = 1; dur <= nodeMaxDur; dur++)
				{
				//***********************************//
					for (int lab = 0; lab < this->nActualLabs; lab++)
					{
						//***** for CRF_StdSegStateNode *****//
						//cur_state = fst->AddState();
						//num_new_states++;
						//***********************************//

						//***** CRF_StdSegStateNode_WithoutDurLab *****//
						if (createNewNodes)
						{
							cur_state = fst->AddState();
							num_new_states++;
						}
						//***********************************//

					//***** CRF_StdSegStateNode_WithoutDurLab *****//
					//for (int dur = 1; dur <= nodeMaxDur; dur++)
					//{
					//***********************************//

						int prev_numAvailLabs = this->nodeList->at(nodeCnt-dur)->getNumAvailLabs();
						for (int prev_lab = 0; prev_lab < prev_numAvailLabs; prev_lab++) {
							// for CRF_StdSegStateNode
							//float value = -1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
							// for CRF_StdSegStateNode_WithoutDurLab
							float value = -1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,lab,dur);

							// just for debugging
//							cout << "nodeCnt=" << nodeCnt << ", dur=" << dur << ", lab=" << lab << ", prev_lab=" << prev_lab << endl;
//							cout << "Arc value = -1*this->nodeList->at(" << nodeCnt << ")->getFullTransValue(prev_lab=" << prev_lab << ", lab=" << lab << ", dur=" << dur << ")=" << value << endl;

							int prev_state = this->nodeStartStates->at(nodeCnt-dur) + prev_lab;
							fst->AddArc(prev_state,Arc(cur_lab+1,cur_lab+1,value,cur_state));

							// just for debugging
//							cout << "AddArc(" << prev_state << ",Arc(" << cur_lab+1 << "," << cur_lab+1 << "," << value << "," << cur_state << "));" << endl;
						}
						cur_lab++;

						//***** CRF_StdSegStateNode_WithoutDurLab *****//
						cur_state++;
						//***********************************//
					}
					//***** CRF_StdSegStateNode_WithoutDurLab *****//
					cur_state -= this->nActualLabs;
					createNewNodes = false;
					//***********************************//
				}

				QNUInt32 nextNodeCnt = nodeCnt + 1;
				QNUInt32 nextNodeStartState = this->nodeStartStates->at(nodeCnt) + num_new_states;

				// just for debugging
//				cout << "To add start state " << nextNodeStartState << " to node " << nextNodeCnt << endl;

				assert(nextNodeCnt <= this->nodeStartStates->size());
				if (nextNodeCnt < this->nodeStartStates->size()) {
					this->nodeStartStates->at(nextNodeCnt) = nextNodeStartState;
				} else {
					this->nodeStartStates->push_back(nextNodeStartState);
				}

				// just for debugging
//				cout << "end of a regular node: " << nodeCnt << endl;
			}
			if (align) {
				if (this->labs_width == 0)
				{
					string errstr="CRF_LatticeBuilder::buildLattice() caught exception: The label stream is found NULL or the label width is found 0 under the align mode.";
					throw runtime_error(errstr);
				}
				QNUInt32 lab=this->nodeList->at(nodeCnt)->getLabel()+1;
				if (firstLab or (lab != curLab)) {
					int prevLabState=curLabState;
					curLabState=labFst->AddState();
					labFst->AddArc(prevLabState,Arc(lab,lab,0,curLabState));
					labFst->AddArc(curLabState,Arc(lab,lab,0,curLabState)); // Add self loop
					firstLab=false;
					curLab=lab;
				}
			}

			// just for debugging
//			cout << endl;

			if (prevNodes != NULL)
				delete [] prevNodes;
			nodeCnt++;
		}
	} while (ftr_count > 0);

	this->nodeList->setNodeCount(nodeCnt);
	double Zx = 0;
	//double Zx=this->nodeList->at(nodeCnt-1)->computeAlphaSum();
	if (norm) { Zx = -1*this->nodeList->at(nodeCnt-1)->computeAlphaSum(); }
	//cout << "Zx: " << Zx << endl;
	int final_state = fst->AddState();
	int lastNodeNumAvailLabs = this->nodeList->at(nodeCnt-1)->getNumAvailLabs();
	for (int prev_lab = 0; prev_lab < lastNodeNumAvailLabs; prev_lab++) {
		int prev_state = this->nodeStartStates->at(nodeCnt-1) + prev_lab;
		//fst->AddArc(prev_state,Arc(0,0,0,final_state));
		fst->AddArc(prev_state,Arc(0,0,Zx,final_state));

		// just for debugging
//		cout << "Zx=" << Zx << ", final_state=" << final_state << endl;
//		cout << "AddArc(" << prev_state << ",Arc(" << 0 << "," << 0 << "," << Zx << "," << final_state << "));" << endl;

	}
	fst->SetFinal(final_state,0);
	if (align) {
		labFst->SetFinal(curLabState,0);
	}
#endif

	return seq_len;
}

/*
 * CRF_LatticeBuilder::nStateBuildLattice
 *
 * Input: *fst - pointer to VectorFst where result should be stored
 *        align - true if resulting fst should be aligned to input labels
 *        *labFst - pointer to a VectorFst where label lattice should be stored
 *                  (for alignment)
 *        norm - true if the normalization constant should be computed and
 *               placed on the final arc of the lattice
 *
 * Returns: number of observation nodes in the sequence from the input stream
 *
 * Reads the current sequence of feature observations from the input feature
 * stream and uses them to construct an OpenFst lattice for CRF best-path
 * processing, specially constructed to pass through n states, where "n" is
 * set in the CRF_Model given to the CRF_LatticeBuilder at instantiation.
 * If the align flag is set, the resulting lattice is composed with
 * an fst created from the label sequence and only the paths that match the
 * label sequence taken from the input feature stream is returned.
 *
 * Much of this code follows the form of buildLattice (above), with the addition
 * of the n-state restrictions on the allowed paths.
 */
template <class Arc> int CRF_LatticeBuilder::nStateBuildLattice(VectorFst<Arc>* fst,
									  bool align,
									  VectorFst<Arc>*labFst,
									  bool norm) {

	QNUInt32 nStates = this->crf->getFeatureMap()->getNumStates();
	// Returns the best path through the current segment
	QNUInt32 ftr_count;
//	StdVectorFst* fst = new StdVectorFst();
//	StdVectorFst* labFst = new StdVectorFst();
	int startState=0;
	int labStartState=0;
	QNUInt32 curLab=0;
	int curLabState=labStartState;
	bool firstLab=true;
	fst->AddState();   // 1st state will be state 0 (returned by AddState)
	fst->SetStart(startState);  // arg is state ID
	if (align) {
		labFst->AddState(); // 1st state will be state 0 (resturned by AddState);
		labFst->SetStart(labStartState); // arg is stateID
	}


	int seq_len=0;
	QNUInt32 nodeCnt=0;
	do {
		ftr_count=ftr_strm->read(this->bunch_size,ftr_buf,lab_buf);

		for (QNUInt32 i=0; i<ftr_count; i++) {
			float* new_buf = new float[this->num_ftrs];
			for (QNUInt32 j=0; j<this->num_ftrs; j++) {
				int idx=i*this->num_ftrs+j;
				new_buf[j]=ftr_buf[idx];
			}
			this->nodeList->set(nodeCnt,new_buf,num_ftrs,this->lab_buf[i],this->crf);
			seq_len++;
			float value=this->nodeList->at(nodeCnt)->computeTransMatrix();
			double scale;
			double* prev_alpha;
			if (nodeCnt == 0) {
				prev_alpha=this->alpha_base;
				scale=this->nodeList->at(nodeCnt)->computeFirstAlpha(prev_alpha);
			}
			else {
				prev_alpha=this->nodeList->at(nodeCnt-1)->getAlpha();
				scale=this->nodeList->at(nodeCnt)->computeAlpha(prev_alpha);
			}
			if (nodeCnt==startState) {
				// Add arcs from the startState to each possible STARTING label
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					// The above should be modified so that we can only start in a start state, but for
					// now ignore this and let our dictionary take care of it (this makes computing the
					// previous state easier later
					float value=-1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
					int cur_state=fst->AddState();
					fst->AddArc(startState,Arc(cur_lab+1,cur_lab+1,value,cur_state));
				}
			}
			else {
				int cur_time=nodeCnt+1;
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					int cur_state=fst->AddState();
					if (cur_lab % nStates ==0) {
						// We're in a start state - add arcs from all possible previous end states
						for (int prev_lab=nStates-1; prev_lab<this->num_labs; prev_lab+=nStates) {
							float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
							int prev_state=(this->num_labs)*(cur_time-2)+(prev_lab+1);
							fst->AddArc(prev_state,Arc(cur_lab+1,cur_lab+1,value,cur_state));
						}
						// Special case handling - if we have more than one state we have to explicitly
						// put a self loop in
						if (nStates >1) {
							float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(cur_lab,cur_lab);
							int prev_state=(this->num_labs)*(cur_time-2)+(cur_lab+1);
							fst->AddArc(prev_state,Arc(cur_lab+1,cur_lab+1,value,cur_state));
						}
					}
					else {
						// We're not in a start state - all we need are arcs from the previous label
						// and arcs from the previous self label
						int prev_lab=cur_lab-1;
						float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						int prev_state=this->num_labs*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,Arc(cur_lab+1,cur_lab+1,value,cur_state));
						prev_lab=cur_lab;
						value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						prev_state=this->num_labs*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,Arc(cur_lab+1,cur_lab+1,value,cur_state));
					}
				}
			}
			if (align) {
				QNUInt32 lab=this->nodeList->at(nodeCnt)->getLabel()+1;
				if (firstLab or (lab != curLab)) {
					int prevLabState=curLabState;
					curLabState=labFst->AddState();
					labFst->AddArc(prevLabState,Arc(lab,lab,0,curLabState));
					labFst->AddArc(curLabState,Arc(lab,lab,0,curLabState)); // Add self loop
					firstLab=false;
					curLab=lab;
				}
			}
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	this->nodeList->setNodeCount(nodeCnt);
	double Zx=0;
	//double Zx=-1*this->nodeList->at(nodeCnt-1)->computeAlphaSum();
	if (norm) { Zx=-1*this->nodeList->at(nodeCnt-1)->computeAlphaSum(); }
	int final_state = fst->AddState();
	for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
		int prev_state=(this->num_labs)*(nodeCnt-1)+(prev_lab+1);
		//fst->AddArc(prev_state,StdArc(0,0,0,final_state));
		fst->AddArc(prev_state,Arc(0,0,Zx,final_state));
	}
	fst->SetFinal(final_state,0);
	//StdVectorFst* final_result=new StdVectorFst();

	if (align ) {
		labFst->SetFinal(curLabState,0);
	}

	return seq_len;
}

#endif /*CRF_LATTICEBUILDER_H_*/
