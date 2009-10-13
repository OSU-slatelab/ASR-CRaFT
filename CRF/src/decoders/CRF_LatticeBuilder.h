#ifndef CRF_LATTICEBUILDER_H_
#define CRF_LATTICEBUILDER_H_

#include "fst/lib/fstlib.h"
#include "../CRF.h"
#include "../CRF_Model.h"
#include "../io/CRF_FeatureStream.h"
#include "../nodes/CRF_StateVector.h"

/*#include "CRF_StdStateVectorLog.h"
#include "CRF_StdNStateVectorLog.h"*/

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
	double* alpha_base;
	virtual void _computeGamma(vector<double> *stategamma,vector<double> *transgamma,Fst<LogArc>& fst,int nstates);

public:
	CRF_LatticeBuilder(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in);
	virtual ~CRF_LatticeBuilder();
	virtual StdVectorFst* testBuild();
	virtual StdVectorFst* buildLattice();
	template <class Arc> int buildLattice(VectorFst<Arc>*fst, bool align=false, VectorFst<Arc>*alignFst=NULL,bool norm=true);
	template <class Arc> int nStateBuildLattice(VectorFst<Arc>*fst, bool align=false, VectorFst<Arc>*alignFst=NULL, bool norm=true);
	virtual void computeAlignedAlphaBeta(Fst<LogArc>& fst, int nstates);
	virtual StdVectorFst* bestPath(bool align);
	virtual StdVectorFst* bestPath_old(bool align);
	virtual StdVectorFst* LMBestPath(bool align, StdFst* lmFst);
	virtual int getAlignmentGammas(vector<double> *denominatorStateGamma,
									vector<double> *numeratorStateGamma,
									vector<double> *denominatorTransGamma,
									vector<double> *numeratorTransGamma);
	virtual int getAlignmentGammasNState(vector<double> *denominatorStateGamma,
									vector<double> *numeratorStateGamma,
									vector<double> *denominatorTransGamma,
									vector<double> *numeratorTransGamma);
	virtual StdVectorFst* nStateBuildLattice();
	virtual StdVectorFst* nStateBuildLattice(StdVectorFst* labFst);
	virtual StdVectorFst* nStateBestPath(bool align);
	virtual StdVectorFst* nStateLMBestPath(bool align, StdFst* lmFst);
	virtual StdVectorFst* nStateBestPath_old(bool align);
	virtual CRF_StateVector* getNodeList();
};

// These template functions are included here for compile-time generation

// this version requires you to create fst and labFst
template <class Arc> int CRF_LatticeBuilder::buildLattice(VectorFst<Arc>* fst,
									  bool align,
									  VectorFst<Arc>*labFst,
									  bool norm) {
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

	do {
		ftr_count=ftr_strm->read(this->bunch_size,ftr_buf,lab_buf);

		for (QNUInt32 i=0; i<ftr_count; i++) {
			float* new_buf = new float[this->num_ftrs];
			for (QNUInt32 j=0; j<this->num_ftrs; j++) {
				int idx=i*this->num_ftrs+j;
				new_buf[j]=ftr_buf[idx];
			}
			this->nodeList->set(nodeCnt,new_buf,num_ftrs,this->lab_buf[i],this->crf);
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
	//double Zx=this->nodeList->at(nodeCnt-1)->computeAlphaSum();
	if (norm) { Zx=-1*this->nodeList->at(nodeCnt-1)->computeAlphaSum(); }
	//cout << "Zx: " << Zx << endl;
	int final_state = fst->AddState();
	for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
		int prev_state=(this->num_labs)*(nodeCnt-1)+(prev_lab+1);
		//fst->AddArc(prev_state,Arc(0,0,0,final_state));
		fst->AddArc(prev_state,Arc(0,0,Zx,final_state));
	}
	fst->SetFinal(final_state,0);
	if (align) {
		labFst->SetFinal(curLabState,0);
	}
	return seq_len;
}

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
