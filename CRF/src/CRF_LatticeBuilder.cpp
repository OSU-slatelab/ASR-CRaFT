#include "CRF_LatticeBuilder.h"

CRF_LatticeBuilder::CRF_LatticeBuilder(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in)
	: crf(crf_in),
	  ftr_strm(ftr_strm_in)
{
	if (this->crf->getFeatureMap()->getNumStates() == 1) {
		this->nodeList = new CRF_StdStateVectorLog();
	}
	else {
		this->nodeList = new CRF_StdNStateVectorLog();
	}
	this->bunch_size=1;
	this->num_ftrs=this->ftr_strm->num_ftrs();
	this->num_labs=this->crf->getNLabs();
	this->ftr_buf = new float[num_ftrs*bunch_size];
	this->lab_buf = new QNUInt32[bunch_size];
}

CRF_LatticeBuilder::~CRF_LatticeBuilder()
{
}

StdVectorFst* CRF_LatticeBuilder::testBuild()
{
	// A vector FST is a general mutable FST
	StdVectorFst* fst = new StdVectorFst();

	// Adds state 0 to the initially empty FST and make it the start state.
	fst->AddState();   // 1st state will be state 0 (returned by AddState)
	fst->SetStart(0);  // arg is state ID

	// Adds two arcs exiting state 0.
	// Arc constructor args: ilabel, olabel, weight, dest state ID.
	fst->AddArc(0, StdArc(1, 1, 0.5, 1));  // 1st arg is src state ID
	fst->AddArc(0, StdArc(2, 2, 1.5, 1));

	// Adds state 1 and its arc.
	fst->AddState();
	fst->AddArc(1, StdArc(3, 3, 2.5, 2));

	// Adds state 2 and set its final weight.
	fst->AddState();
	fst->SetFinal(2, 3.5);  // 1st arg is state ID, 2nd arg weight
	return fst;
}

StdVectorFst* CRF_LatticeBuilder::buildLattice()
{
		// Returns the best path through the current segment
	QNUInt32 ftr_count;
	StdVectorFst* fst = new StdVectorFst();
	int startState=0;
	fst->AddState();   // 1st state will be state 0 (returned by AddState)
	fst->SetStart(startState);  // arg is state ID


	int seq_len=0;
	QNUInt32 nodeCnt=0;
	//cout << "Starting Viterbi Best Path: " << endl;
	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		ftr_count=ftr_strm->read(this->bunch_size,ftr_buf,lab_buf);

		//cout << "Feature Count: " << ftr_count << endl;
		for (QNUInt32 i=0; i<ftr_count; i++) {
			//cout << "\tLabel: " << lab_buf[i] << "\tFeas:";
			// Now, separate the bunch into individual frames
			float* new_buf = new float[this->num_ftrs];
			for (QNUInt32 j=0; j<this->num_ftrs; j++) {
				int idx=i*this->num_ftrs+j;
				new_buf[j]=ftr_buf[idx];
				//cout << " " << new_buf[j];
			}
			//cout << endl;
			// Store the current frame/label information in a sequence node
			//	* sequence nodes create a doubly-linked list, with the previous node known at creation time
			//next_seq = new CRF_Seq(new_buf,this->num_ftrs,lab_buf[i],crf->getNLabs(),cur_seq);
			this->nodeList->set(nodeCnt,new_buf,num_ftrs,this->lab_buf[i],this->crf);
			seq_len++;
			/*seq_len++;
			if (seq_head == NULL)
			{
				seq_head=next_seq;
			}
			if (cur_seq != NULL) {
				cur_seq->setNext(next_seq);
			}
			cur_seq=next_seq;*/
			//float value=this->gb->computeTransMatrix(cur_seq,crf);
			//float value=this->gb->computeTransMatrixLog(cur_seq);
			float value=this->nodeList->at(nodeCnt)->computeTransMatrix();
			if (nodeCnt==startState) {
				// Add arcs from the startState to each possible label
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					float value=-1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
					int cur_state=fst->AddState();
					fst->AddArc(startState,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
				}
			}
			else {
				int cur_time=nodeCnt+1;
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					int cur_state=fst->AddState();
					for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
						float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						int prev_state=(this->num_labs)*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
					}
				}
			}
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	int final_state = fst->AddState();
	for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
		int prev_state=(this->num_labs)*(nodeCnt-1)+(prev_lab+1);
		fst->AddArc(prev_state,StdArc(0,0,0,final_state));
	}
	fst->SetFinal(final_state,0);

	return fst;

}

StdVectorFst* CRF_LatticeBuilder::bestPath(bool align)
{
		// Returns the best path through the current segment
	QNUInt32 ftr_count;
	StdVectorFst* fst = new StdVectorFst();
	StdVectorFst* labFst = new StdVectorFst();
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
	//cout << "Starting Viterbi Best Path: " << endl;
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
			if (nodeCnt==startState) {
				// Add arcs from the startState to each possible label
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					float value=-1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
					int cur_state=fst->AddState();
					fst->AddArc(startState,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
				}
			}
			else {
				int cur_time=nodeCnt+1;
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					int cur_state=fst->AddState();
					for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
						float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						int prev_state=(this->num_labs)*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
					}
				}
			}
			if (align) {
				QNUInt32 lab=this->nodeList->at(nodeCnt)->getLabel()+1;
				if (firstLab or (lab != curLab)) {
					int prevLabState=curLabState;
					curLabState=labFst->AddState();
					labFst->AddArc(prevLabState,StdArc(lab,lab,0,curLabState));
					labFst->AddArc(curLabState,StdArc(lab,lab,0,curLabState)); // Add self loop
					firstLab=false;
					curLab=lab;
				}
			}

			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	int final_state = fst->AddState();
	for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
		int prev_state=(this->num_labs)*(nodeCnt-1)+(prev_lab+1);
		fst->AddArc(prev_state,StdArc(0,0,0,final_state));
	}
	fst->SetFinal(final_state,0);
	StdVectorFst* final_result=new StdVectorFst();

	if (align) {
		labFst->SetFinal(curLabState,0);
		StdComposeFst* result=new StdComposeFst(*fst,*labFst);
		ShortestPath(*result,final_result,1);
		delete result;
	}
	else {
		ShortestPath(*fst,final_result,1);
	}
	Project(final_result,PROJECT_OUTPUT);
	RmEpsilon(final_result);
	TopSort(final_result);

	delete fst;
	delete labFst;

	return final_result;

}

StdVectorFst* CRF_LatticeBuilder::LMBestPath(bool align, StdFst* lmFst)
{
	StdVectorFst* fst = this->buildLattice();
	StdVectorFst* final_result=new StdVectorFst();

	StdComposeFst* result=new StdComposeFst(*fst,*lmFst);
	ShortestPath(*result,final_result,1);
	delete result;
	Project(final_result,PROJECT_OUTPUT);
	RmEpsilon(final_result);
	TopSort(final_result);
	delete fst;

	return final_result;
}


StdVectorFst* CRF_LatticeBuilder::nStateBuildLattice()
{
	QNUInt32 nStates = this->crf->getFeatureMap()->getNumStates();
		// Returns the best path through the current segment
	QNUInt32 ftr_count;
	StdVectorFst* fst = new StdVectorFst();
	int startState=0;
	fst->AddState();   // 1st state will be state 0 (returned by AddState)
	fst->SetStart(startState);  // arg is state ID


	int seq_len=0;
	QNUInt32 nodeCnt=0;
	//cout << "Starting Viterbi Best Path: " << endl;
	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		ftr_count=ftr_strm->read(this->bunch_size,ftr_buf,lab_buf);

		//cout << "Feature Count: " << ftr_count << endl;
		for (QNUInt32 i=0; i<ftr_count; i++) {
			//cout << "\tLabel: " << lab_buf[i] << "\tFeas:";
			// Now, separate the bunch into individual frames
			float* new_buf = new float[this->num_ftrs];
			for (QNUInt32 j=0; j<this->num_ftrs; j++) {
				int idx=i*this->num_ftrs+j;
				new_buf[j]=ftr_buf[idx];
				//cout << " " << new_buf[j];
			}
			//cout << endl;
			// Store the current frame/label information in a sequence node
			//	* sequence nodes create a doubly-linked list, with the previous node known at creation time
			//next_seq = new CRF_Seq(new_buf,this->num_ftrs,lab_buf[i],crf->getNLabs(),cur_seq);
			this->nodeList->set(nodeCnt,new_buf,num_ftrs,this->lab_buf[i],this->crf);
			seq_len++;
			/*seq_len++;
			if (seq_head == NULL)
			{
				seq_head=next_seq;
			}
			if (cur_seq != NULL) {
				cur_seq->setNext(next_seq);
			}
			cur_seq=next_seq;*/
			//float value=this->gb->computeTransMatrix(cur_seq,crf);
			//float value=this->gb->computeTransMatrixLog(cur_seq);
			float value=this->nodeList->at(nodeCnt)->computeTransMatrix();
			if (nodeCnt==startState) {
				// Add arcs from the startState to each possible STARTING label
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					// The above should be modified so that we can only start in a start state, but for
					// now ignore this and let our dictionary take care of it (this makes computing the
					// previous state easier later
					float value=-1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
					int cur_state=fst->AddState();
					fst->AddArc(startState,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
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
							fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						}
						// Special case handling - if we have more than one state we have to explicitly
						// put a self loop in
						if (nStates >1) {
							float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(cur_lab,cur_lab);
							int prev_state=(this->num_labs)*(cur_time-2)+(cur_lab+1);
							fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						}
					}
					else {
						// We're not in a start state - all we need are arcs from the previous label
						// and arcs from the previous self label
						int prev_lab=cur_lab-1;
						float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						int prev_state=this->num_labs*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						prev_lab=cur_lab;
						value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						prev_state=this->num_labs*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
					}
				}
			}
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	int final_state = fst->AddState();
	for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
		int prev_state=(this->num_labs)*(nodeCnt-1)+(prev_lab+1);
		fst->AddArc(prev_state,StdArc(0,0,0,final_state));
	}
	fst->SetFinal(final_state,0);

	return fst;

}

StdVectorFst* CRF_LatticeBuilder::nStateBuildLattice(StdVectorFst* labFst)
{
	// Note: Passing a NULL value in for labFst says not to create an alignment label Fst

	QNUInt32 nStates = this->crf->getFeatureMap()->getNumStates();
	// Returns the best path through the current segment
	QNUInt32 ftr_count;
	StdVectorFst* fst = new StdVectorFst();
//	StdVectorFst* labFst = new StdVectorFst();
	int startState=0;
	int labStartState=0;
	QNUInt32 curLab=0;
	int curLabState=labStartState;
	bool firstLab=true;
	fst->AddState();   // 1st state will be state 0 (returned by AddState)
	fst->SetStart(startState);  // arg is state ID
	if (labFst != NULL) {
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
			if (nodeCnt==startState) {
				// Add arcs from the startState to each possible STARTING label
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					// The above should be modified so that we can only start in a start state, but for
					// now ignore this and let our dictionary take care of it (this makes computing the
					// previous state easier later
					float value=-1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
					int cur_state=fst->AddState();
					fst->AddArc(startState,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
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
							fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						}
						// Special case handling - if we have more than one state we have to explicitly
						// put a self loop in
						if (nStates >1) {
							float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(cur_lab,cur_lab);
							int prev_state=(this->num_labs)*(cur_time-2)+(cur_lab+1);
							fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						}
					}
					else {
						// We're not in a start state - all we need are arcs from the previous label
						// and arcs from the previous self label
						int prev_lab=cur_lab-1;
						float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						int prev_state=this->num_labs*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						prev_lab=cur_lab;
						value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						prev_state=this->num_labs*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
					}
				}
			}
			if (labFst != NULL) {
				QNUInt32 lab=this->nodeList->at(nodeCnt)->getLabel()+1;
				if (firstLab or (lab != curLab)) {
					int prevLabState=curLabState;
					curLabState=labFst->AddState();
					labFst->AddArc(prevLabState,StdArc(lab,lab,0,curLabState));
					labFst->AddArc(curLabState,StdArc(lab,lab,0,curLabState)); // Add self loop
					firstLab=false;
					curLab=lab;
				}
			}
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	int final_state = fst->AddState();
	for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
		int prev_state=(this->num_labs)*(nodeCnt-1)+(prev_lab+1);
		fst->AddArc(prev_state,StdArc(0,0,0,final_state));
	}
	fst->SetFinal(final_state,0);
	StdVectorFst* final_result=new StdVectorFst();

	if (labFst != NULL ) {
		labFst->SetFinal(curLabState,0);
	}

	return fst;

}

StdVectorFst* CRF_LatticeBuilder::nStateBestPath(bool align)
{
	StdVectorFst* labFst;
	if (align) {
		labFst = new StdVectorFst();
	}
	else {
		labFst = NULL;
	}
	StdVectorFst* fst = this->nStateBuildLattice(labFst);

	StdVectorFst* final_result=new StdVectorFst();

	if (align) {

		StdComposeFst* result=new StdComposeFst(*fst,*labFst);
		ShortestPath(*result,final_result,1);
		delete result;
	}
	else {
		ShortestPath(*fst,final_result,1);
	}
//	Project(final_result,PROJECT_OUTPUT);
	RmEpsilon(final_result);
	TopSort(final_result);
	if (labFst != NULL ) {
		delete labFst;
	}
	delete fst;

	return final_result;

}

StdVectorFst* CRF_LatticeBuilder::nStateLMBestPath(bool align, StdFst* lmFst)
{
	StdVectorFst* fst = this->nStateBuildLattice(NULL);
	StdVectorFst* final_result=new StdVectorFst();

	StdComposeFst* result=new StdComposeFst(*fst,*lmFst);
	ShortestPath(*result,final_result,1);
	delete result;
	Project(final_result,PROJECT_OUTPUT);
	RmEpsilon(final_result);
	TopSort(final_result);
	delete fst;

	return final_result;
}

StdVectorFst* CRF_LatticeBuilder::nStateBestPath_old(bool align)
{
	QNUInt32 nStates = this->crf->getFeatureMap()->getNumStates();
		// Returns the best path through the current segment
	QNUInt32 ftr_count;
	StdVectorFst* fst = new StdVectorFst();
	StdVectorFst* labFst = new StdVectorFst();
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
			if (nodeCnt==startState) {
				// Add arcs from the startState to each possible STARTING label
				for (int cur_lab=0; cur_lab<this->num_labs; cur_lab++) {
					// The above should be modified so that we can only start in a start state, but for
					// now ignore this and let our dictionary take care of it (this makes computing the
					// previous state easier later
					float value=-1*this->nodeList->at(nodeCnt)->getStateValue(cur_lab);
					int cur_state=fst->AddState();
					fst->AddArc(startState,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
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
							fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						}
						// Special case handling - if we have more than one state we have to explicitly
						// put a self loop in
						if (nStates >1) {
							float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(cur_lab,cur_lab);
							int prev_state=(this->num_labs)*(cur_time-2)+(cur_lab+1);
							fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						}
					}
					else {
						// We're not in a start state - all we need are arcs from the previous label
						// and arcs from the previous self label
						int prev_lab=cur_lab-1;
						float value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						int prev_state=this->num_labs*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
						prev_lab=cur_lab;
						value=-1*this->nodeList->at(nodeCnt)->getFullTransValue(prev_lab,cur_lab);
						prev_state=this->num_labs*(cur_time-2)+(prev_lab+1);
						fst->AddArc(prev_state,StdArc(cur_lab+1,cur_lab+1,value,cur_state));
					}
				}
			}
			if (align) {
				QNUInt32 lab=this->nodeList->at(nodeCnt)->getLabel()+1;
				if (firstLab or (lab != curLab)) {
					int prevLabState=curLabState;
					curLabState=labFst->AddState();
					labFst->AddArc(prevLabState,StdArc(lab,lab,0,curLabState));
					labFst->AddArc(curLabState,StdArc(lab,lab,0,curLabState)); // Add self loop
					firstLab=false;
					curLab=lab;
				}
			}
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	int final_state = fst->AddState();
	for (int prev_lab=0; prev_lab<this->num_labs; prev_lab++) {
		int prev_state=(this->num_labs)*(nodeCnt-1)+(prev_lab+1);
		fst->AddArc(prev_state,StdArc(0,0,0,final_state));
	}
	fst->SetFinal(final_state,0);
	StdVectorFst* final_result=new StdVectorFst();


	if (align) {
		labFst->SetFinal(curLabState,0);
		StdComposeFst* result=new StdComposeFst(*fst,*labFst);
		ShortestPath(*result,final_result,1);
		delete result;
	}
	else {
		ShortestPath(*fst,final_result,1);
	}
	Project(final_result,PROJECT_OUTPUT);
	RmEpsilon(final_result);
	TopSort(final_result);
	delete labFst;
	delete fst;

	return final_result;

}
