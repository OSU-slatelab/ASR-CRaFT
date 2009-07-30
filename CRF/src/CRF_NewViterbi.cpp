#include "CRF_NewViterbi.h"

CRF_NewViterbi::CRF_NewViterbi(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in)
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

CRF_NewViterbi::~CRF_NewViterbi()
{
}

CRF_LabelPath* CRF_NewViterbi::bestPath()
{

	if (this->crf->getFeatureMap()->getNumStates() != 1) {
		//return this->nStateBestPath();
	}
	// Returns the best path through the current segment
	QNUInt32 ftr_count;

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
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);

	//cout << "Reached end of utterance read" << endl;

	nodeCnt--; // Correct for the fact that we add one at the end of the loop
	// Transition matrices have been computed and we know how long the sequence is.
	// Now we can perform the Viterbi algorithm over the transition matrices

	typedef float* floatPtr;
	typedef double* doublePtr;
	typedef int* intPtr;

	//floatPtr* delta = new floatPtr[seq_len+1];
	doublePtr* delta = new doublePtr[seq_len+1];
	intPtr* psi = new intPtr[seq_len+1];

	for (int idx=0; idx<seq_len+1; idx++) {
		delta[idx] = new double[this->num_labs];
		psi[idx] = new int[this->num_labs];
		for (size_t clab=0; clab<this->num_labs; clab++) {
			delta[idx][clab]=0;
			psi[idx][clab]=0;
		}
	}

	QNUInt32 idx=0;

	// Initialization
	for (size_t clab=0; clab<this->num_labs; clab++) {
		//delta[0][j]=log(R[j]);
		//delta[0][j]=Mstate[j];
		delta[0][clab]=this->nodeList->at(idx)->getStateValue(clab);
		psi[0][clab]=-1;
	}

	// Recursion
	//int idx=1;
	idx=1;
	while (idx <= nodeCnt+1) {
	//while (idx < seq_len) {
		for (size_t clab=0; clab<this->num_labs; clab++) {
			for (size_t plab=0; plab<this->num_labs; plab++) {
				float dval=this->nodeList->at(idx-1)->getFullTransValue(plab,clab);
				if (delta[idx][clab] < delta[idx-1][plab]+dval) {
					delta[idx][clab] = delta[idx-1][plab]+dval;
					psi[idx][clab] = plab;
				}
			}
		}
		//cur_seq=cur_seq->getNext();
		idx++;
	}

	CRF_LabelPath* bestPath = new CRF_LabelPath(seq_len);
	QNUInt32* bp = bestPath->getLabelPath();

	int best=0;

	float bestValue=delta[seq_len][0];


	for (size_t clab=0; clab<this->num_labs; clab++) {
		//cout << "Label: " << j << " Delta value " << delta[seq_len][j] << " Best: " << bestValue << endl;
		if (delta[seq_len][clab]>bestValue) {
			bestValue=delta[seq_len][clab];
			best=clab;
		}
	}



	bp[seq_len-1]=best;
	//cout << "Seq Len: " << seq_len << " Best: " << best << endl;

	for (int t=seq_len-2; t>=0; t--) {
		//cout << "Fine in loop: " << t << " BP[t+1]: " << bp[t+1] << " PSI: " << psi[t+1][bp[t+1]] << endl;
		bp[t]=psi[t+1][bp[t+1]]; // Assign to best path the backpointer from the successive best path
	}

	for (int idx=0; idx<seq_len+1; idx++) {
		delete[] delta[idx];
		delete[] psi[idx];
	}

/*	cur_seq=seq_head;
	while (cur_seq != NULL) {
		last_seq=cur_seq;
		cur_seq=cur_seq->getNext();
		delete last_seq;
	}*/

	delete[] delta;
	delete[] psi;

	return bestPath;

}

CRF_LabelPath* CRF_NewViterbi::bestPathVec()
{
	//if (this->crf->getFeatureMap()->getNumStates() != 1) {
		return this->nStateBestPathVec();
	//}

	vector <vector <int> > psi;
	vector <vector <float> > delta;

	QNUInt32 ftr_count;

	int seq_len=0;
	QNUInt32 nodeCnt=0;
	//cout << "Starting Viterbi Best Path: " << endl;
	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		//cout << "Reading frame " << nodeCnt << endl;
		ftr_count=ftr_strm->read(this->bunch_size,ftr_buf,lab_buf);
		//cout << "Frame read " << nodeCnt << " " << ftr_count << endl;

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
			this->nodeList->set(nodeCnt,new_buf,num_ftrs,this->lab_buf[i],this->crf);
			float value=this->nodeList->at(nodeCnt)->computeTransMatrix();

			psi.push_back(vector <int> (this->num_labs,0));
			delta.push_back(vector <float> (this->num_labs,0));
			//cout << "Delta Psi Handling" << endl;
			if (nodeCnt == 0) {
				// First value handled somewhat differently
				for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
					delta[nodeCnt][clab]=this->nodeList->at(nodeCnt)->getStateValue(clab);
					psi[nodeCnt][clab]=-1;
				}
			}
			else {
				// Normal handling
				for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
					for (QNUInt32 plab=0; plab<this->num_labs; plab++) {
						float dval=this->nodeList->at(nodeCnt-1)->getFullTransValue(plab,clab);
						//float dval=0;
						if (delta[nodeCnt][clab] < delta[nodeCnt-1][plab]+dval) {
							delta[nodeCnt][clab] = delta[nodeCnt-1][plab]+dval;
							psi[nodeCnt][clab] = plab;
						}
					}
				}
			}
			//cout << "Delta Psi Handled" << endl;
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	psi.push_back(vector <int> (this->num_labs,0));
	delta.push_back(vector <float> (this->num_labs,0));
	for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
		for (QNUInt32 plab=0; plab<this->num_labs; plab++) {
			float dval=this->nodeList->at(nodeCnt-1)->getFullTransValue(plab,clab);
			//float dval=0;
			if (delta[nodeCnt][clab] < delta[nodeCnt-1][plab]+dval) {
				delta[nodeCnt][clab] = delta[nodeCnt-1][plab]+dval;
				psi[nodeCnt][clab] = plab;
			}
		}
	}
	//cout << "Delta & Psi computed " << endl;
	//cout << "Finding Best Path" << endl;

	CRF_LabelPath* bestPath = new CRF_LabelPath(nodeCnt);
	QNUInt32* bp = bestPath->getLabelPath();

	int best=0;

	float bestValue=delta[nodeCnt][0];

	for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
		if (delta[nodeCnt][clab]>bestValue) {
			bestValue=delta[nodeCnt][clab];
			best=clab;
		}
	}


	bp[nodeCnt-1]=best;

	for (int t=nodeCnt-2; t>=0; t--) {
		bp[t]=psi[t+1][bp[t+1]]; // Assign to best path the backpointer from the successive best path
	}
	//cout << "Best Path Found" << endl;
	return bestPath;

}


CRF_LabelPath* CRF_NewViterbi::nStateBestPathVec()
{

	QNUInt32 nStates = this->crf->getFeatureMap()->getNumStates();

	vector <vector <int> > psi;
	vector <vector <float> > delta;

	QNUInt32 ftr_count;

	int seq_len=0;
	QNUInt32 nodeCnt=0;
	//cout << "Starting Viterbi N-State Best Path: " << endl;
	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		//cout << "Reading frame " << nodeCnt << endl;
		ftr_count=ftr_strm->read(this->bunch_size,ftr_buf,lab_buf);
		//cout << "Frame read " << nodeCnt << " " << ftr_count << endl;

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
			this->nodeList->set(nodeCnt,new_buf,num_ftrs,this->lab_buf[i],this->crf);
			float value=this->nodeList->at(nodeCnt)->computeTransMatrix();

			psi.push_back(vector <int> (this->num_labs,LOG0));
			delta.push_back(vector <float> (this->num_labs,LOG0));
			if (nodeCnt == 0) {
				// First value handled somewhat differently
				// For the start states we set the value normally, for the non-start states leave it
				//  as LOG0
				for (QNUInt32 clab=0; clab<this->num_labs; clab+=nStates) {
					delta[nodeCnt][clab]=this->nodeList->at(nodeCnt)->getStateValue(clab);
					psi[nodeCnt][clab]=-1;
				}
			}
			else {
				// Normal handling
				for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
					// First check for clab as a start-state.  If so, handle end-state to start state trans
					if (clab % nStates == 0) {
						// Note that if nStates == 1 then the above is always true and the below
						//  is computed for all possible transitions
						for (QNUInt32 plab=nStates-1; plab<this->num_labs; plab+=nStates) {
							float dval=this->nodeList->at(nodeCnt-1)->getFullTransValue(plab,clab);
							//float dval=0;
							if (delta[nodeCnt][clab] < delta[nodeCnt-1][plab]+dval) {
								delta[nodeCnt][clab] = delta[nodeCnt-1][plab]+dval;
								psi[nodeCnt][clab] = plab;
							}
						}
						// Special case handling - if nStates == 1, then the above loop handles self
						//   transitions.  But for nStates>1 it does not.
						if (nStates>1) {
							QNUInt32 plab=clab;
							float dval=this->nodeList->at(nodeCnt-1)->getFullTransValue(plab,clab);
							if (delta[nodeCnt][clab] < delta[nodeCnt-1][plab]+dval) {
								delta[nodeCnt][clab] = delta[nodeCnt-1][plab]+dval;
								psi[nodeCnt][clab] = plab;
							}
						}
					}
					else {
						// Otherwise handle it as an internal state - we only care about transitions from
						// a previous label and self transitions.
						// Previous label transitions
						QNUInt32 plab=clab-1;
						float dval=this->nodeList->at(nodeCnt-1)->getFullTransValue(plab,clab);
						if (delta[nodeCnt][clab] < delta[nodeCnt-1][plab]+dval) {
							delta[nodeCnt][clab] = delta[nodeCnt-1][plab]+dval;
							psi[nodeCnt][clab] = plab;
						}
						// Self transitions
						plab=clab;
						dval=this->nodeList->at(nodeCnt-1)->getFullTransValue(plab,clab);
						if (delta[nodeCnt][clab] < delta[nodeCnt-1][plab]+dval) {
							delta[nodeCnt][clab] = delta[nodeCnt-1][plab]+dval;
							psi[nodeCnt][clab] = plab;
						}
					}
				}
			}
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);
	psi.push_back(vector <int> (this->num_labs,0));
	delta.push_back(vector <float> (this->num_labs,0));
	for (QNUInt32 clab=0; clab<this->num_labs; clab++) {
		for (QNUInt32 plab=0; plab<this->num_labs; plab++) {
			float dval=this->nodeList->at(nodeCnt-1)->getFullTransValue(plab,clab);
			//float dval=0;
			if (delta[nodeCnt][clab] < delta[nodeCnt-1][plab]+dval) {
				delta[nodeCnt][clab] = delta[nodeCnt-1][plab]+dval;
				psi[nodeCnt][clab] = plab;
			}
		}
	}

	//cout << "Delta & Psi computed " << endl;
	//cout << "Finding Best Path" << endl;

	CRF_LabelPath* bestPath = new CRF_LabelPath(nodeCnt);
	QNUInt32* bp = bestPath->getLabelPath();

	int best=0;

	// We only need to check for the best value of an end-state - we've got to stop at an end state
	//float bestValue=delta[nodeCnt][0];
	float bestValue=delta[nodeCnt][nStates-1];

	for (QNUInt32 clab=nStates-1; clab<this->num_labs; clab=clab+nStates) {
		if (delta[nodeCnt][clab]>bestValue) {
			bestValue=delta[nodeCnt][clab];
			best=clab;
		}
	}


	bp[nodeCnt-1]=best;

	for (int t=nodeCnt-2; t>=0; t--) {
		bp[t]=psi[t+1][bp[t+1]]; // Assign to best path the backpointer from the successive best path
	}
	//cout << "Best Path Found" << endl;
	return bestPath;

}


CRF_LabelPath* CRF_NewViterbi::alignPath()
{
	// Returns the best path through the current segment, aligned to the hardtarget label sequence

	size_t ftr_count;


	QNUInt32 seq_len=0;
	QNUInt32 nodeCnt=0;
	//cout << "Starting Viterbi Best Path: " << endl;

	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		ftr_count=ftr_strm->read(this->bunch_size,this->ftr_buf,this->lab_buf);

		//cout << "Feature Count: " << ftr_count << endl;
		for (QNUInt32 i=0; i<ftr_count; i++) {
			//cout << "\tLabel: " << lab_buf[i] << "\tFeas:";
			// Now, separate the bunch into individual frames

			float* new_buf = new float[this->num_ftrs];
			for (QNUInt32 j=0; j<this->num_ftrs; j++) {
				int idx=i*this->num_ftrs+j;
				new_buf[j]=this->ftr_buf[idx];
				//cout << " " << new_buf[j];
			}
			//cout << endl;
			// Store the current frame/label information in a sequence node
			//	* sequence nodes create a doubly-linked list, with the previous node known at creation time
			//next_seq = new CRF_Seq(new_buf,this->num_ftrs,this->lab_buf[i],crf->getNLabs(),cur_seq);
			this->nodeList->set(nodeCnt,new_buf,num_ftrs,this->lab_buf[i],this->crf);
			seq_len++;
			/*if (seq_head == NULL)
			{
				seq_head=next_seq;
			}
			if (cur_seq != NULL) {
				cur_seq->setNext(next_seq);
			}
			cur_seq=next_seq;
			//float value=this->gb->computeTransMatrix(cur_seq,crf);
			float value=this->gb->computeTransMatrixLog(cur_seq);*/
			float value=this->nodeList->at(nodeCnt)->computeTransMatrix();
			nodeCnt++;
		}
	} while (ftr_count >= this->bunch_size);

	nodeCnt--; //Correct for adding 1 at the end of the above sequence

	// Alpha values have been computed and we know how long the sequence is.
	// Next, we need to get the underlying label sequence.  Per-frame labels are stored in the sequence nodes,
	// so we need to collapse this into a higher-level sequence.

	QNUInt32* labPath = new QNUInt32[seq_len]; // Maximum size is the length of the sequence

	//cur_seq=seq_head;
	//labPath[0]=cur_seq->getLab();

	QNUInt32 cur_idx=0;
	labPath[cur_idx]=this->nodeList->at(cur_idx)->getLabel();
	//cur_seq=cur_seq->getNext();

	for (QNUInt32 idx=1; idx<seq_len; idx++) {
		//QNUInt32 cur_lab=cur_seq->getLab();
		QNUInt32 cur_lab=this->nodeList->at(idx)->getLabel();
		if (cur_lab != labPath[cur_idx]) {
			// Labels are different, increment and store
			cur_idx++;
			labPath[cur_idx]=cur_lab;
		}
		//cur_seq=cur_seq->getNext();
	}

	QNUInt32 tot_labs=cur_idx+1;

	// Now we can perform the Viterbi algorithm over the stored values
	typedef float* floatPtr;
	typedef double* doublePtr;
	typedef int* intPtr;

	doublePtr* delta = new doublePtr[seq_len+1];
	doublePtr* hold = new doublePtr[seq_len+1];
	intPtr* psi = new intPtr[seq_len+1];

	for (QNUInt32 idx=0; idx<seq_len+1; idx++) {
		delta[idx] = new double[tot_labs];
		psi[idx] = new int[tot_labs];
		hold[idx] = new double[tot_labs];
		for (size_t clab=0; clab<tot_labs; clab++) {
			delta[idx][clab]=0;
			psi[idx][clab]=0;
			hold[idx][clab]=0;
		}
	}

/*	cur_seq=seq_head;
	alpha=cur_seq->getAlpha();
	Mtrans=cur_seq->getMtrans();
	Mstate=cur_seq->getMstate();*/


	//delta[0][0]=Mstate[labPath[0]];
	delta[0][0]=this->nodeList->at(0)->getStateValue(labPath[0]);
	psi[0][0]=-1;
	hold[0][0]=delta[0][0];

	// Recursion
	QNUInt32 idx=1;
	//while (cur_seq != NULL) {
	while (idx <= nodeCnt+1) {
		/*alpha=cur_seq->getAlpha();
		Mtrans=cur_seq->getMtrans();
		Mstate=cur_seq->getMstate();*/
		QNUInt32 maxLabel;
		// We want to limit this so that we start at the upper corner
		// and work downwards
		if (idx < tot_labs) {
			maxLabel=idx;
		}
		else {
			maxLabel=tot_labs;
		}
		for (size_t j=0; j<maxLabel; j++) {
			// Each Label can make only one of two transitions - either a transition to itself OR a transition to
			// the next label in the sequence.

			// First compute the self transition
			QNUInt32 prev=j;
			QNUInt32 prevIdx = labPath[prev];
			//float dval=Mstate[prevIdx]+Mtrans[prevIdx*this->num_labs+prevIdx];
			float dval=this->nodeList->at(idx-1)->getFullTransValue(prevIdx,prevIdx);
			if (delta[idx][prev] < delta[idx-1][prev]+dval) {
				//cout << "DVAL: " << dval << endl;
				delta[idx][prev] = delta[idx-1][prev]+dval;
				psi[idx][prev] = prev;
				hold[idx][prev] = dval;
			}

			// Next, if we are not in the last state we can make a transition to the next state
			// if it's better than the self transition
			if (j < (maxLabel-1)) {
				QNUInt32 next=j+1;
				QNUInt32 nextIdx=labPath[next];
				//dval = Mstate[nextIdx]+Mtrans[prevIdx*this->num_labs+nextIdx];
				dval=this->nodeList->at(idx-1)->getFullTransValue(prevIdx,nextIdx);

				if (delta[idx][next] < delta[idx-1][prev]+dval) {
					delta[idx][next] = delta[idx-1][prev]+dval;
					hold[idx][next] = dval;
					psi[idx][next]=prev;
				}
			}
		}
		//cur_seq=cur_seq->getNext();
		idx++;
	}


	CRF_LabelPath* bestPath = new CRF_LabelPath(seq_len);
	QNUInt32* bp = bestPath->getLabelPath();

	CRF_LabelPath* alignPath = new CRF_LabelPath(seq_len);
	QNUInt32* ap = alignPath->getLabelPath();

	// Unlike standard Viterbi, we know which state we have to end in
	//  start from the backpointer there and see the best path that
	//  leads us there.

	ap[seq_len-1]=tot_labs-1; //Index of our terminal position is the last label in the sequence
	bp[seq_len-1]=labPath[tot_labs-1];  // bestPath actually stores the associated label number

	for (int t=seq_len-2; t>=0; t--) {
		//cout << "Fine in loop: " << t << " BP[t+1]: " << bp[t+1] << " PSI: " << psi[t+1][ap[t+1]];
		ap[t]=psi[t+1][ap[t+1]]; // Assign to best aligned path the backpointer from the successive best aligned path
		//cout << " AP[t+1]: " << ap[t+1] << " AP[t]: " << ap[t] << endl;
		// Check to make sure that we are either making a self transition OR moving back one step
		//  * This shouldn't be necessary - in the loop above we're restricting to these two transitions, but
		//    just in case we're wrong, we'll put some error checking here
		if ((ap[t] != ap[t+1]) && (ap[t] != ap[t+1]-1)) {
			cerr << "ERROR: Skipped Label - moved from label at pos " << ap[t+1] << " to label at pos " << ap[t] << endl;
			for (int frm=t; frm<=seq_len; frm++) {
				cerr << frm << "\t" << ap[frm] << "\t" << labPath[ap[frm]] << "\n";
			}
			exit(-1);
		}
		bp[t]=labPath[ap[t]];  // Assign to best path the associated label
	}

	if (ap[0] != 0) {
		// Again, this check just makes sure that we ended at the starting point.  This shouldn't be necessary
		// either, but just in case there's a bug in the logic above
		cerr << "ERROR: Label at timestep 0 is label at position " << ap[0] << " instead of label at position 0" << endl;
		exit(-1);
	}

	for (QNUInt32 idx=0; idx<seq_len+1; idx++) {
		delete[] delta[idx];
		delete[] psi[idx];
		delete[] hold[idx];
	}

	/*cur_seq=seq_head;
	while (cur_seq != NULL) {
		last_seq=cur_seq;
		cur_seq=cur_seq->getNext();
		delete last_seq;
	}*/

	delete[] delta;
	delete[] psi;
	delete[] hold;
	delete[] labPath;
	delete alignPath;

	return bestPath;
}
