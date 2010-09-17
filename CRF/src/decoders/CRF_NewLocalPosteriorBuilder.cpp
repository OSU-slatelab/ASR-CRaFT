/*
 * CRF_NewLocalPosteriorBuilder.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_NewLocalPosteriorBuilder.h"



/*
 * CRF_NewLoclPosteriorBuilder constructor
 *
 * Input: *crf_in - pointer to the CRF model to be used for building lattices
 * 	      norm - boolean flag to control normalization of output posteriors (defaults to true)
 *
 */
CRF_NewLocalPosteriorBuilder::CRF_NewLocalPosteriorBuilder(CRF_Model* crf_in, bool norm)
	: crf(crf_in),
	  normalize(norm)
{
	//this->nodeList=new CRF_StdStateVectorLog();
	this->ftr_buf=NULL;
	this->lab_buf=NULL;
	this->alpha_base = new double[this->crf->getNLabs()];
	for (QNUInt32 i=0; i<this->crf->getNLabs(); i++) {
		this->alpha_base[i]=0.0;
	}
	this->nodeList = new CRF_StateVector();
	this->ownsNodeList=true;
}

/*
 * CRF_NewLoclPosteriorBuilder destructor
 *
 */

CRF_NewLocalPosteriorBuilder::~CRF_NewLocalPosteriorBuilder()
{
	if (this->ftr_buf != NULL) { delete [] ftr_buf;}
	if (this->lab_buf != NULL) { delete [] lab_buf;}
	if (this->alpha_base != NULL) { delete [] alpha_base;}
	// should delete nodeList?
}

/*
 * CRF_NewLocalPosteriorBuilder::buildFtrSeq
 *
 * Input: *ftr_strm - pointer to feature stream of input features
 *
 * Returns: vector of nodes containing local posterior values
 *
 * Reads the current sequence of feature observations from the input feature
 * stream and uses them to construct a vector of nodes of local posterior values
 * calculated using the CRF Model used to instantiate the PosteriorBuilder object.
 * Posterior values are computed using the forward-backward algorithm in a manner
 * similar to that used for training the CRF via graident descent.
 *
 */
CRF_StateVector* CRF_NewLocalPosteriorBuilder::buildFtrSeq(CRF_FeatureStream* ftr_strm)
{
	size_t bunch_size = 3;
	size_t num_ftrs=ftr_strm->num_ftrs();

	if (this->ftr_buf==NULL) { // First pass through, initialize the buffers
		this->ftr_buf = new float[num_ftrs*bunch_size];
		this->lab_buf = new QNUInt32[bunch_size];
	}

	//CRF_StateVector* nodeList = new CRF_StdStateVectorLog();

	size_t ftr_count;

	int seq_cnt=0;
	QNUInt32 nodeCnt=0;
	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		ftr_count=ftr_strm->read(bunch_size,ftr_buf,lab_buf);

		for (QNUInt32 i=0; i<ftr_count; i++) {
			// Now, separate the bunch into individual frames
			float* new_buf = new float[num_ftrs];
			for (QNUInt32 j=0; j<num_ftrs; j++) {
				int idx=i*num_ftrs+j;
				new_buf[j]=ftr_buf[idx];
				//cout << " " << new_buf[j];
			}
			//cout << endl;
			// Store the current frame/label information in a sequence node
			//	* sequence nodes create a doubly-linked list, with the previous node known at creation time
			//  * new_buf will be deleted when this sequence object gets deleted
			//next_seq = new CRF_Seq(new_buf,num_ftrs,lab_buf[i],this->num_labs,cur_seq);
			nodeList->set(nodeCnt,new_buf,lab_buf[i],this->crf->getNLabs(),this->crf);

			//double value=this->computeTransMatrix(cur_seq,crf);
			//double value=this->computeTransMatrixLog(cur_seq); // logspace computation
			double value=nodeList->at(nodeCnt)->computeTransMatrix();
			double scale;
			double* prev_alpha;
			if (nodeCnt == 0) {
				prev_alpha=this->alpha_base;
				scale=nodeList->at(nodeCnt)->computeFirstAlpha(prev_alpha);
			}
			else {
				prev_alpha=nodeList->at(nodeCnt-1)->getAlpha();
				scale=nodeList->at(nodeCnt)->computeAlpha(prev_alpha);
			}
			//scale=nodeList->at(nodeCnt)->computeAlpha(prev_alpha);
			nodeCnt++;
			seq_cnt++;
		}

	} while (ftr_count >= bunch_size);

	nodeList->setNodeCount(nodeCnt);
	nodeCnt--;//Correct for the fact that we add 1 to the nodeCnt at the end of the above loop...
			  // So the last index is actually nodeCount-1

	QNUInt32 lastNode=nodeCnt;
	double Zx=nodeList->at(lastNode)->computeAlphaSum();

	double norm_const;

	if (this->normalize) {
		norm_const=Zx;
	}
	else {
		norm_const=0.0;
	}

	bool stop=false;

	int s_counter=0;
	//while (cur_seq != NULL) {
	while (!stop) {

		double* beta = nodeList->at(nodeCnt)->getBeta();
		if (nodeCnt==lastNode) {
			nodeList->at(nodeCnt)->setTailBeta();
		}
		else {
			// We compute the beta value for the node following our current one, and store the result
			// as the beta for our current node (as per the equations).
			nodeList->at(nodeCnt+1)->computeBeta(beta,nodeList->at(nodeCnt)->getAlphaScale());
		}
		//cout << "Normalizing constant Zx: " << Zx << "\t" << norm_const << endl;
		double* alpha_beta = nodeList->at(nodeCnt)->computeAlphaBeta(norm_const);
		double alpha_beta_tot=0.0;
		for (QNUInt32 clab=0; clab<this->crf->getNLabs(); clab++) {
			try {
				if (this->normalize) {
					//double ab= alpha_beta[clab];
					//cout << "AB: " << ab << endl;
				alpha_beta_tot += expE(alpha_beta[clab]);
				}
				else {
					double ab = alpha_beta[clab]-Zx;
					//cout << "AB: " << ab << endl;
					alpha_beta_tot += expE(ab);
				}
			}
			catch (overflow_error &e) {
				string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq caught exception: "+string(e.what())+" while computing alpha_beta";
				throw runtime_error(errstr);
				return(NULL);
			}
		}
		//cout << "Alpha Beta Tot: " << alpha_beta_tot << endl;
		if ((alpha_beta_tot >1.1))  {
			cout << "Total: " << alpha_beta_tot << endl;
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq: Probability sums greater than 1.0";
			throw runtime_error(errstr);
			return(NULL);
		}
		if (alpha_beta_tot < 0.9) {
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq: Probability sums less than 1.0";
			throw runtime_error(errstr);
			return(NULL);
		}
		if (nodeCnt==0) { stop=true; } // unconventional while loop driver due to unsigned int type
		nodeCnt--;
		s_counter++;
	}

	return nodeList;
}

/*
 * CRF_NewLocalPosteriorBuilder::buildFtrSeqNState
 *
 * Input: *ftr_strm - pointer to feature stream of input features
 *
 * Returns: vector of nodes containing local posterior values
 *
 * As above, but builds the local posterior vector using multi-state constraints
 * in the alpha-beta computation.
 *
 */
CRF_StateVector* CRF_NewLocalPosteriorBuilder::buildFtrSeqNState(CRF_FeatureStream* ftr_strm)
{
	size_t bunch_size = 3;
	size_t num_ftrs=ftr_strm->num_ftrs();

	//QNUInt32 nStates = this->crf->getFeatureMap()->getNumStates();

	if (this->ftr_buf==NULL) { // First pass through, initialize the buffers
		this->ftr_buf = new float[num_ftrs*bunch_size];
		this->lab_buf = new QNUInt32[bunch_size];
	}

	//CRF_StateVector* nodeList = new CRF_StdStateVectorLog();

	size_t ftr_count;

	int seq_cnt=0;
	QNUInt32 nodeCnt=0;
	do {
		// First, read in the next training value from the file
		//	We can read in a "bunch" at a time, then separate them into individual frames
		ftr_count=ftr_strm->read(bunch_size,ftr_buf,lab_buf);

		for (QNUInt32 i=0; i<ftr_count; i++) {
			// Now, separate the bunch into individual frames
			float* new_buf = new float[num_ftrs];
			for (QNUInt32 j=0; j<num_ftrs; j++) {
				int idx=i*num_ftrs+j;
				new_buf[j]=ftr_buf[idx];
				//cout << " " << new_buf[j];
			}
			//cout << endl;
			// Store the current frame/label information in a sequence node
			//	* sequence nodes create a doubly-linked list, with the previous node known at creation time
			//  * new_buf will be deleted when this sequence object gets deleted
			//next_seq = new CRF_Seq(new_buf,num_ftrs,lab_buf[i],this->num_labs,cur_seq);
			nodeList->set(nodeCnt,new_buf,lab_buf[i],this->crf->getNLabs(),this->crf);

			//double value=this->computeTransMatrix(cur_seq,crf);
			//double value=this->computeTransMatrixLog(cur_seq); // logspace computation
			double value=nodeList->at(nodeCnt)->computeTransMatrix();
			double scale;
			double* prev_alpha;
			if (nodeCnt == 0) {
				prev_alpha=this->alpha_base;
			}
			else {
				prev_alpha=nodeList->at(nodeCnt-1)->getAlpha();
			}
			scale=nodeList->at(nodeCnt)->computeAlpha(prev_alpha);
			nodeCnt++;
			seq_cnt++;
		}

	} while (ftr_count >= bunch_size);

	nodeList->setNodeCount(nodeCnt);
	nodeCnt--;//Correct for the fact that we add 1 to the nodeCnt at the end of the above loop...
			  // So the last index is actually nodeCount-1

	QNUInt32 lastNode=nodeCnt;
	double Zx=nodeList->at(lastNode)->computeAlphaSum();

	bool stop=false;

	int s_counter=0;
	//while (cur_seq != NULL) {
	while (!stop) {

		double* beta = nodeList->at(nodeCnt)->getBeta();
		if (nodeCnt==lastNode) {
			nodeList->at(nodeCnt)->setTailBeta();
		}
		else {
			// We compute the beta value for the node following our current one, and store the result
			// as the beta for our current node (as per the equations).
			nodeList->at(nodeCnt+1)->computeBeta(beta,nodeList->at(nodeCnt)->getAlphaScale());
		}
		double* alpha_beta = nodeList->at(nodeCnt)->computeAlphaBeta(Zx);
		double alpha_beta_tot=0.0;
		for (QNUInt32 clab=0; clab<this->crf->getNLabs(); clab++) {
			try {
				alpha_beta_tot += expE(alpha_beta[clab]);
			}
			catch (overflow_error &e) {
				string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq caught exception: "+string(e.what())+" while computing alpha_beta";
				throw runtime_error(errstr);
				return(NULL);
			}
		}
		//cout << "Alpha Beta Tot: " << alpha_beta_tot << endl;
		if ((alpha_beta_tot >1.1))  {
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq: Probability sums greater than 1.0";
			throw runtime_error(errstr);
			return(NULL);
		}
		if (alpha_beta_tot < 0.9) {
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq: Probability sums less than 1.0";
			throw runtime_error(errstr);
			return(NULL);
		}
		if (nodeCnt==0) { stop=true; } // unconventional while loop driver due to unsigned int type
		nodeCnt--;
		s_counter++;
	}

	return nodeList;
}


/*
void CRF_NewLocalPosteriorBuilder::computeAlphaBeta(CRF_StateVector *nodelist)
{
	if (this->ownsNodeList) {
		delete this->nodeList;
	}
	this->nodeList=nodelist;
	this->ownsNodeList=false;

	for (int nodeCnt=0;nodeCnt<nodeList->getNodeCount();nodeCnt++){
		double value=nodeList->at(nodeCnt)->computeTransMatrix();
		double scale;
		double* prev_alpha;
		if (nodeCnt == 0) {
			prev_alpha=this->alpha_base;
			scale=nodeList->at(nodeCnt)->computeFirstAlpha(prev_alpha);
		}
		else {
			prev_alpha=nodeList->at(nodeCnt-1)->getAlpha();
			scale=nodeList->at(nodeCnt)->computeAlpha(prev_alpha);
		}

	}


	QNUInt32 lastNode=nodeList->getNodeCount()-1;
	double Zx=nodeList->at(lastNode)->computeAlphaSum();

	double norm_const;

	if (this->normalize) {
		norm_const=Zx;
	}
	else {
		norm_const=0.0;
	}

	bool stop=false;

	int s_counter=0;
	//while (cur_seq != NULL) {
	int nodeCnt=lastNode;
	while (!stop) {

		double* beta = nodeList->at(nodeCnt)->getBeta();
		if (nodeCnt==lastNode) {
			nodeList->at(nodeCnt)->setTailBeta();
		}
		else {
			// We compute the beta value for the node following our current one, and store the result
			// as the beta for our current node (as per the equations).
			nodeList->at(nodeCnt+1)->computeBeta(beta,nodeList->at(nodeCnt)->getAlphaScale());
		}
		//cout << "Normalizing constant Zx: " << Zx << "\t" << norm_const << endl;
		double* alpha_beta = nodeList->at(nodeCnt)->computeAlphaBeta(norm_const);
		double alpha_beta_tot=0.0;
		for (QNUInt32 clab=0; clab<this->crf->getNLabs(); clab++) {
			try {
				if (this->normalize) {
					double ab= alpha_beta[clab];
					//cout << "AB: " << ab << endl;
				alpha_beta_tot += expE(alpha_beta[clab]);
				}
				else {
					double ab = alpha_beta[clab]-Zx;
					//cout << "AB: " << ab << endl;
					alpha_beta_tot += expE(ab);
				}
			}
			catch (overflow_error &e) {
				string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq caught exception: "+string(e.what())+" while computing alpha_beta";
				throw runtime_error(errstr);
				//return(NULL);
			}
		}
		//cout << "Alpha Beta Tot: " << alpha_beta_tot << endl;
		if ((alpha_beta_tot >1.1))  {
			cout << "Total: " << alpha_beta_tot << endl;
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq: Probability sums greater than 1.0";
			throw runtime_error(errstr);
			//return(NULL);
		}
		if (alpha_beta_tot < 0.9) {
			string errstr="CRF_LocalPosteriorBuilder::buildFtrSeq: Probability sums less than 1.0";
			throw runtime_error(errstr);
			//return(NULL);
		}
		if (nodeCnt==0) { stop=true; } // unconventional while loop driver due to unsigned int type
		nodeCnt--;
		s_counter++;
	}

}

*/
