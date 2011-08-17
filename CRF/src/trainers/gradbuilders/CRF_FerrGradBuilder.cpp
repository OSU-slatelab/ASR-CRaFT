/*
 * CRF_FerrGradBuilder.cpp
 *
 *  Created on: Oct 2, 2009
 *      Author: morrijer
 */

#include "CRF_FerrGradBuilder.h"
#include "../../nodes/CRF_StdStateNode.h"
#include "../../nodes/CRF_StateVector.h"
#include "fst/fstlib.h"

typedef StdArc::StateId StateId;


CRF_FerrGradBuilder::CRF_FerrGradBuilder(CRF_Model* crf_in)
: CRF_GradBuilder(crf_in){
	// TODO Auto-generated constructor stub
	this->lambda_len = crf_in->getLambdaLen();
	this->ExpF = new double[this->lambda_len];
	this->lb = NULL;
	this->nodeList = NULL;
}

CRF_FerrGradBuilder::~CRF_FerrGradBuilder() {
	// TODO Auto-generated destructor stub
}

double CRF_FerrGradBuilder::buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out)
{
	double logli=0.0;
	int logflag=(int) *Zx_out;
	//cerr << "Entering BuildGradient" << endl;
	if (this->lb == NULL ) {
		//cerr << "Creating LatticeBuilder" << endl;
		// Slight hack so that we don't need to re-engineer LatticeBuilder
		this->lb = new CRF_LatticeBuilder(ftr_strm,this->crf);
		this->nodeList=this->lb->getNodeList();
	}
	QNUInt32 lambda_len = this->crf->getLambdaLen();
	double* lambda = this->crf->getLambda();
	for (QNUInt32 i=0; i<lambda_len; i++) {
		this->ExpF[i]=0.0;
	}

	VectorFst<StdArc>* phn_lat=new StdVectorFst();
	VectorFst<StdArc>* lab_lat=new StdVectorFst();

	//cerr << "Building lattice" << endl;
	if (this->crf->getFeatureMap()->getNumStates()==1) {
		lb->buildLattice(phn_lat,true,lab_lat);
	}
	else {
		lb->nStateBuildLattice(phn_lat,true,lab_lat);
	}
	//cerr << "Lattice Built" << endl;
	VectorFst<StdArc>* shortest_fst = new VectorFst<StdArc>();
	ShortestPath(*phn_lat,shortest_fst,1);
	Project(shortest_fst,PROJECT_OUTPUT);
	RmEpsilon(shortest_fst);
	TopSort(shortest_fst);

	QNUInt32 num_labels = shortest_fst->NumStates() - 1;
	QNUInt32* bestarr = new QNUInt32[num_labels];
	//Gets the initial state; if kNoState => empty FST.
	StateId initial_state = shortest_fst->Start();
	//	Iterates over the FSTs states.
	QNUInt32 frm_no=0;
	for (StateIterator<StdFst> siter(*shortest_fst); !siter.Done(); siter.Next())
	{
			StateId state_id = siter.Value();
			//Iterates over state state_id's arcs.
			for (ArcIterator<StdFst> aiter(*shortest_fst, state_id); !aiter.Done(); aiter.Next())
		{
				const StdArc &arc = aiter.Value();
				bestarr[frm_no]=arc.olabel-1;
				frm_no++;
		}
	}

	// labarr now contains the best path through the utterance based on the current model

	delete shortest_fst;
	shortest_fst = new VectorFst<StdArc>();
	ComposeFst<StdArc>* result=new ComposeFst<StdArc>(*phn_lat,*lab_lat);
	ShortestPath(*result,shortest_fst,1);
	delete result;
	Project(shortest_fst,PROJECT_OUTPUT);
	RmEpsilon(shortest_fst);
	TopSort(shortest_fst);

	QNUInt32 true_num_labels = shortest_fst->NumStates() - 1;
	//vector<StdArc::Weight>* distance = new vector<StdArc::Weight>[true_num_labels];
	vector<StdArc::Weight> distance[true_num_labels];
	ShortestDistance(*shortest_fst,distance);

	logli=distance->at(true_num_labels).Value();

	//delete distance;

	QNUInt32* truearr = new QNUInt32[true_num_labels];
	//Gets the initial state; if kNoState => empty FST.
	initial_state = shortest_fst->Start();
	//	Iterates over the FSTs states.
	frm_no=0;
	for (StateIterator<StdFst> siter(*shortest_fst); !siter.Done(); siter.Next())
	{
			StateId state_id = siter.Value();
			//Iterates over state state_id's arcs.
			for (ArcIterator<StdFst> aiter(*shortest_fst, state_id); !aiter.Done(); aiter.Next())
		{
				const StdArc &arc = aiter.Value();
				truearr[frm_no]=arc.olabel-1;
				frm_no++;
		}
	}

	// true arr now contains the force aligned best path through our current model

	if (true_num_labels != num_labels) {
		string errstr="CRF_FerrGradBuilder had mismatched bestarr and truearr lengths";
		throw runtime_error(errstr);
	}
	QNUInt32 nodeCnt=this->nodeList->getNodeCount();

	QNUInt32 plab=0;
	for (QNUInt32 i=0; i<num_labels; i++) {
		float* ftr_buf = this->nodeList->at(i)->getFtrBuffer();
		QNUInt32 tlab=this->nodeList->at(i)->getLabel();
		//if (counter % 10 == 0) {
		//cerr << num_labels << " " << i << " TrueArr[i]: " << truearr[i] << "BestArr[i]: " << bestarr[i] << endl;
		//}
		this->crf->getFeatureMap()->computeStateExpF(ftr_buf,lambda,ExpF,NULL,1,truearr[i]+1,truearr[i]);
		//this->crf->getFeatureMap()->computeStateExpF(ftr_buf,lambda,ExpF,NULL,1,tlab+1,tlab);
		this->crf->getFeatureMap()->computeStateExpF(ftr_buf,lambda,ExpF,NULL,-1,bestarr[i]+1,bestarr[i]);
		if (i!=0) {
			this->crf->getFeatureMap()->computeTransExpF(ftr_buf,lambda,ExpF,NULL,1,truearr[i-1]+1,truearr[i]+1,truearr[i-1],truearr[i]);
			//this->crf->getFeatureMap()->computeTransExpF(ftr_buf,lambda,ExpF,NULL,1,plab+1,tlab+1,plab,tlab);
			this->crf->getFeatureMap()->computeTransExpF(ftr_buf,lambda,ExpF,NULL,-1,bestarr[i-1]+1,bestarr[i]+1,bestarr[i-1],bestarr[i]);
		}
		plab=tlab;

	}

	//cerr << "Computing ExpF" << endl;
	double max_val=0.0;
	double min_val=9999;
	for (QNUInt32 i=0; i<lambda_len; i++) {
		grad[i]=this->ExpF[i];
		if (max_val<grad[i]) { max_val=grad[i];}
		if (min_val>grad[i]) { min_val=grad[i];}
	}
	if (logflag ==0) { cout << "MAX: " << max_val << " MIN: " << min_val << endl; }
/*	try {
		logLi-=logE(Zx);
	}
	catch (exception &e) {
		string errstr="CRF_NewGradBuilder::buildGradient caught exception: "+string(e.what())+" while taking log of scale factor";
		throw runtime_error(errstr);
		return(-1);
	}*/
	*Zx_out=0.0;

	delete phn_lat;
	delete lab_lat;
	delete shortest_fst;
	delete [] truearr;
	delete [] bestarr;
	//delete distance;

	//nodeList.clear();
	return logli;
}
