/*
 * CRF_NewGradBuilderLogSoft.cpp
 *
 *  Created on: Oct 6, 2009
 *      Author: morrijer
 */

#include "CRF_NewGradBuilderSoft.h"
#include "fst/lib/fstlib.h"

CRF_NewGradBuilderSoft::CRF_NewGradBuilderSoft(CRF_Model* crf_in)
	: CRF_GradBuilder(crf_in)
{
	for (QNUInt32 i=0; i<this->num_labs; i++) {
		this->alpha_base[i]=0.0;
	}
	this->alphaAlignBase.assign(this->num_labs,0);
	// nodelist may come from base constructor
	if (this->nodeList) {
		delete this->nodeList;
	}
	this->lb=NULL;
	this->nodeList=NULL;
}

CRF_NewGradBuilderSoft::~CRF_NewGradBuilderSoft() {
	// TODO Auto-generated destructor stub
}

double CRF_NewGradBuilderSoft::buildGradient(CRF_FeatureStream* ftr_strm, double* grad, double* Zx_out)
{
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

	VectorFst<LogArc> phn_lat;
	VectorFst<LogArc> lab_lat;

	//cerr << "Building lattice" << endl;
	int nstates=0;
	if (this->crf->getFeatureMap()->getNumStates()==1) {
		nstates=lb->buildLattice(&phn_lat,true,&lab_lat,false);
	}
	else {
		nstates=lb->nStateBuildLattice(&phn_lat,true,&lab_lat,false);
	}

	ComposeFst<LogArc> result(phn_lat,lab_lat);

	this->lb->computeAlignedAlphaBeta(result,nstates);

	double logLi = 0.0;

	for (QNUInt32 i=0; i<lambda_len; i++) {
		this->ExpF[i]=0.0;
	}

	QNUInt32 nodeCnt=this->nodeList->getNodeCount()-1;
	// Alphas are already computed, so we just need the beta values
	QNUInt32 lastNode=nodeCnt;
	double Zx=this->nodeList->at(lastNode)->computeAlphaSum();
	double Zx_aligned=this->nodeList->at(lastNode)->computeAlphaAlignedSum();

	bool stop=false;
	while (!stop) {
		double* beta = this->nodeList->at(nodeCnt)->getBeta();
		vector<double>* betaAligned=this->nodeList->at(nodeCnt)->getBetaAligned();
		logLi-=this->nodeList->at(nodeCnt)->getAlphaScale();
		if (nodeCnt==lastNode) {
			this->nodeList->at(nodeCnt)->setTailBeta();
		}
		else {
			// We compute the beta value for the node following our current one, and store the result
			// as the beta for our current node (as per the equations).
			this->nodeList->at(nodeCnt+1)->computeBeta(beta,this->nodeList->at(nodeCnt)->getAlphaScale());
		}
		double* prev_alpha;
		vector<double>* prev_alphaAligned;
		QNUInt32 prev_lab;
		if (nodeCnt>0) {
			prev_alpha=this->nodeList->at(nodeCnt-1)->getAlpha();
			prev_lab = this->nodeList->at(nodeCnt-1)->getLabel();
			prev_alphaAligned = this->nodeList->at(nodeCnt-1)->getAlphaAligned();
		}
		else {
			prev_alpha=this->alpha_base;
			prev_lab=this->num_labs+1;
			prev_alphaAligned = &(this->alphaAlignBase);
		}
		//double cur_alpha_sum = this->nodeList->at(nodeCnt)->computeAlphaSum(); //*DEBUG*//
		logLi += this->nodeList->at(nodeCnt)->computeExpF(this->ExpF, grad, Zx, prev_alpha, prev_lab);
		// Here we replace the ExpF calculation with our new soft ExpF calculation

		//cout << "\t" << nodeCnt << ":\tLogLi is now: " << logLi << "\tAlpha Sum: " << cur_alpha_sum << endl;
		if (nodeCnt==0) { stop=true;} // nodeCnt is unsigned, so we can't do the obvious loop control here
		nodeCnt--;
	}

	for (QNUInt32 i=0; i<lambda_len; i++) {
		grad[i]-=this->ExpF[i];
	}
	*Zx_out=Zx;
	//logLi-=Zx;


	//nodeList.clear();
	return logLi;
}
