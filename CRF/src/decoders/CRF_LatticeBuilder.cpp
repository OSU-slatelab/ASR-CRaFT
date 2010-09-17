/*
 * CRF_LatticeBuilder.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_LatticeBuilder.h"

/*
 * CRF_LatticeBuilder constructor
 *
 * Input: *ftr_stream_in - pointer to input stream of features
 *        *crf_in - pointer to the CRF model to be used for building lattices
 *
 */

CRF_LatticeBuilder::CRF_LatticeBuilder(CRF_FeatureStream* ftr_strm_in, CRF_Model* crf_in)
	: crf(crf_in),
	  ftr_strm(ftr_strm_in)
{
	this->nodeList= new CRF_StateVector();
	this->bunch_size=1;
	this->num_ftrs=this->ftr_strm->num_ftrs();
	this->num_labs=this->crf->getNLabs();
	this->ftr_buf = new float[num_ftrs*bunch_size];
	this->lab_buf = new QNUInt32[bunch_size];
	this->alpha_base = new double[this->num_labs];
	for (QNUInt32 i=0; i<this->num_labs; i++) {
		this->alpha_base[i]=0.0;
	}
}

/*
 * CRF_LatticeBuilder destructor
 *
 */
CRF_LatticeBuilder::~CRF_LatticeBuilder()
{
	delete [] this->ftr_buf;
	delete [] this->lab_buf;
	delete [] this->alpha_base;
	delete this->nodeList;
}

/*
 * CRF_LatticeBuilder::buildLattice
 *
 * Input: none
 *
 * Returns: OpenFst lattice using the current sequence of observations in the
 *          input stream
 *
 * Wrapper function for the templated buildLattice function.  Builds and returns
 * the most common kind of lattice for CRF processing.
 */

StdVectorFst* CRF_LatticeBuilder::buildLattice()
{
	VectorFst<StdArc> *fst=new StdVectorFst();
	VectorFst<StdArc> *alignFst=NULL;
	this->buildLattice(fst,false,alignFst);
	return fst;
}



/*
 * CRF_LatticeBuilder::getAlignmentGammas
 *
 * Input: *denominatorStateGamma
 *        *numeratorStateGamma
 *        *denominatorTransGamma
 *        *numeratorTransGamma
 *
 *        Each of these are pointers to vectors of doubles used to store the
 *        computation of the gamma values for each state of the CRF computation.
 *        The length of these vectors will be:
 *          (number of labels)*(length of sequence)
 *        for the state gammas and
 *          (number of labels)*(number of labels)*(length of sequence)
 *        for the transition gammas
 *
 * Returns: number of states in the denominator lattice.
 *
 * Reads the current sequence of feature observations from the input feature
 * stream and uses them to construct an OpenFst lattice.  Uses this lattice to
 * compute the gamma values of the foward-backward computation for each timestep
 * (i.e. the normalized alpha*beta values for each label for each sequence element).
 *
 * Denominator values are computed using the whole lattice.  Numerator values
 * are computed using the denominator lattice aligned with a lattice constructed
 * from the labels in the feature stream.
 *
 */

int CRF_LatticeBuilder::getAlignmentGammas(vector<double> *denominatorStateGamma,
											vector<double> *numeratorStateGamma,
											vector<double> *denominatorTransGamma,
											vector<double> *numeratorTransGamma) {

	// normally, you'd want to use StdArc, but we use LogArc so we can get
	// alphas/betas.  Critically, the aligner is unweighted, so that being
	// in the LogArc semiring doesn't matter.
	VectorFst<LogArc> denominator;
	VectorFst<LogArc> aligner;

	//cout << "In getAlignmentGammas" << endl;
	int nstates=this->buildLattice(&denominator,(numeratorStateGamma!=NULL),&aligner);
	//cout << "...(lattice)" << endl;

	// do the denominator first
	if (denominatorStateGamma != NULL) {
		//cout << "Computing gamma for denominator" << endl;
		_computeGamma(denominatorStateGamma,denominatorTransGamma,denominator,nstates);
		//cout << "Done computing gamma" << endl;
        //cout << "...(denominator gamma)" << endl;
	}

	if (numeratorStateGamma) {
		ComposeFst<LogArc> numerator(denominator,aligner);
		//cout << "start: " << numerator.Start() << " numarcs: " << numerator.NumArcs(numerator.Start());
		//Project(numerator,PROJECT_OUTPUT);
		//RmEpsilon(numerator);
		//TopSort(numerator);

		_computeGamma(numeratorStateGamma,numeratorTransGamma,numerator,nstates);
		//cout << "...(numerator gamma)" << endl;
	}

	return nstates;
}

/*
 * CRF_LatticeBuilder::getAlignmentGammasNState
 *
 * Input: *denominatorStateGamma
 *        *numeratorStateGamma
 *        *denominatorTransGamma
 *        *numeratorTransGamma
 *
 *        Each of these are pointers to vectors of doubles used to store the
 *        computation of the gamma values for each state of the CRF computation.
 *        The length of these vectors will be:
 *          (number of labels)*(length of sequence)
 *        for the state gammas and
 *          (number of labels)*(number of labels)*(length of sequence)
 *        for the transition gammas
 *
 * Returns: number of states in the denominator lattice.
 *
 * Performs the same function as getAlignmentGammas above, but with builds the
 * initial lattice with a multi-state restriction on the paths that the CRF is
 * allowed to follow.
 *
 */
int CRF_LatticeBuilder::getAlignmentGammasNState(vector<double> *denominatorStateGamma,
											vector<double> *numeratorStateGamma,
											vector<double> *denominatorTransGamma,
											vector<double> *numeratorTransGamma) {

	// normally, you'd want to use StdArc, but we use LogArc so we can get
	// alphas/betas.  Critically, the aligner is unweighted, so that being
	// in the LogArc semiring doesn't matter.
	VectorFst<LogArc> denominator;
	VectorFst<LogArc> aligner;

	//cout << "In getAlignmentGammas" << endl;
	int nstates=this->nStateBuildLattice(&denominator,(numeratorStateGamma!=NULL),&aligner);
	//cout << "...(lattice)" << endl;

	// do the denominator first
	if (denominatorStateGamma != NULL) {
		//cout << "Computing gamma for denominator" << endl;
		_computeGamma(denominatorStateGamma,denominatorTransGamma,denominator,nstates);
		//cout << "Done computing gamma" << endl;
        //cout << "...(denominator gamma)" << endl;
	}

	if (numeratorStateGamma) {
		ComposeFst<LogArc> numerator(denominator,aligner);
		//cout << "start: " << numerator.Start() << " numarcs: " << numerator.NumArcs(numerator.Start());
		//Project(numerator,PROJECT_OUTPUT);
		//RmEpsilon(numerator);
		//TopSort(numerator);

		_computeGamma(numeratorStateGamma,numeratorTransGamma,numerator,nstates);
		//cout << "...(numerator gamma)" << endl;
	}

	return nstates;
}



void CRF_LatticeBuilder::_computeGamma(vector<double> *stateGamma, vector<double> *transGamma, Fst<LogArc> &fst,int nstates) {
	vector<LogWeight> alpha;
	vector<LogWeight> beta;

	//cout << "In computeGamma" << endl;
	// when shortest distance is done with log semiring, we get alphas
	ShortestDistance(fst,&alpha,false);
	//cout << "Got alphas, size " << alpha.size() << endl;
	//for(int i=0;i<alpha.size();i++) {
	//	cout << alpha[i] << " ";
	//}
	// run it in reverse and you get betas
	ShortestDistance(fst,&beta,true);
	//cout << "Got betas, size " << beta.size() << endl;
	//for(int i=0;i<alpha.size();i++) {
	//	cout << beta[i] << " ";
	//}

	//cout << "nstates=" << nstates << endl;

	// set up the vector
	stateGamma->reserve(nstates*this->num_labs);
	stateGamma->assign(nstates*this->num_labs,CRF_LogMath::LOG0);

	if (transGamma != NULL) {
		transGamma->reserve(nstates*this->num_labs*this->num_labs);
		transGamma->assign(nstates*this->num_labs*this->num_labs,
									  CRF_LogMath::LOG0);
	}
	//cout << "Set up vectors" << endl;
	//cout << "number of states: " << nstates << endl;

	vector<int> positions(alpha.size(),0);

	for (StateIterator<Fst<LogArc> > siter(fst);!siter.Done(); siter.Next()) {

		LogArc::StateId s=siter.Value();
		int crfstate=positions[s];
		int stateoffset=crfstate*this->num_labs;

		// ok maybe this isn't the best way.  This is to cure the problem of
		// having a null transition on the last arc, but crucially it assumes
		// that this arc has no cost.
		if (crfstate>=nstates)
			continue;

		for (ArcIterator< Fst<LogArc> > aiter(fst,s);!aiter.Done(); aiter.Next()) {
			const LogArc &arc = aiter.Value();

			positions[arc.nextstate]=crfstate+1;
			double alpha_s=(double)(alpha[s].Value());
			double beta_next=(double)(beta[arc.nextstate].Value());

			double gammaarc=-1*(alpha_s+arc.weight.Value()+beta_next);
			int label=arc.olabel-1;

			// should't see these but just in case
			if (label<0)
				continue;

			stateGamma->at(stateoffset+label)=
				CRF_LogMath::logAdd(stateGamma->at(stateoffset+label),
									gammaarc);

#ifdef DEBUG
			cout << "state:" << s << " pos:" << positions[s] << " alpha:" << alpha_s <<
						" weight: " << arc.weight.Value()
						<< " label: " << label
						<< " gammaindex: " << stateoffset+label
						<< " next:" << arc.nextstate << " beta:" << beta_next
						<< " gammaarc: " << gammaarc << " totalgamma: " << stateGamma->at(stateoffset+label)
						<< endl;
#endif

			// if there are transition gammas needed, then look over pairs of
			// labels
			if (transGamma != NULL) {
				int transoffset=(stateoffset+label)*this->num_labs;
				for (ArcIterator< Fst<LogArc> >  aiter2(fst,arc.nextstate);
					 !aiter2.Done();
					 aiter2.Next()) {
					const LogArc &arc2 = aiter2.Value();
					positions[arc2.nextstate]=crfstate+1;
					double gammatransarc=-1*(alpha_s+
											 arc.weight.Value()+
											 arc2.weight.Value()+
											 beta[arc2.nextstate].Value());
					int label2=arc2.olabel-1;
					transGamma->at(transoffset+label2)=
						CRF_LogMath::logAdd(transGamma->at(transoffset+label2),
											gammatransarc);
				}
			}

		}
	}
	//cout << "computed gammas" << endl;
	// now normalize to get probability
	int offset=0;
	for (int i=0;i<nstates;i++) {
		double total=CRF_LogMath::LOG0;
		for (int j=0;j<this->num_labs;j++) {
			total=CRF_LogMath::logAdd(total,stateGamma->at(offset+j));
			//cout << "pos: " << i << " label: " << j << " unnormgamma: " << stateGamma->at(offset+j) << endl;
		}
		for (int j=0;j<this->num_labs;j++) {
			stateGamma->at(offset+j)-=total;
			//cout << "pos: " << i << " label: " << j << " gamma: " << stateGamma->at(offset+j) << endl;
		}

		if (transGamma != NULL) {
			total=CRF_LogMath::LOG0;
			for(int j=0;j<this->num_labs;j++) {
				int transoffset=(offset+j)*this->num_labs;
				for(int k=0;k<this->num_labs;k++) {
					total=CRF_LogMath::logAdd(total,
										  transGamma->at(transoffset+k));
				}
			}
			for(int j=0;j<this->num_labs;j++) {
				int transoffset=(offset+j)*this->num_labs;
				for(int k=0;k<this->num_labs;k++) {
					transGamma->at(transoffset+k)-=total;
				}
			}
		}
		offset+=this->num_labs;
	}
	//cout << "normalized gammas" << endl;

}

/*
 * CRF_LatticeBuilder::getNodeList
 *
 *  Accessor function for the nodeList vector generated during buildLattice
 *  processing.  Useful if further processing is needed on the nodes, rather
 *  than on the lattice.
 */

CRF_StateVector * CRF_LatticeBuilder::getNodeList() {
	return this->nodeList;
}

/*
 * CRF_LatticeBuilder::computeAlignedAlphaBeta
 *
 * Input: &fst - previously built CRF lattice in OpenFst format
 *        nstates - number of states in the OpenFst formatted fst object
 *
 *  Performs a function similar to getAlignmentGammas above, but instead of
 *  getting the gamma values, stores the alpha and beta values for each step
 *  separately for further processing.  Alpha and beta values are stored in
 *  the appropriate node in nodeList, which can be accessed using the
 *  getNodeList function above.
 */

void CRF_LatticeBuilder::computeAlignedAlphaBeta(Fst<LogArc>& fst, int nstates) {
	vector<LogWeight> alpha;
	vector<LogWeight> beta;

	// when shortest distance is done with log semiring, we get alphas
	ShortestDistance(fst,&alpha,false);
	// run in reverse to get betas
	ShortestDistance(fst,&beta,true);

	// next we store the alphas and the betas in their proper places in the nodeList
	// Because we don't know the structure of the input fst a priori, we need to build this
	// based on the arc labels in the lattice and store the values from the larger alpha and
	// beta arrays in the right places

	vector<int> positions(alpha.size(),0);
	vector<int> labels(alpha.size(),0);
	vector<double>* curAlpha=NULL;
	vector<double>* curBeta=NULL;
	for (StateIterator<Fst<LogArc> > siter(fst);!siter.Done(); siter.Next()) {

		LogArc::StateId s=siter.Value();
		int crfstate=positions[s];

		// ok maybe this isn't the best way.  This is to cure the problem of
		// having a null transition on the last arc, but crucially it assumes
		// that this arc has no cost.
		if (crfstate>=nstates)
			continue;

		// To conform with our standard method of alphas and betas, nodelist 0 contains the
		// data from the first observation, so we need to make exceptions for it.
		// (the first state should have an alpha of 0 and a beta of shortest path anyway...)
		if (crfstate!=0) {
			curAlpha=this->nodeList->at(crfstate-1)->getAlphaAligned();
			curBeta=this->nodeList->at(crfstate-1)->getBetaAligned();
		}
		for (ArcIterator< Fst<LogArc> > aiter(fst,s);!aiter.Done(); aiter.Next()) {
			const LogArc &arc = aiter.Value();
			// The labels array contains what label the state at the end of this arc should have
			labels[arc.nextstate]=arc.olabel-1;
			// Positions array tells us what timestep the state at the end of this arc occurs at
			positions[arc.nextstate]=crfstate+1;

			if (crfstate!=0) {
				// skip the first state because there are no alpha or beta arrays here
				curAlpha->at(labels[s])=(double)(alpha[s].Value());
				curBeta->at(labels[s])=(double)(beta[s].Value());
			}
		}
	}
}

