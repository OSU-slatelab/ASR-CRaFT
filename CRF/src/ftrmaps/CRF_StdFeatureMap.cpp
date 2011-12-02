/*
 * CRF_StdFeatureMap.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_StdFeatureMap.h"

/*
 * CRF_StdFeatureMap constructor
 *
 * Input: nlabs - number of possible labels in the CRF
 *        nfeas - number of features per label in the CRF
 *
 */
CRF_StdFeatureMap::CRF_StdFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas)
	: CRF_FeatureMap(nlabs,nfeas)
{
	this->stateFeatureIdxCache = new QNUInt32[nlabs];
	this->transFeatureIdxCache = new QNUInt32[nlabs*nlabs];

	this->recalc();

}

/*
 * CRF_StdFeatureMap constructor
 *
 * Input: *cnf - struct containing feature mapping information
 *
 * Note:  This constructor is preferred for new code.
 *
 */
CRF_StdFeatureMap::CRF_StdFeatureMap(CRF_FeatureMap_config* cnf)
  :CRF_FeatureMap(cnf)
{
	this->stateFeatureIdxCache = new QNUInt32[config->numLabs];
	this->transFeatureIdxCache = new QNUInt32[config->numLabs*config->numLabs];

	this->recalc();
}

/*
 * CRF_StdFeatureMap destructor
 *
 *
 */
CRF_StdFeatureMap::~CRF_StdFeatureMap()
{
	delete[] this->stateFeatureIdxCache;
	delete[] this->transFeatureIdxCache;
}

/*
 * CRF_StdFeatureMap::computeStateArrayValue
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        clab - label to compute the CRF value for
 *
 * Returns: computed state value for the label clab given the set of features (ftr_buf)
 *          and the current CRF labmda values (lambda)
 */
double CRF_StdFeatureMap::computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab)
{
	double stateValue=0.0;
	QNUInt32 lc = this->stateFeatureIdxCache[clab];
	if (config->useStateFtrs) {
		for (QNUInt32 fidx=config->stateFidxStart; fidx<=config->stateFidxEnd; fidx++)
		{
			stateValue+=ftr_buf[fidx]*lambda[lc];

			// just for debugging
			//cout << "stateValue+=ftr_buf[" << fidx << "](" << ftr_buf[fidx] << ")*lambda[" << lc << "](" << lambda[lc] << ")=" << stateValue << endl;

			lc++;
		}
	}
	if (config->useStateBias) {
		stateValue+=lambda[lc]*config->stateBiasVal;

		// just for debugging
		//cout << "stateValue+=config->stateBiasVal(" << config->stateBiasVal << ")*lambda[" << lc << "](" << lambda[lc] << ")=" << stateValue << endl;

		lc++;
	}
	return stateValue;
}

/*
 * CRF_StdFeatureMap::computeTransMatrixValue
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        plab - the previous label to compute the CRF value from
 *        clab - the current label to compute the CRF value for
 *
 * Returns: computed state value for the transition between labels plab and clab
 *           given the set of features (ftr_buf) and the current CRF labmda values (lambda)
  */
double CRF_StdFeatureMap::computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab)
{
	double transMatrixValue=0.0;
	QNUInt32 lc = this->transFeatureIdxCache[plab*config->numLabs+clab];
	if (config->useTransFtrs) {
		for (QNUInt32 fidx=config->transFidxStart; fidx<=config->transFidxEnd; fidx++)
		{
			transMatrixValue+=ftr_buf[fidx]*lambda[lc];
			lc++;
		}
	}
	if (config->useTransBias) {
		transMatrixValue+=lambda[lc]*config->transBiasVal;
		lc++;
	}
	return transMatrixValue;
}

/*
 * CRF_StdFeatureMap::computeStateExpF
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        *ExpF - vector to store expected values for each state feature
 *        *grad - vector used to store gradient values
 *        alpha_beta - gamma value of the current label (clab)
 *        t_clab - true label (used in training)
 *        clab - label to compute Expected state feature values for
 *        compute_grad - flag to control whether grad vector is updated or not
 *
 * Returns: log likelihood of the state features
 *
 * Fills the ExpF vector and gradient vector with their values using the pre-computed
 * alpha_beta (probability of being in the state with the label clab).
 *
 */
double CRF_StdFeatureMap::computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab, bool compute_grad)
{
	double logLi=0.0;
	QNUInt32 lc = this->stateFeatureIdxCache[clab];
	if (config->useStateFtrs) {
#ifdef VECTOR_ROUTINES
	// Worked this out but then remembered that
    // ExpF, grad, and lambda are doubles so this code won't work
    // We're keeping this here to remind us to implement the following functions
    //
    // craft_mulacc_vfdvd_vd
	// craft_add_vfvd_vd
    // craft_mulsum_vfvd_d

    int veclen=this->stateFidxEnd - this->stateFidxStart;
	qn_mulacc_vffvf_vf(veclen,ftr_buf,alpha_beta,&ExpF[lc],&ExpF[lc]);
	if (t_clab == clab) {
		qn_add_vfvf_vf(veclen,ftr_buf,&grad[lc],&grad[lc]);
		logLi += qn_mulsum_vfvf_f(veclen,ftr_buf,&lambda[lc]);
	}
	lc+=veclen;

#else

		for (QNUInt32 fidx=config->stateFidxStart; fidx<=config->stateFidxEnd; fidx++)
		{
			ExpF[lc]+=alpha_beta*ftr_buf[fidx];

			//just for debugging
//			if (lc == 0)
//			{
//				cout << "StateFeature for lab=" << clab << ": ExpF[" << lc << "]+=" << alpha_beta << "*" << ftr_buf[fidx] << "=" << ExpF[lc] << endl;
//			}

			if (compute_grad && (t_clab == clab)) {
				grad[lc]+=ftr_buf[fidx];

				//just for debugging
//				if (lc == 0)
//				{
//					cout << "StateFeature for lab=" << clab << ": ExpF[" << lc << "]+=" << alpha_beta << "*" << ftr_buf[fidx] << "=" << ExpF[lc] << endl;
//				}

				logLi += lambda[lc]*ftr_buf[fidx];
			}
			lc++;
		}
#endif
	}

	if (config->useStateBias) {
		ExpF[lc]+=alpha_beta*config->stateBiasVal;

		//just for debugging
//		cout << "StateBias for lab=" << clab << ": ExpF[" << lc << "]+=" << alpha_beta << "*" << config->stateBiasVal << "=" << ExpF[lc] << endl;

		if (compute_grad && (t_clab == clab)) {
			grad[lc]+=config->stateBiasVal;

			//just for debugging
//			cout << "StateBias for lab=" << clab << ": grad[" << lc << "]+=" << config->stateBiasVal << "=" << grad[lc] << endl;

			logLi += lambda[lc]*config->stateBiasVal;
		}
		lc++;
	}
	return logLi;
}

/*
 * CRF_StdFeatureMap::computeTransExpF
 *
 * Input: *ftr_buf - vector of observed feature values for computation
 *        *lambda - lambda vector from the CRF
 *        *ExpF - vector to store expected values for each state feature
 *        *grad - vector used to store gradient values
 *        alpha_beta - gamma value of the current transition (plab->clab)
 *        t_plab - true previous label (used in training)
 *        t_clab - true label (used in training)
 *        plab - previous label used for computing expected transition values
 *        clab - label to compute expected transition values
 *        compute_grad - flag to control whether grad vector is updated or not
 *
 * Returns: log likelihood of the state features
 *
 * Fills the ExpF vector and gradient vector with their values using the pre-computed
 * alpha_beta (probability of being in the transition between plab and clab).
 *
 */
double CRF_StdFeatureMap::computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad)
{
	double logLi=0.0;
	QNUInt32 lc = this->transFeatureIdxCache[plab*config->numLabs+clab];

	if (config->useTransFtrs) {
		for (QNUInt32 fidx=config->transFidxStart; fidx<=config->transFidxEnd; fidx++)
		{
			ExpF[lc]+=alpha_beta*ftr_buf[fidx];
			if (compute_grad && (clab==t_clab) && (plab==t_plab)) {
				grad[lc]+=ftr_buf[fidx];
				logLi += lambda[lc]*ftr_buf[fidx];
			}
			lc++;
		}
	}
	if (config->useTransBias) {
		ExpF[lc] += alpha_beta*config->transBiasVal;

		//just for debugging
//		cout << "TransBias for (plab=" << plab << ",clab=" << clab << "): ExpF[" << lc << "]+=" << alpha_beta << "*" << config->transBiasVal << "=" << ExpF[lc] << endl;

		if (compute_grad && (clab==t_clab) && (plab==t_plab)) {
			grad[lc]+=config->transBiasVal;

			//just for debugging
//			cout << "TransBias for (plab=" << plab << ",clab=" << clab << "): grad[" << lc << "]+=" << config->transBiasVal << "=" << grad[lc] << endl;

			logLi+=lambda[lc]*config->transBiasVal;
		}
		lc++;
	}

	return logLi;
}


QNUInt32 CRF_StdFeatureMap::getNumStates()
{
	return config->numStates;
}

/*
 * CRF_StdFeatureMap::getNumStateFuncs
 *
 * Input: clab - label under examination
 *
 * Returns: total number of state functions defined for the label clab
 *
 * In this case, all labels have the same number of state functions
 */
QNUInt32 CRF_StdFeatureMap::getNumStateFuncs(QNUInt32 clab)
{
	return this->numStateFuncs;
}

/*
 * CRF_StdFeatureMap::getNumTransFuncs
 *
 * Input: plab - previous label under examination
 *        clab - label under examination
 *
 * Returns: total number of transition functions defined for the transition plab->clab
 *
 * For CRF_StdFeatureMap this is the same for all transitions.
 */
QNUInt32 CRF_StdFeatureMap::getNumTransFuncs(QNUInt32 plab, QNUInt32 clab)
{
	return this->numTransFuncs;
}

/*
 * CRF_StdFeatureMap::computeStateFeatureIdx
 *
 * Input: clab - label under examination
 *        fno - feature number under examination
 *
 * Returns: the index into the lambda vector for the input feature identified
 *           by fno and the label identified by clab.
 *
 * This is a private function used to create the idxCache.  To access an index, use
 * getStateFeatureIdx instead.
 *
 * EFL: I rewrote this because I don't like the loop in the second part
 * JJM: Internally this function only initializes the lookup table, so it should only get called once
 *      what really needs to be done is that this needs to be moved to a private function with a
 *      new name and replaced with public functions that look at the cached indices and do the correct
 *      thing with them, rather than recomputing them on the fly
 * JJM: This function is now "computeStateFeatureIdx", and is a private function.
 */
#ifndef OLDCODE
QNUInt32 CRF_StdFeatureMap::computeStateFeatureIdx(QNUInt32 clab, QNUInt32 fno)
{
	QNUInt32 retVal=0;
	if (clab == 0) { retVal = fno; }
	else {
		if (config->numStates==1) {
			retVal=clab*(this->numStateFuncs + config->numLabs*this->numTransFuncs);
		}
		else {
			for (QNUInt32 i=0; i<clab; i++) {
				// For each label, add in the state features
				retVal+=this->numStateFuncs;
				// Next, if this is a start state, add in trans features from all previous end
				// states AND A SELF TRANSITION
				if (i % config->numStates == 0) {
					retVal+=this->numActualLabels*this->numTransFuncs + 1;
				}
				else {
					// Otherwise, add in a self transition and a transition from the previous label
					retVal+=2*this->numTransFuncs;
				}
			}
		}
	}
	retVal+=fno;  //Ryan, I think this has problems when fno > 0 since it adds fno twice.
	return retVal;
}
#else
QNUInt32 CRF_StdFeatureMap::getStateFeatureIdx(QNUInt32 clab, QNUInt32 fno) {
	if (this->numStates==1) {
		return clab*(this->numStateFuncs + this->numLabs*this->numTransFuncs)+fno;
	} else {
		// n-state model: need to compute number of initial states prior to this label
		QNUInt32 n_initials=(clab+this->numStates-1)/this->numStates;
		QNUInt32 n_noninitials=clab-n_initials;
		// every state carries numStateFuncs
		// initial state transitions: add in trans features from all previous end
		//                            states AND A SELF TRANSITION
		// non-initial state transitions: add in a self transition and a transition from the previous label
		return (clab-1)*this->numStateFuncs +
			n_initials*(this->numActualLabels*this->numTransFuncs+1) +
			n_noninitials*(2*this->numTransFuncs) +
			fno;
	}
}
#endif

/*
 * CRF_StdFeatureMap::computeTransFeatureIdx
 *
 * Input: clab - label under examination
 *        plab - previous label under examination
 *        fno - feature number under examination
 *
 * Returns: the index into the lambda vector for the input feature identified
 *           by fno and the transition identified by plab->clab.
 *
 * This is a private function used to create the idxCache.  To access an index, use
 * getTransFeatureIdx instead.
 *
 */
QNUInt32 CRF_StdFeatureMap::computeTransFeatureIdx(QNUInt32 clab, QNUInt32 plab, QNUInt32 fno)
{
	QNUInt32 retVal=0;
	if (config->numStates == 1) {
		retVal=clab*(this->numStateFuncs + config->numLabs*this->numTransFuncs) + this->numStateFuncs + plab*this->numTransFuncs + fno;
	}
	else {
		retVal=this->getStateFeatureIdx(clab);
		retVal+=this->numStateFuncs;
		// We're now in the right spot if plab == clab.  Otherwise we have to increment
		if (plab != clab) {
			retVal+=1; // Increment to the next one
			//We're now in the right spot if we have a non-start state clab
			if (clab % config->numStates == 0) {
				QNUInt32 real_plab = (plab+1)/config->numStates-1;
				retVal+=real_plab*this->numTransFuncs;
			}
		}
	}
	retVal+=fno;   //Ryan, I think this has problems when fno > 0 since it adds fno twice.
	return retVal;
}

/*
 * CRF_StdFeatureMap::getStateFeatureIdx
 *
 * Input: clab - label under examination
 *        fno - feature number under examination
 *
 * Returns: index into the lambda vector for the input feature identified by fno
 *            and the label identified by clab.
 */
QNUInt32 CRF_StdFeatureMap::getStateFeatureIdx(QNUInt32 clab, QNUInt32 fno) {
	return this->stateFeatureIdxCache[clab]+fno;
}

/*
 * CRF_StdFeatureMap::getTransFeatureIdx
 *
 * Input: clab - label under examination
 *        plab - previous label under examination
 *        fno - feature number under examination
 *
 * Returns: index into the lambda vector for the input feature identified by fno
 *            and the transition identified by plab->clab
 */
QNUInt32 CRF_StdFeatureMap::getTransFeatureIdx(QNUInt32 clab, QNUInt32 plab, QNUInt32 fno) {
	return this->transFeatureIdxCache[plab*config->numLabs+clab]+fno;
}

/*
 * CRF_StdFeatureMap::getStateBiasIdx
 *
 * Input: clab - label under examination
 *
 * Returns: index into the lambda vector for the state bias feature for the label clab
 */
QNUInt32 CRF_StdFeatureMap::getStateBiasIdx(QNUInt32 clab) {
	//return getStateFeatureIdx(clab,this->numFtrFuncs-1);

	//Changed by Ryan, I think it should be this:
	return getStateFeatureIdx(clab,this->numStateFuncs-1);
}

/*
 * CRF_StdFeatureMap::getTransBiasIdx
 *
 * Input: clab - label under examination
 *        plab - previous label under examination
 *
 * Returns: index into the lambda vector for the state bias feature for the transition
 *          plab->clab
 */
QNUInt32 CRF_StdFeatureMap::getTransBiasIdx(QNUInt32 clab, QNUInt32 plab) {
	return getTransFeatureIdx(clab,plab,this->numTransFuncs-1);
}

/*
 * CRF_StdFeatureMap::recalc
 *
 * Recomputes internal values.  Use after using a mutator function.
 * Returns: total number of state feature functions in the feature map.
 */
QNUInt32 CRF_StdFeatureMap::recalc()
{
	this->numActualLabels = config->numLabs/config->numStates;
	cout << "ACTUAL LABELS COMPUTED: " << this->numActualLabels << endl;
	if (this->numActualLabels * config->numStates != config->numLabs) {
		string errstr="CRF_StdFeatureMap created exception: Invalid state/label combination while computing transitions";
		throw runtime_error(errstr);
	}
	if (config->numStates==1) {
		this->transMult=this->numActualLabels*this->numActualLabels;
	}
	else {
		this->transMult = this->numActualLabels*this->numActualLabels+config->numLabs + config->numLabs-this->numActualLabels;
	}
	// end->start transitions + diagonal self transitions + offDiagonal transitions
	this->numFtrFuncs=0;
	this->numStateFuncs=0;
	this->numTransFuncs=0;
	if (config->useStateFtrs) {
		this->numFtrFuncs+=(config->stateFidxEnd - config->stateFidxStart + 1 )*config->numLabs; // Count start at zero
		this->numStateFuncs+=(config->stateFidxEnd - config->stateFidxStart +1 );
	}
	if (config->useStateBias) {
		this->numFtrFuncs += config->numLabs;
		this->numStateFuncs+=1;
	}
	if (config->useTransFtrs) {
		//this->numFtrFuncs+=(this->transFidxEnd - this->transFidxStart + 1 )*this->numLabs*this->numLabs; //Count still starts at zero
		this->numFtrFuncs+=(config->transFidxEnd - config->transFidxStart + 1 )*this->transMult; //Count still starts at zero
		this->numTransFuncs+=(config->transFidxEnd - config->transFidxStart +1);
	}
	if (config->useTransBias) {
		//this->numFtrFuncs += this->numLabs*this->numLabs;
		this->numFtrFuncs += this->transMult;
		this->numTransFuncs+=1;
	}

	for (QNUInt32 clab=0; clab<config->numLabs; clab++) {
		this->stateFeatureIdxCache[clab]=this->computeStateFeatureIdx(clab);
		for (QNUInt32 plab=0; plab<config->numLabs; plab++) {
			this->transFeatureIdxCache[plab*config->numLabs+clab]=this->computeTransFeatureIdx(clab,plab);
		}
	}

	return this->numFtrFuncs;
}

/*
 * CRF_FeatureMap::getMapDescriptor
 *
 * Returns: a string that describes the feature map for the input lambda index lambdaNum.
 */
string CRF_StdFeatureMap::getMapDescriptor(QNUInt32 lambdaNum) {
	std::stringstream ss;
	QNUInt32 stateOffset=((config->useStateFtrs)?config->numFeas:0)+
						 ((config->useStateBias)?1:0)+
					     ((config->useTransFtrs)?config->numFeas*config->numLabs:0)+
					     ((config->useTransBias)?config->numLabs:0);
	QNUInt32 stateNum=lambdaNum/stateOffset;
	QNUInt32 remainder=lambdaNum-(stateNum*stateOffset);
	if (config->useStateFtrs) {
		if (remainder < config->numFeas) {
			ss<< "State" << stateNum << ":ftr" << remainder;
			return ss.str();
		} else {
			remainder-=config->numFeas;
		}
	}
	if (config->useStateBias) {
		if (remainder == 0) {
			ss << "StateBias" << stateNum;
			return ss.str();
		} else {
			remainder--;
		}
	}
	// not sure about this!
	if (config->useTransFtrs) {
		if (remainder < config->numFeas*config->numLabs) {
			QNInt32 fromState=remainder/config->numFeas;
			QNInt32 ftr=remainder-(fromState*config->numFeas);
			ss << "Trans" << fromState << ">" << stateNum << ":ftr" << ftr;
			return ss.str();
		} else {
			remainder-=(config->numFeas*config->numLabs);
		}
	}
	if (config->useTransBias) {
		ss << "TransBias" << remainder << ">" << stateNum;
		return ss.str();
	}

	return string("Out of range");
}

/*
 * CRF_StdFeatureMap::accumulateFeatures
 *
 * Fills the vector *accumulator with features from  *ftr_buf for the label clab.
 *   Used in training the CRF via AIS.
 */
void CRF_StdFeatureMap::accumulateFeatures(float *ftr_buf,double *accumulator,QNUInt32 clab) {
	int offset=clab*this->numFtrFuncs;
	for (QNUInt32 ii=0;ii<this->numFtrFuncs;ii++) {
		accumulator[offset+ii]+=ftr_buf[ii];
	}
}


/*
 * Added by Ryan
 *
 * CRF_StdFeatureMap::tieGradient
 *
 * Input: *grad - vector of gradient values
 *        maxDur - the maximum duration for the all phone-duration labels
 *        durFtrStart - the start index of the binary-coded features which
 *                      are not for tying
 *
 * Tie gradients for tied parameters. Currently this function only applies to single
 * state model. And the tied parameters are hardcoded as the ones for state features
 * and transition features for the phone-duration labels with the same phone, e.g.,
 * parameters are tied for labels ah-1, ah-2, ah-3, ..., ah-maxDur, except for the
 * duration features, which are binary-coded features with index starting from
 * durFtrStart.
 *
 *
 */
//void CRF_StdFeatureMap::tieGradient(double* grad, QNUInt32 maxDur, QNUInt32 durFtrStart)
void CRF_StdFeatureMap::tieGradient(double* grad)
{

	//cout << "Tying Parameter. maxDur=" << maxDur << " durFtrStart=" << durFtrStart << endl;

	if (config->numStates != 1)
	{
		string errstr="CRF_StdFeatureMap::tieGradient() threw exception: Currently the function of parameter tying only applies to single state models.";
		throw runtime_error(errstr);
	}

	if (config->maxDur <= 0)
		return;

	//assert(config->numLabs % maxDur == 0);
	//assert(config->numLabs == config->maxDur * 48);
	if (config->numLabs % config->maxDur != 0) {
		string errstr="CRF_StdFeatureMap::tieGradient() threw exception: numLabs and maxDur do not correspond.";
		throw runtime_error(errstr);
	}
	QNUInt32 nPureLabs = config->numLabs / config->maxDur;
	assert(nPureLabs == 48);    //this statement is just for debugging: currently the number of labels is 48. It should be commented out later.

	QNUInt32 durFtrEnd = config->durFtrStart + config->maxDur;  // exclusive end
	if (config->useStateFtrs || config->useStateBias)
	{
		//QNUInt32 stateParamStep = getStateFeatureIdx(1) - getStateFeatureIdx(0);
		QNUInt32 stateParamStep = (getStateFeatureIdx(1) - getStateFeatureIdx(0)) * nPureLabs;
		QNUInt32 clab = 0;
		//while (clab < config->numLabs)
		while (clab < nPureLabs)
		{
			for (QNUInt32 fid = 0; fid < numStateFuncs; fid++)
			{
				if (fid >= config->durFtrStart && fid < durFtrEnd)
					continue;
				QNUInt32 start = getStateFeatureIdx(clab,fid);
				tieGradForSingleParam(grad, start, stateParamStep, config->maxDur);
			}
			//clab += config->maxDur;
			clab++;
		}
	}

	if (config->useTransFtrs || config->useTransBias)
	{
		//QNUInt32 clabParamStep = getTransFeatureIdx(1,0) - getTransFeatureIdx(0,0);
		QNUInt32 clabParamStep = (getTransFeatureIdx(1,0) - getTransFeatureIdx(0,0)) * nPureLabs;
		QNUInt32 clab = 0;
		//while (clab < config->numLabs)
		while (clab < nPureLabs)
		{
			for (QNUInt32 plab = 0; plab < config->numLabs; plab++)
			{
				for (QNUInt32 fid = 0; fid < numTransFuncs; fid++)
				{
					if (fid >= config->durFtrStart && fid < durFtrEnd)
						continue;
					QNUInt32 start = getTransFeatureIdx(clab,plab,fid);
					tieGradForSingleParam(grad, start, clabParamStep, config->maxDur);
				}
			}
			//clab += config->maxDur;
			clab++;
		}
		//QNUInt32 plabParamStep = getTransFeatureIdx(0,1) - getTransFeatureIdx(0,0);
		QNUInt32 plabParamStep = (getTransFeatureIdx(0,1) - getTransFeatureIdx(0,0)) * nPureLabs;
		QNUInt32 plab = 0;
		//while (plab < config->numLabs)
		while (plab < nPureLabs)
		{
			for (QNUInt32 clab = 0; clab < config->numLabs; clab++)
			{
				for (QNUInt32 fid = 0; fid < numTransFuncs; fid++)
				{
					if (fid >= config->durFtrStart && fid < durFtrEnd)
						continue;
					QNUInt32 start = getTransFeatureIdx(clab,plab,fid);
					tieGradForSingleParam(grad, start, plabParamStep, config->maxDur);
				}
			}
			//plab += config->maxDur;
			plab++;
		}
	}

}

/*
 * Added by Ryan
 *
 * CRF_StdFeatureMap::tieGradForSingleParam
 *
 * Input: *grad - vector of gradient values
 *        numTiedParam - the number of parameters to be tied for a single value
 *        start - the start index of the series of parameters
 *        step - the difference of the indice of each two consecutive tied parameters
 *               in the series
 *
 * The subroutine for tying the gradients of a series of parameters to a single parameter.
 * The index of parameters in the series increases with the same step value.
 *
 */
void CRF_StdFeatureMap::tieGradForSingleParam(double* grad, QNUInt32 start, QNUInt32 step, QNUInt32 numTiedParam)
{
//	if (step <= 0)
//		return;
	assert(step > 0);
	double gradSum = 0.0;
	QNUInt32 gradID = start;
	for (QNUInt32 i = 0; i < numTiedParam; i++)
	{
		gradSum += grad[gradID];
		gradID += step;
	}
	gradID = start;
	for (QNUInt32 i = 0; i < numTiedParam; i++)
	{
		grad[gradID] = gradSum;
		gradID += step;
	}
}
