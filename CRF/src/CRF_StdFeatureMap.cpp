#include "CRF_StdFeatureMap.h"

CRF_StdFeatureMap::CRF_StdFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas)
	: CRF_FeatureMap(nlabs,nfeas),
	  numLabs(nlabs),
	  numFeas(nfeas)
{
	//this->numFtrFuncs=nlabs*nfeas + nlabs + nlabs*nlabs;
	this->useStateBias=true;
	this->useTransBias=true;
	this->useStateFtrs=true;
	this->useTransFtrs=false;
	this->stateFidxStart=0;
	this->stateFidxEnd=nfeas-1;
	this->transFidxEnd=nlabs+1;
	this->transFidxStart=0;
	this->numStates=1;
	this->stateBiasVal=1.0;
	this->transBiasVal=1.0;
	this->stateFeatureIdxCache = new QNUInt32[nlabs];
	this->transFeatureIdxCache = new QNUInt32[nlabs*nlabs];

	this->recalc();

}

CRF_StdFeatureMap::~CRF_StdFeatureMap()
{
	delete[] this->stateFeatureIdxCache;
	delete[] this->transFeatureIdxCache;
}

double CRF_StdFeatureMap::computeRi(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 clab)
{
	double Ri=0.0;
	QNUInt32 tmp_lc = this->stateFeatureIdxCache[clab];
//	QNUInt32 tmp_lc = this->getStateFeatureIdx(clab);
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=this->stateFidxStart; fidx<=this->stateFidxEnd; fidx++)
		{
			Ri+=ftr_buf[fidx]*lambda[tmp_lc];
			tmp_lc++;
		}
	}
	if (this->useStateBias) {
		Ri+=lambda[tmp_lc]*this->stateBiasVal;
		tmp_lc++;
	}
	lc=tmp_lc;
	return Ri;
}

double CRF_StdFeatureMap::computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab)
{
	double stateValue=0.0;
	QNUInt32 lc = this->stateFeatureIdxCache[clab];
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=this->stateFidxStart; fidx<=this->stateFidxEnd; fidx++)
		{
			stateValue+=ftr_buf[fidx]*lambda[lc];
			lc++;
		}
	}
	if (this->useStateBias) {
		stateValue+=lambda[lc]*this->stateBiasVal;
		lc++;
	}
	return stateValue;
}

double CRF_StdFeatureMap::computeMij(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 plab, QNUInt32 clab)
{
	double Mi=0.0;
	QNUInt32 tmp_lc = this->transFeatureIdxCache[plab*this->numLabs+clab];
	//QNUInt32 tmp_lc = this->getTransFeatureIdx(clab,plab);
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=this->transFidxStart; fidx<=this->transFidxEnd; fidx++)
		{
			Mi+=ftr_buf[fidx]*lambda[tmp_lc];
			tmp_lc++;
		}
	}
	if (this->useTransBias) {
		Mi+=lambda[tmp_lc]*this->transBiasVal;
		tmp_lc++;
	}
	lc=tmp_lc;
	return Mi;
}

double CRF_StdFeatureMap::computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab)
{
	double transMatrixValue=0.0;
	QNUInt32 lc = this->transFeatureIdxCache[plab*this->numLabs+clab];
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=this->transFidxStart; fidx<=this->transFidxEnd; fidx++)
		{
			transMatrixValue+=ftr_buf[fidx]*lambda[lc];
			lc++;
		}
	}
	if (this->useTransBias) {
		transMatrixValue+=lambda[lc]*this->transBiasVal;
		lc++;
	}
	return transMatrixValue;
}

double CRF_StdFeatureMap::computeExpFState(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab)
{
	double logLi=0.0;
	QNUInt32 tmp_lc = this->stateFeatureIdxCache[clab];
	//QNUInt32 tmp_lc = this->getStateFeatureIdx(clab);
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=this->stateFidxStart; fidx<=this->stateFidxEnd; fidx++)
		{
			ExpF[tmp_lc]+=alpha_beta*ftr_buf[fidx];
			if (match) {
				grad[tmp_lc]+=ftr_buf[fidx];
				logLi += lambda[tmp_lc]*ftr_buf[fidx];
			}
			tmp_lc++;
		}
	}
	if (this->useStateBias) {
		ExpF[tmp_lc]+=alpha_beta*this->stateBiasVal;
		if (match) {
			grad[tmp_lc]+=this->stateBiasVal;
			logLi += lambda[tmp_lc]*this->stateBiasVal;
		}
		tmp_lc++;
	}
	lc=tmp_lc;
	return logLi;
}

double CRF_StdFeatureMap::computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab)
{
	double logLi=0.0;
	QNUInt32 tmp_lc = this->transFeatureIdxCache[plab*this->numLabs+clab];
	//QNUInt32 tmp_lc = this->getTransFeatureIdx(clab,plab);
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=this->transFidxStart; fidx<=this->transFidxEnd; fidx++)
		{
			ExpF[tmp_lc]+=alpha_beta*ftr_buf[fidx];
			if (match) {
				grad[tmp_lc]+=ftr_buf[fidx];
				logLi += lambda[tmp_lc]*ftr_buf[fidx];
			}
			tmp_lc++;
		}
	}
	if (this->useTransBias) {
		ExpF[tmp_lc] += alpha_beta*this->transBiasVal;
		if (match) {
			grad[tmp_lc]+=this->transBiasVal;
			logLi+=lambda[tmp_lc]*this->transBiasVal;

		}
		tmp_lc++;
	}
	lc=tmp_lc;
	return logLi;
}

double CRF_StdFeatureMap::computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab, bool compute_grad)
{
	double logLi=0.0;
	QNUInt32 lc = this->stateFeatureIdxCache[clab];
	if (this->useStateFtrs) {
#ifdef VECTOR_ROUTINES
	// Eric and Jeremy worked this out but then remembered that
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

		for (QNUInt32 fidx=this->stateFidxStart; fidx<=this->stateFidxEnd; fidx++)
		{
			ExpF[lc]+=alpha_beta*ftr_buf[fidx];
			if (compute_grad && (t_clab == clab)) {
				grad[lc]+=ftr_buf[fidx];
				logLi += lambda[lc]*ftr_buf[fidx];
			}
			lc++;
		}
#endif
	}

	if (this->useStateBias) {
		ExpF[lc]+=alpha_beta*this->stateBiasVal;
		if (compute_grad && (t_clab == clab)) {
			grad[lc]+=this->stateBiasVal;
			logLi += lambda[lc]*this->stateBiasVal;
		}
		lc++;
	}
	return logLi;
}

double CRF_StdFeatureMap::computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab, bool compute_grad)
{
	double logLi=0.0;
	QNUInt32 lc = this->transFeatureIdxCache[plab*this->numLabs+clab];

	if (this->useTransFtrs) {
		for (QNUInt32 fidx=this->transFidxStart; fidx<=this->transFidxEnd; fidx++)
		{
			ExpF[lc]+=alpha_beta*ftr_buf[fidx];
			if (compute_grad && (clab==t_clab) && (plab==t_plab)) {
				grad[lc]+=ftr_buf[fidx];
				logLi += lambda[lc]*ftr_buf[fidx];
			}
			lc++;
		}
	}
	if (this->useTransBias) {
		ExpF[lc] += alpha_beta*this->transBiasVal;
		if (compute_grad && (clab==t_clab) && (plab==t_plab)) {
			grad[lc]+=this->transBiasVal;
			logLi+=lambda[lc]*this->transBiasVal;
		}
		lc++;
	}

	return logLi;
}

void CRF_StdFeatureMap::setStateFtrRange(QNUInt32 st, QNUInt32 end)
{
	this->stateFidxStart=st;
	this->stateFidxEnd=end;
	if (this->stateFidxEnd - this->stateFidxStart +1 >0 ) {
		this->useStateFtrs=true;
	}
	this->recalc();
}

void CRF_StdFeatureMap::setTransFtrRange(QNUInt32 st, QNUInt32 end)
{
	this->transFidxStart=st;
	this->transFidxEnd=end;
	if (this->transFidxEnd - this->transFidxStart +1 > 0) {
		this->useTransFtrs=true;
	}
	this->recalc();
}


void CRF_StdFeatureMap::setNumStates(QNUInt32 ns)
{
	this->numStates=ns;
	this->recalc();
}

void CRF_StdFeatureMap::setUseStateBias(bool useState)
{
	this->useStateBias=useState;
	this->recalc();
}

void CRF_StdFeatureMap::setUseTransBias(bool useTrans)
{
	this->useTransBias=useTrans;
	this->recalc();
}

void CRF_StdFeatureMap::setUseStateFtrs(bool useState)
{
	this->useStateFtrs = useState;
	this->recalc();
}

void CRF_StdFeatureMap::setUseTransFtrs(bool useTrans)
{
	this->useTransFtrs = useTrans;
	if (this->transFidxEnd == this->numLabs+1) { // By default if it isn't set, use the same features as the state features
		this->transFidxStart=this->stateFidxStart;
		this->transFidxEnd=this->stateFidxEnd;
	}
	this->recalc();
}

void CRF_StdFeatureMap::setStateBiasVal(double stateBias)
{
	this->stateBiasVal=stateBias;
}

void CRF_StdFeatureMap::setTransBiasVal(double transBias)
{
	this->transBiasVal=transBias;
}

QNUInt32 CRF_StdFeatureMap::getNumStates()
{
	return this->numStates;
}

QNUInt32 CRF_StdFeatureMap::getNumStateFuncs(QNUInt32 clab)
{
/*	QNUInt32 retVal=0;
	if (this->useStateFtrs)
	{
		retVal+=this->stateFidxEnd - this->stateFidxStart +1;
	}
	if (this->useStateBias)
	{
		retVal+=1;
	}
	return retVal;*/
	return this->numStateFuncs;
}

QNUInt32 CRF_StdFeatureMap::getNumTransFuncs(QNUInt32 plab, QNUInt32 clab)
{
/*	QNUInt32 retVal=0;
	if (this->useTransFtrs)
	{
		retVal+=this->transFidxEnd - this->transFidxStart +1;
	}
	if (this->useTransBias)
	{
		retVal+=1;
	}
	return retVal;*/
	return this->numTransFuncs;
}

// EFL: I rewrote this because I don't like the loop in the second part
#ifdef OLDCODE
QNUInt32 CRF_StdFeatureMap::getStateFeatureIdx(QNUInt32 clab, QNUInt32 fno)
{
	QNUInt32 retVal=0;
	if (clab == 0) { retVal = fno; }
	else {
		if (this->numStates==1) {
			retVal=clab*(this->numStateFuncs + this->numLabs*this->numTransFuncs);
		}
		else {
			for (QNUInt32 i=0; i<clab; i++) {
				// For each label, add in the state features
				retVal+=this->numStateFuncs;
				// Next, if this is a start state, add in trans features from all previous end
				// states AND A SELF TRANSITION
				if (i % this->numStates == 0) {
					retVal+=this->numActualLabels*this->numTransFuncs + 1;
				}
				else {
					// Otherwise, add in a self transition and a transition from the previous label
					retVal+=2*this->numTransFuncs;
				}
			}
		}
	}
	retVal+=fno;
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

QNUInt32 CRF_StdFeatureMap::getTransFeatureIdx(QNUInt32 clab, QNUInt32 plab, QNUInt32 fno)
{
	QNUInt32 retVal=0;
	if (this->numStates == 1) {
		retVal=clab*(this->numStateFuncs + this->numLabs*this->numTransFuncs) + this->numStateFuncs + plab*this->numTransFuncs + fno;
	}
	else {
		retVal=this->getStateFeatureIdx(clab);
		retVal+=this->numStateFuncs;
		// We're now in the right spot if plab == clab.  Otherwise we have to increment
		if (plab != clab) {
			retVal+=1; // Increment to the next one
			//We're now in the right spot if we have a non-start state clab
			if (clab % this->numStates == 0) {
				QNUInt32 real_plab = (plab+1)/this->numStates-1;
				retVal+=real_plab*this->numTransFuncs;
			}
		}
	}
	retVal+=fno;
	return retVal;
}

QNUInt32 CRF_StdFeatureMap::getStateBiasIdx(QNUInt32 clab) {
	return getStateFeatureIdx(clab,this->numFtrFuncs-1);
}

QNUInt32 CRF_StdFeatureMap::getTransBiasIdx(QNUInt32 clab, QNUInt32 plab) {
	return getTransFeatureIdx(clab,plab,this->numTransFuncs-1);
}

QNUInt32 CRF_StdFeatureMap::recalc()
{
	this->numActualLabels = this->numLabs/this->numStates;
	cout << "ACTUAL LABELS COMPUTED: " << this->numActualLabels << endl;
	if (this->numActualLabels * this->numStates != this->numLabs) {
		string errstr="CRF_StdFeatureMap created exception: Invalid state/label combination while computing transitions";
		throw runtime_error(errstr);
	}
	if (this->numStates==1) {
		this->transMult=this->numActualLabels*this->numActualLabels;
	}
	else {
		this->transMult = this->numActualLabels*this->numActualLabels+this->numLabs + this->numLabs-this->numActualLabels;
	}
	// end->start transitions + diagonal self transitions + offDiagonal transitions
	this->numFtrFuncs=0;
	this->numStateFuncs=0;
	this->numTransFuncs=0;
	if (this->useStateFtrs) {
		this->numFtrFuncs+=(this->stateFidxEnd - this->stateFidxStart + 1 )*this->numLabs; // Count start at zero
		this->numStateFuncs+=(this->stateFidxEnd - this->stateFidxStart +1 );
	}
	if (this->useStateBias) {
		this->numFtrFuncs += this->numLabs;
		this->numStateFuncs+=1;
	}
	if (this->useTransFtrs) {
		//this->numFtrFuncs+=(this->transFidxEnd - this->transFidxStart + 1 )*this->numLabs*this->numLabs; //Count still starts at zero
		this->numFtrFuncs+=(this->transFidxEnd - this->transFidxStart + 1 )*this->transMult; //Count still starts at zero
		this->numTransFuncs+=(this->transFidxEnd - this->transFidxStart +1);
	}
	if (this->useTransBias) {
		//this->numFtrFuncs += this->numLabs*this->numLabs;
		this->numFtrFuncs += this->transMult;
		this->numTransFuncs+=1;
	}

	for (QNUInt32 clab=0; clab<this->numLabs; clab++) {
		this->stateFeatureIdxCache[clab]=this->getStateFeatureIdx(clab);
		for (QNUInt32 plab=0; plab<this->numLabs; plab++) {
			this->transFeatureIdxCache[plab*this->numLabs+clab]=this->getTransFeatureIdx(clab,plab);
		}
	}

	return this->numFtrFuncs;
}

string CRF_StdFeatureMap::getMapDescriptor(QNUInt32 lambdaNum) {
	std::stringstream ss;
	QNUInt32 stateOffset=((this->useStateFtrs)?numFeas:0)+
						 ((this->useStateBias)?1:0)+
					     ((this->useTransFtrs)?numFeas*numLabs:0)+
					     ((this->useTransBias)?numLabs:0);
	QNUInt32 stateNum=lambdaNum/stateOffset;
	QNUInt32 remainder=lambdaNum-(stateNum*stateOffset);
	if (this->useStateFtrs) {
		if (remainder < numFeas) {
			ss<< "State" << stateNum << ":ftr" << remainder;
			return ss.str();
		} else {
			remainder-=numFeas;
		}
	}
	if (this->useStateBias) {
		if (remainder == 0) {
			ss << "StateBias" << stateNum;
			return ss.str();
		} else {
			remainder--;
		}
	}
	// not sure about this!
	if (this->useTransFtrs) {
		if (remainder < numFeas*numLabs) {
			QNInt32 fromState=remainder/numFeas;
			QNInt32 ftr=remainder-(fromState*numFeas);
			ss << "Trans" << fromState << ">" << stateNum << ":ftr" << ftr;
			return ss.str();
		} else {
			remainder-=(numFeas*numLabs);
		}
	}
	if (this->useTransBias) {
		ss << "TransBias" << remainder << ">" << stateNum;
		return ss.str();
	}

	return string("Out of range");
}

void CRF_StdFeatureMap::accumulateFeatures(float *ftr_buf,double *accumulator,QNUInt32 clab) {
	int offset=clab*this->numFtrFuncs;
	for (int ii=0;ii<this->numFtrFuncs;ii++) {
		accumulator[offset+ii]+=ftr_buf[ii];
	}
}
