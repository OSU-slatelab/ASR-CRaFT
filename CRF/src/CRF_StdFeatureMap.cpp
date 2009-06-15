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
	this->recalc();
}

CRF_StdFeatureMap::~CRF_StdFeatureMap()
{
}

double CRF_StdFeatureMap::computeRi(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 clab)
{
	double Ri=0.0;
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=this->stateFidxStart; fidx<=this->stateFidxEnd; fidx++)
		{
			Ri+=ftr_buf[fidx]*lambda[lc];
			lc++;
		}
	}
	if (this->useStateBias) {
		Ri+=lambda[lc]*this->stateBiasVal;
		lc++;
	}
	return Ri;
}

double CRF_StdFeatureMap::computeMij(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 plab, QNUInt32 clab)
{
	double Mi=0.0;
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=this->transFidxStart; fidx<=this->transFidxEnd; fidx++)
		{
			Mi+=ftr_buf[fidx]*lambda[lc];
			lc++;
		}
	}
	if (this->useTransBias) {
		Mi+=lambda[lc]*this->transBiasVal;
		lc++;
	}
	return Mi;
}

double CRF_StdFeatureMap::computeExpFState(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab)
{
	double logLi=0.0;
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=this->stateFidxStart; fidx<=this->stateFidxEnd; fidx++)
		{
			ExpF[lc]+=alpha_beta*ftr_buf[fidx];
			if (match) {
				grad[lc]+=ftr_buf[fidx];
				logLi += lambda[lc]*ftr_buf[fidx];
			}
			lc++;
		}
	}
	if (this->useStateBias) {
		ExpF[lc]+=alpha_beta*this->stateBiasVal;
		if (match) {
			grad[lc]+=this->stateBiasVal;
			logLi += lambda[lc]*this->stateBiasVal;
		}
		lc++;
	}
	return logLi;
}

double CRF_StdFeatureMap::computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab)
{
	double logLi=0.0;
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=this->transFidxStart; fidx<=this->transFidxEnd; fidx++)
		{
			ExpF[lc]+=alpha_beta*ftr_buf[fidx];
			if (match) {
				grad[lc]+=ftr_buf[fidx];
				logLi += lambda[lc]*ftr_buf[fidx];
			}
			lc++;
		}
	}
	if (this->useTransBias) {
		ExpF[lc] += alpha_beta*this->transBiasVal;
		if (match) {
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
	QNUInt32 retVal=0;
	if (this->useStateFtrs)
	{
		retVal+=this->stateFidxEnd - this->stateFidxStart +1;
	}
	if (this->useStateBias)
	{
		retVal+=1;
	}
	return retVal;
}

QNUInt32 CRF_StdFeatureMap::getNumTransFuncs(QNUInt32 plab, QNUInt32 clab)
{
	QNUInt32 retVal=0;
	if (this->useTransFtrs)
	{
		retVal+=this->transFidxEnd - this->transFidxStart +1;
	}
	if (this->useTransBias)
	{
		retVal+=1;
	}
	return retVal;
}


QNUInt32 CRF_StdFeatureMap::recalc()
{
	QNUInt32 actualLabels = this->numLabs/this->numStates;
	cout << "ACTUAL LABELS COMPUTED: " << actualLabels << endl;
	if (actualLabels * this->numStates != this->numLabs) {
		string errstr="CRF_StdFeatureMap created exception: Invalid state/label combination while computing transitions";
		throw runtime_error(errstr);
	}
	QNUInt32 transMult;
	if (this->numStates==1) {
		transMult=actualLabels*actualLabels;
	}
	else {
		transMult = actualLabels*actualLabels+this->numLabs + this->numLabs-actualLabels;
	}
	// end->start transitions + diagonal self transitions + offDiagonal transitions
	this->numFtrFuncs=0;
	if (this->useStateFtrs) {
		this->numFtrFuncs+=(this->stateFidxEnd - this->stateFidxStart + 1 )*this->numLabs; // Count start at zero
	}
	if (this->useStateBias) {
		this->numFtrFuncs += this->numLabs;
	}
	if (this->useTransFtrs) {
		//this->numFtrFuncs+=(this->transFidxEnd - this->transFidxStart + 1 )*this->numLabs*this->numLabs; //Count still starts at zero
		this->numFtrFuncs+=(this->transFidxEnd - this->transFidxStart + 1 )*transMult; //Count still starts at zero
	}
	if (this->useTransBias) {
		//this->numFtrFuncs += this->numLabs*this->numLabs;
		this->numFtrFuncs += transMult;
	}
	return this->numFtrFuncs;
}

