#include "CRF_StdSparseFeatureMapOld.h"

CRF_StdSparseFeatureMapOld::CRF_StdSparseFeatureMapOld(QNUInt32 nlabs, QNUInt32 nfeas)
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
	this->recalc();
}

CRF_StdSparseFeatureMapOld::~CRF_StdSparseFeatureMapOld()
{
}

// In all of the following functions, lc contains a pointer to the last observed lambda value, plus one (i.e. the "current"
//  lambda value if we were going sequentially).  To maintain compatibility, each function should exit with lc equal to the
//  "next" lambda value to be observed in sequence.

double CRF_StdSparseFeatureMapOld::computeRi(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 clab)
{
	double Ri=0.0;
	QNUInt32 tmp_lc;
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=0; fidx<this->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= this->stateFidxStart) && (new_idx <= this->stateFidxEnd)) {
				tmp_lc=lc+new_idx;
				//cout << "S: Computing for new_idx: " << new_idx << " value " << ftr_buf[fidx+1] << " tmp_lc " << tmp_lc << endl;
				Ri+=ftr_buf[fidx+1]*lambda[tmp_lc];
			}
		}
	}
	if (this->useStateBias) {
		tmp_lc=lc+this->numStateFuncs-1;
		//cout << "SB: Computing for new_idx: XX" << " value " << 1 << " tmp_lc " << tmp_lc << " C: " << clab << endl;
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs-1; //Advance to end of state functions, then backup
		Ri+=lambda[tmp_lc];
	}
	lc+=this->numStateFuncs;
	return Ri;
}

double CRF_StdSparseFeatureMapOld::computeMij(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 plab, QNUInt32 clab)
{
	// Here, lc is the start of the trans feature weights for the combination of plab,clab
	double Mi=0.0;
	QNUInt32 tmp_lc;
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=0; fidx<this->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= this->transFidxStart) && (new_idx<=this->transFidxEnd)) {
				tmp_lc=lc+new_idx;
				//cout << "T: Computing for new_idx: " << new_idx << " value " << ftr_buf[fidx+1] << " tmp_lc " << tmp_lc << endl;
				//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+plab*this->numTransFuncs+new_idx;
				Mi+=ftr_buf[fidx+1]*lambda[tmp_lc];
				//lc++;
			}
		}
	}
	if (this->useTransBias) {
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+(plab+1)*this->numTransFuncs-1;
		tmp_lc=lc+this->numTransFuncs-1;
		//cout << "TB: Computing for new_idx: XX"  << " value " << 1 << " tmp_lc " << tmp_lc << " P: " << plab << " C: "<< clab << endl;
		Mi+=lambda[tmp_lc];
		//lc++;
	}
	lc+=this->numTransFuncs;
	return Mi;
}

double CRF_StdSparseFeatureMapOld::computeExpFState(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab)
{
	double logLi=0.0;
	QNUInt32 tmp_lc;
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=0; fidx<this->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= this->stateFidxStart) && (new_idx <= this->stateFidxEnd)) {
				tmp_lc=lc+new_idx;
				//lc=clab*(this->numStateFuncs+this->numTransFuncs)+new_idx;
				ExpF[tmp_lc]+=alpha_beta*ftr_buf[fidx+1];
				if (match) {
					grad[tmp_lc]+=ftr_buf[fidx+1];
					logLi += lambda[tmp_lc]*ftr_buf[fidx+1];
				}
				//lc++;
			}
		}
	}
	if (this->useStateBias) {
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs-1; //Advance to end of state functions, then backup
		tmp_lc=lc+this->numStateFuncs-1;
		ExpF[tmp_lc]+=alpha_beta;
		if (match) {
			grad[tmp_lc]+=1;
			logLi += lambda[tmp_lc];
		}
		//lc++;
	}
	lc+=this->numStateFuncs;
	return logLi;
}

double CRF_StdSparseFeatureMapOld::computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab)
{
	double logLi=0.0;
	QNUInt32 tmp_lc;
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=0; fidx<this->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= this->transFidxStart) && (new_idx<=this->transFidxEnd)) {
				//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+plab*this->numTransFuncs+new_idx;
				tmp_lc=lc+new_idx;
				ExpF[tmp_lc]+=alpha_beta*ftr_buf[fidx+1];
				if (match) {
					grad[tmp_lc]+=ftr_buf[fidx+1];
					logLi += lambda[tmp_lc]*ftr_buf[fidx+1];
				}
				//lc++;
			}
		}
	}
	if (this->useTransBias) {
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+(plab+1)*this->numTransFuncs-1;
		tmp_lc=lc+this->numTransFuncs-1;
		ExpF[tmp_lc] += alpha_beta;
		if (match) {
			grad[tmp_lc]+=1;
			logLi+=lambda[tmp_lc];
		}
		//lc++;
	}
	lc+=this->numTransFuncs;
	return logLi;
}

void CRF_StdSparseFeatureMapOld::setStateFtrRange(QNUInt32 st, QNUInt32 end)
{
	this->stateFidxStart=st;
	this->stateFidxEnd=end;
	if (this->stateFidxEnd - this->stateFidxStart +1 >0 ) {
		this->useStateFtrs=true;
	}
	this->recalc();
}

void CRF_StdSparseFeatureMapOld::setTransFtrRange(QNUInt32 st, QNUInt32 end)
{
	this->transFidxStart=st;
	this->transFidxEnd=end;
	if (this->transFidxEnd - this->transFidxStart +1 > 0) {
		this->useTransFtrs=true;
	}
	this->recalc();
}


void CRF_StdSparseFeatureMapOld::setNumStates(QNUInt32 ns)
{
	this->numStates=ns;
	this->recalc();
}

void CRF_StdSparseFeatureMapOld::setUseStateBias(bool useState)
{
	this->useStateBias=useState;
	this->recalc();
}

void CRF_StdSparseFeatureMapOld::setUseTransBias(bool useTrans)
{
	this->useTransBias=useTrans;
	this->recalc();
}

void CRF_StdSparseFeatureMapOld::setUseStateFtrs(bool useState)
{
	this->useStateFtrs = useState;
	this->recalc();
}

void CRF_StdSparseFeatureMapOld::setUseTransFtrs(bool useTrans)
{
	this->useTransFtrs = useTrans;
	if (this->transFidxEnd == this->numLabs+1) { // By default if it isn't set, use the same features as the state features
		this->transFidxStart=this->stateFidxStart;
		this->transFidxEnd=this->stateFidxEnd;
	}
	this->recalc();
}

void CRF_StdSparseFeatureMapOld::setStateBiasVal(double stateBias)
{
	this->stateBiasVal=stateBias;
}

void CRF_StdSparseFeatureMapOld::setTransBiasVal(double transBias)
{
	this->transBiasVal=transBias;
}

QNUInt32 CRF_StdSparseFeatureMapOld::getNumStateFuncs(QNUInt32 clab)
{
	return this->numStateFuncs;
}

QNUInt32 CRF_StdSparseFeatureMapOld::getNumTransFuncs(QNUInt32 plab, QNUInt32 clab)
{
	return this->numTransFuncs;
}


QNUInt32 CRF_StdSparseFeatureMapOld::getNumStates()
{
	return this->numStates;
}

QNUInt32 CRF_StdSparseFeatureMapOld::recalc()
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
	this->numStateFuncs=0;
	this->numTransFuncs=0;
	if (this->useStateFtrs) {
		this->numFtrFuncs+=(this->stateFidxEnd - this->stateFidxStart + 1 )*this->numLabs; // Count start at zero
		this->numStateFuncs+=(this->stateFidxEnd - this->stateFidxStart +1 );
	}
	if (this->useStateBias) {
		this->numFtrFuncs += this->numLabs;
		this->numStateFuncs += 1;
	}
	if (this->useTransFtrs) {
		//this->numFtrFuncs+=(this->transFidxEnd - this->transFidxStart + 1 )*this->numLabs*this->numLabs; //Count still starts at zero
		this->numFtrFuncs+=(this->transFidxEnd - this->transFidxStart + 1 )*transMult; //Count still starts at zero
		this->numTransFuncs += (this->transFidxEnd - this->transFidxStart +1);
	}
	if (this->useTransBias) {
		//this->numFtrFuncs += this->numLabs*this->numLabs;
		this->numFtrFuncs += transMult;
		this->numTransFuncs += 1;
	}
	return this->numFtrFuncs;
}

