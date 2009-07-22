#include "CRF_StdSparseFeatureMap.h"

CRF_StdSparseFeatureMap::CRF_StdSparseFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas)
	: CRF_StdFeatureMap(nlabs,nfeas)
{
	this->recalc();
}

CRF_StdSparseFeatureMap::~CRF_StdSparseFeatureMap()
{
}

// In all of the following functions, lc contains a pointer to the last observed lambda value, plus one (i.e. the "current"
//  lambda value if we were going sequentially).  To maintain compatibility, each function should exit with lc equal to the
//  "next" lambda value to be observed in sequence.

double CRF_StdSparseFeatureMap::computeRi(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 clab)
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

double CRF_StdSparseFeatureMap::computeMij(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 plab, QNUInt32 clab)
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

double CRF_StdSparseFeatureMap::computeStateArrayValue(float* ftr_buf, double* lambda, QNUInt32 clab)
{
	double stateValue=0.0;
	QNUInt32 lc=this->stateFeatureIdxCache[clab];
	//QNUInt32 tmp_lc;
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=0; fidx<this->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= this->stateFidxStart) && (new_idx <= this->stateFidxEnd)) {
				//tmp_lc=lc+new_idx;
				//cout << "S: Computing for new_idx: " << new_idx << " value " << ftr_buf[fidx+1] << " tmp_lc " << tmp_lc << endl;
				stateValue+=ftr_buf[fidx+1]*lambda[lc+new_idx];
			}
		}
	}
	if (this->useStateBias) {
		//tmp_lc=lc+this->numStateFuncs-1;
		//cout << "SB: Computing for new_idx: XX" << " value " << 1 << " tmp_lc " << tmp_lc << " C: " << clab << endl;
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs-1; //Advance to end of state functions, then backup
		stateValue+=lambda[lc+this->numStateFuncs-1];
	}
	return stateValue;
}

double CRF_StdSparseFeatureMap::computeTransMatrixValue(float* ftr_buf, double* lambda, QNUInt32 plab, QNUInt32 clab)
{
	// Here, lc is the start of the trans feature weights for the combination of plab,clab
	double transMatrixValue=0.0;
	QNUInt32 lc=this->transFeatureIdxCache[plab*this->numLabs+clab];
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=0; fidx<this->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= this->transFidxStart) && (new_idx<=this->transFidxEnd)) {
				//cout << "T: Computing for new_idx: " << new_idx << " value " << ftr_buf[fidx+1] << " tmp_lc " << tmp_lc << endl;
				//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+plab*this->numTransFuncs+new_idx;
				transMatrixValue+=ftr_buf[fidx+1]*lambda[lc+new_idx];
				//lc++;
			}
		}
	}
	if (this->useTransBias) {
		//lc=clab*(this->numStateFuncs+this->numTransFuncs)+this->numStateFuncs+(plab+1)*this->numTransFuncs-1;
		//cout << "TB: Computing for new_idx: XX"  << " value " << 1 << " tmp_lc " << tmp_lc << " P: " << plab << " C: "<< clab << endl;
		transMatrixValue+=lambda[lc+this->numTransFuncs-1];
		//lc++;
	}
	return transMatrixValue;
}

double CRF_StdSparseFeatureMap::computeExpFState(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 clab)
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

double CRF_StdSparseFeatureMap::computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab)
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

double CRF_StdSparseFeatureMap::computeStateExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_clab, QNUInt32 clab)
{
	double logLi=0.0;
	QNUInt32 lc = this->stateFeatureIdxCache[clab];
	QNUInt32 tmp_lc;
	if (this->useStateFtrs) {
		for (QNUInt32 fidx=0; fidx<this->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= this->stateFidxStart) && (new_idx <= this->stateFidxEnd)) {
				tmp_lc=lc+new_idx;
				ExpF[tmp_lc]+=alpha_beta*ftr_buf[fidx+1];
				if (t_clab == clab) {
					grad[tmp_lc]+=ftr_buf[fidx+1];
					logLi += lambda[tmp_lc]*ftr_buf[fidx+1];
				}
			}
		}
	}
	if (this->useStateBias) {
		tmp_lc=lc+this->numStateFuncs-1;
		ExpF[tmp_lc]+=alpha_beta;
		if (t_clab==clab) {
			grad[tmp_lc]+=1;
			logLi += lambda[tmp_lc];
		}
	}
	return logLi;
}

double CRF_StdSparseFeatureMap::computeTransExpF(float* ftr_buf, double* lambda, double* ExpF, double* grad, double alpha_beta, QNUInt32 t_plab, QNUInt32 t_clab, QNUInt32 plab, QNUInt32 clab)
{
	double logLi=0.0;
	QNUInt32 lc = this->transFeatureIdxCache[plab*this->numLabs+clab];
	QNUInt32 tmp_lc;
	if (this->useTransFtrs) {
		for (QNUInt32 fidx=0; fidx<this->numFeas; fidx=fidx+2)
		{
			QNUInt32 new_idx=(QNUInt32) ftr_buf[fidx];
			if ((new_idx >= this->transFidxStart) && (new_idx<=this->transFidxEnd)) {
				tmp_lc=lc+new_idx;
				ExpF[tmp_lc]+=alpha_beta*ftr_buf[fidx+1];
				if ((clab==t_clab) && (plab==t_plab)) {
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
		if ((clab==t_clab) && (plab==t_plab)) {
			grad[tmp_lc]+=1;
			logLi+=lambda[tmp_lc];
		}
	}
	return logLi;
}


// Below removed because it is now exactly the same as the code in StdFeatureMap
/*QNUInt32 CRF_StdSparseFeatureMap::recalc()
{
	this->numActualLabels = this->numLabs/this->numStates;
	cout << "ACTUAL LABELS COMPUTED: " << this->numActualLabels << endl;
	if (this->numActualLabels * this->numStates != this->numLabs) {
		string errstr="CRF_StdFeatureMap created exception: Invalid state/label combination while computing transitions";
		throw runtime_error(errstr);
	}
	QNUInt32 transMult;
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
}*/

