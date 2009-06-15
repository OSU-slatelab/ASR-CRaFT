#include "CRF_StdTransFeatureMap.h"

CRF_StdTransFeatureMap::CRF_StdTransFeatureMap(QNUInt32 nlabs, QNUInt32 nfeas)
	: CRF_StdFeatureMap(nlabs,nfeas)
{
	this->numFtrFuncs=nlabs*nfeas + nlabs + nlabs*nlabs*nfeas + nlabs*nlabs;
	this->transFidxStart=0;
	this->transFidxEnd=nfeas-1;
}

CRF_StdTransFeatureMap::~CRF_StdTransFeatureMap()
{
}

double CRF_StdTransFeatureMap::computeMij(float* ftr_buf, double* lambda, QNUInt32& lc, QNUInt32 plab, QNUInt32 clab)
{
	double Mi=0.0;
//	for (QNUInt32 fidx=0; fidx<this->numFeas; fidx++) 
	for (QNUInt32 fidx=this->transFidxStart; fidx<=this->transFidxEnd; fidx++)
	{
		Mi+=ftr_buf[fidx]*lambda[lc];
		lc++;
	}
	if (this->useTransBias) {
		Mi+=lambda[lc];
		lc++;
	}
	return Mi;
}

double CRF_StdTransFeatureMap::computeExpFTrans(float* ftr_buf, double* lambda, QNUInt32& lc, double* ExpF, double* grad, double alpha_beta, bool match, QNUInt32 plab, QNUInt32 clab)
{
	double logLi=0.0;
//	for (QNUInt32 fidx=0; fidx<this->numFeas; fidx++)
	for (QNUInt32 fidx=this->transFidxStart; fidx<=this->transFidxEnd; fidx++)
	{		
		ExpF[lc]+=alpha_beta*ftr_buf[fidx];
		if (match) {
			grad[lc]+=ftr_buf[fidx];
			logLi += lambda[lc]*ftr_buf[fidx];
		}
		lc++;
	} 
	if (this->useTransBias) {
		ExpF[lc]+=alpha_beta;
		if (match) {
			grad[lc]+=1;
			logLi += lambda[lc];
		}
		lc++;
	}
	return logLi;
}	

void CRF_StdTransFeatureMap::setTransFtrStart(QNUInt32 st)
{
	this->transFidxStart=st;
}

void CRF_StdTransFeatureMap::setTransFtrEnd(QNUInt32 end)
{
	this->transFidxEnd=end;
}	

QNUInt32 CRF_StdTransFeatureMap::init()
{
	this->numFtrFuncs=(this->stateFidxEnd - this->stateFidxStart +1 )*this->numLabs; // Count start at zero
	if (this->useStateBias) {
		this->numFtrFuncs += this->numLabs;
	}
	this->numFtrFuncs+=(this->transFidxEnd - this->transFidxStart + 1 )*this->numLabs*this->numLabs; //Count still starts at zero
	if (this->useTransBias) {
		this->numFtrFuncs += this->numLabs*this->numLabs;
	}
	return this->numFtrFuncs;
	
}
