#include "CRF_StdStateNode.h"

CRF_StdStateNode::CRF_StdStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf)
	: CRF_StateNode(fb, sizeof_fb, lab, crf)
{
	QNUInt32 nLabs=this->crf_ptr->getNLabs();
	this->stateArray = new double[nLabs];
	this->transMatrix = new double[nLabs*nLabs];
	this->alphaArray = new double[nLabs];
	this->betaArray = new double[nLabs];
	this->alphaBetaArray = new double[nLabs];
	this->tempBeta = new double[nLabs];
	this->alphaSize = nLabs;
	this->alphaScale = 0.0;
}

CRF_StdStateNode::~CRF_StdStateNode()
{
	delete [] this->stateArray;
	delete [] this->transMatrix;
	delete [] this->alphaArray;
	delete [] this->betaArray;
	delete [] this->tempBeta;
}

double CRF_StdStateNode::computeTransMatrix()
{
	double result = this->computeTransMatrixLog();
	QNUInt32 nLabs=this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		try {
			this->stateArray[clab]=expE(this->stateArray[clab]);
		}
		catch (exception &e) {
			string errstr="CRF_StdStateNode::computeTransMatrix() caught exception "+string(e.what())+" while taking exp of cell "+stringify(clab);
			throw runtime_error(errstr);
		}
		for (QNUInt32 plab=0; plab<nLabs; plab++) {
			try {
				QNUInt32 idx=plab*nLabs+clab;
				this->transMatrix[idx]=expE(this->transMatrix[idx]);
			}
			catch (exception &e) {
				string errstr="CRF_StdStateNode::computeTransMatrix() caught exception "+string(e.what())+" while taking exp of cell "+stringify(plab)+", "+stringify(clab);
				throw runtime_error(errstr);
			}
		}
	}
	return result;
}


double CRF_StdStateNode::computeTransMatrixLog()
{
	double result=0.0;
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->stateArray[clab]=0;
		for (QNUInt32 plab=0; plab<nLabs; plab++) {
			this->transMatrix[plab*nLabs+clab]=0;
		}
	}
	QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf_ptr->getLambda();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->stateArray[clab]=this->crf_ptr->getFeatureMap()->computeRi(this->ftrBuf,lambda,lc,clab);
		for (QNUInt32 plab=0; plab<nLabs; plab++) {
			QNUInt32 idx=plab*nLabs+clab;
			this->transMatrix[idx]=this->crf_ptr->getFeatureMap()->computeMij(this->ftrBuf,lambda,lc,plab,clab);
		}
	}
	return result;
}

double CRF_StdStateNode::computeAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=0.0; // Reset alphaArray to zero before processing
	}
	// cblas_sgemm(Row|Column, TransA?, TransB?, M, N, K, scale_A, A, lda, B, ldb, scale_B, C, ldc)
	// Where M is the rows of A (if TransA?==no)
	//       N is the columns of B (if TransB?==no)
	//       K is the columns of A and the rows of B (if TransA?==no && TransB? == no)
	//       lda is the linear depth of A (how many elements you must stride to get to the end of a row
	//       ldb is the linear depth of B
	//       ldc is the linear depth of C
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, this->num_labs, this->num_labs, 1.0f, alpha, this->num_labs, Mtrans, this->num_labs, 1.0f, new_alpha,this->num_labs);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, nLabs, nLabs, 1.0f, prev_alpha, nLabs, this->transMatrix, nLabs, 1.0f, this->alphaArray,nLabs);

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=this->alphaArray[clab]*this->stateArray[clab];
		this->alphaScale=this->alphaScale+this->alphaArray[clab];
	}
	if ( (this->alphaScale<1.0) && (this->alphaScale>-1.0)) this->alphaScale=1.0; //Don't inflate scores
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=this->alphaArray[clab]/this->alphaScale;
	}

	return this->alphaScale;

}

double CRF_StdStateNode::computeFirstAlpha(double* prev_alpha)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	this->alphaScale=0.0;
	//for (QNUInt32 clab=0; clab<nLabs; clab++) {
	//	this->alphaArray[clab]=0.0; // Reset alphaArray to zero before processing
	//}

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=this->stateArray[clab];
		this->alphaScale=this->alphaScale+this->alphaArray[clab];
	}
	if ( (this->alphaScale<1.0) && (this->alphaScale>-1.0)) this->alphaScale=1.0; //Don't inflate scores
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaArray[clab]=this->alphaArray[clab]/this->alphaScale;
	}
	return this->alphaScale;

}


double CRF_StdStateNode::computeBeta(double* result_beta, double scale)
{
	// Logic desired:
	//	* Compute beta_i[size of alpha[]+1] to be all 1s
	//	* Multiply M_i[current] by beta_i[current+1] to get beta_i[current]

	QNUInt32 nLabs = this->crf_ptr->getNLabs();
		//for (QNUInt32 clab=0; clab<nLabs; clab++) {
		//	this->betaArray[clab]=0.0;
		//}
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		//this->betaArray[clab]=0.0;
		this->tempBeta[clab]=(this->betaArray[clab]*this->stateArray[clab])/scale;
		result_beta[clab]=0;
	}
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, this->num_labs, 1, this->num_labs, 1.0f, Mtrans, this->num_labs, this->tmp_beta, this->num_labs, 1.0f, new_beta,1);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLabs, 1, nLabs, 1.0f, this->transMatrix, nLabs, this->tempBeta, nLabs, 1.0f, result_beta,1);
	return this->alphaScale;

}

double* CRF_StdStateNode::computeAlphaBeta(double Zx)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->alphaBetaArray[clab]=this->alphaArray[clab]*this->betaArray[clab]*this->alphaScale/Zx;
	}
	return this->alphaBetaArray;
}

void CRF_StdStateNode::setTailBeta()
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		this->betaArray[clab]=1.0/this->alphaScale;
	}
}

double CRF_StdStateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	double logLi=0.0;
	double alpha_beta=0.0;
	QNUInt32 nLabs = this->crf_ptr->getNLabs();

	QNUInt32 lc=0;  // Counter for weight position in crf->lambda[] array...
	double* lambda = this->crf_ptr->getLambda();
	double alpha_beta_tot = 0.0;
	double alpha_beta_trans_tot=0.0;

	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		alpha_beta=this->alphaArray[clab]*this->betaArray[clab]*this->alphaScale/Zx;
		alpha_beta_tot += alpha_beta;
		bool match=(clab==this->label);
		logLi+=this->crf_ptr->getFeatureMap()->computeExpFState(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab);
		if (prev_lab > nLabs) {
			// if prev_lab > nLabs, we're in the first label frame and there are no previous
			// transitions - skip the transition calculation in this case
			// but set the alpha_beta_trans_tot to 1.0 for the check below
			alpha_beta_trans_tot=1.0;
			// We STILL need to cycle through our "lc" values because of how the FeatureMap code works
			// Our feature map is expecting us to have lc point at the next input vector
			//double* tmp_ExpF = new double[this->crf_ptr->getLambdaLen()];
			for (QNUInt32 plab=0; plab<nLabs; plab++) {
				lc+=this->crf_ptr->getFeatureMap()->getNumTransFuncs(plab,clab);
				//this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,NULL,lc,tmp_ExpF,NULL,alpha_beta,false,plab,clab);
			}
			//delete tmp_ExpF;
		}
		else {
			// if prev_lab > nLabs, we're in the first label frame and there are no previous
			// transitions - skip the transition calculation in this case
			for (QNUInt32 plab=0; plab<nLabs; plab++) {
				QNUInt32 idx = plab*nLabs+clab;
				alpha_beta=prev_alpha[plab]*this->transMatrix[idx]*this->stateArray[clab]*this->betaArray[clab]/Zx;
				alpha_beta_trans_tot+=alpha_beta;
				match=((clab==this->label)&&(plab==prev_lab));
				//logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,clab,plab);
				logLi+=this->crf_ptr->getFeatureMap()->computeExpFTrans(this->ftrBuf,lambda,lc,ExpF,grad,alpha_beta,match,plab,clab);
			}
		}
	}

	if ((alpha_beta_tot >1.1))  {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Probability sums greater than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_tot < 0.9) {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Probability sums less than 1.0 "+stringify(alpha_beta_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot > 1.1) {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Trans Probability sums greater than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	else if (alpha_beta_trans_tot < 0.9) {
		string errstr="CRF_StdStateNode::computeExpF() threw exception: Trans Probability sums less than 1.0 "+stringify(alpha_beta_trans_tot);
		throw runtime_error(errstr);
	}
	return logLi;

}

double CRF_StdStateNode::computeAlphaSum()
{
	double Zx = 0.0;
	QNUInt32 nLabs=this->crf_ptr->getNLabs();
	for (QNUInt32 clab=0; clab<nLabs; clab++) {
		Zx+=this->alphaArray[clab];
	}
	return Zx;
}

double CRF_StdStateNode::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	return this->transMatrix[prev_lab*nLabs+cur_lab];
}

double CRF_StdStateNode::getStateValue(QNUInt32 cur_lab)
{
	return this->stateArray[cur_lab];
}

double CRF_StdStateNode::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	QNUInt32 nLabs = this->crf_ptr->getNLabs();
	return this->transMatrix[prev_lab*nLabs+cur_lab]*this->stateArray[cur_lab];
}
