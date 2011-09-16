/*
 * CRF_StateNode.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */

#include "CRF_StateNode.h"

#include "CRF_StdStateNode.h"
#include "CRF_StdNStateNode.h"

/*
 * CRF_StateNode constructor
 *
 * Input: fb - feature buffer for features applicable for this node
 *        sizeof_fb - number of features in buffer fb
 *        lab - label for this node
 *        crf_in - pointer back to CRF model used by this node
 */
CRF_StateNode::CRF_StateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf_in)
	: ftrBuf(fb),
	  ftrBuf_size(sizeof_fb),
	  label(lab),
	  crf_ptr(crf_in)
{
	this->nLabs=this->crf_ptr->getNLabs();
}

/*
 * cRF_StateNode destructor
 */
CRF_StateNode::~CRF_StateNode()
{
	delete [] this->ftrBuf;
}

/*
 * CRF_StateNode::computeTransMatrix
 *
 * Stub function.
 * Should compute the transition matrix for this node based on the crf model and the feature buffer.
 */
double CRF_StateNode::computeTransMatrix()
{
	return 0;
}

/*
 * CRF_StateNode::computeTransMatrixLog
 *
 * Stub function.
 * Should compute the transition matrix in logspace for this node based on the crf model and the
 * feature buffer.
 */
double CRF_StateNode::computeTransMatrixLog()
{
	return 0;
}

/*
 * CRF_StateNode::computeAlpha
 *
 * Input: prev_alpha - alpha vector of the previous node
 *
 * Stub function.
 * Should compute the alpha vector for the forward backward computation for this node.
 */
double CRF_StateNode::computeAlpha(double* prev_alpha)
{
	return 0;
}

// Added by Ryan
/*
 * CRF_StateNode::computeAlpha
 *
 * Input: prev_alpha - vector of alpha vectors of a few previous nodes
 *
 * Stub function.
 * Should compute the alpha vector for the forward backward computation for this node.
 */
double CRF_StateNode::computeAlpha(double** prev_alpha)
{
	return 0;
}


/*
 * CRF_StateNode::computeFirstAlpha
 *
 * Stub function.
 * Should compute the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StateNode::computeFirstAlpha(double* prev_alpha)
{
	return this->computeAlpha(prev_alpha);
}

// Added by Ryan
/*
 * CRF_StateNode::computeFirstAlpha
 *
 * This version is to match the input type of computeAlpha(double** prev_alpha).
 *
 * Stub function.
 * Should compute the alpha vector for this node for the special case where the node is the first
 * node in the sequence.
 */
double CRF_StateNode::computeFirstAlpha(double** prev_alpha)
{
	return this->computeAlpha(prev_alpha);
}

/*
 * CRF_StateNode::computeBeta
 *
 * Inputs: result_beta - vector to store the result of the computation
 *         scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Stub function.
 * Should compute the beta vector for the node before this one and store it in result_beta
 */
double CRF_StateNode::computeBeta(double* result_beta, double scale)
{
	return 0;
}

// Added by Ryan
/*
 * CRF_StateNode::computeBeta
 *
 * Inputs: result_beta - vector of vectors to store the result of beta vectors of a few previous nodes
 *         scale - scaling constant for result_beta array
 *
 * Returns:
 *
 * Stub function.
 * Should compute the beta vector for the node before this one and store it in result_beta
 */
double CRF_StateNode::computeBeta(double** result_beta, double scale)
{
	return 0;
}

/*
 * CRF_StateNode::computeAlphaBeta
 *
 * Inputs: Zx - normalization constant
 *
 * Returns: array of state probabilities based on the precomputed alpha and beta for this node
 *
 * Stub function.
 */
double* CRF_StateNode::computeAlphaBeta(double Zx)
{
	return NULL;
}

/*
 * CRF_StateNode::setTailBeta
 *
 * Stub function.
 * Should set the beta value in this node to the special case for the end of the sequence.
 */
void CRF_StateNode::setTailBeta()
{
}

/*
 * CRF_StateNode::computeExpF
 *
 * Inputs: *ExpF - vector to store expected values of feature functions
 *         *grad - vector to store computed gradient values
 *         Zx - normalization constant
 *         *prev_alpha - vector of alpha values from the previous node (for use in transition feature
 *            ExpF computation
 *         prev_lab - previous node label (transition feature ExpF computation)
 *
 * Returns:
 *
 * Stub function.
 * Should compute gradient and expected values for features in this node and store them in *grad and
 *   *ExpF vectors respectively.  State features and transition features are computed in the same function.
 */
double CRF_StateNode::computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab)
{
	return 0;
}

/*
 * CRF_StateNode::computeExpF
 *
 * Inputs: *ExpF - vector to store expected values of feature functions
 *         *grad - vector to store computed gradient values
 *         Zx - normalization constant
 *         **prev_alpha - vector of vectors of alpha values from a few previous nodes (for use in transition feature
 *            ExpF computation
 *         prev_lab - previous node label (transition feature ExpF computation)
 *
 * Returns:
 *
 * Stub function.
 * Should compute gradient and expected values for features in this node and store them in *grad and
 *   *ExpF vectors respectively.  State features and transition features are computed in the same function.
 */
double CRF_StateNode::computeExpF(double* ExpF, double* grad, double Zx, double** prev_alpha, QNUInt32 prev_lab)
{
	return 0;
}

/*
 * CRF_StateNode::computeSoftExpF
 *
 * Inputs: *ExpF - vector to store expected values of feature functions
 *         *grad - vector to store computed gradient values
 *         Zx - normalization constant
 *         soft_Zx - normalization constant for 'soft' numerator
 *         *prev_alpha - vector of alpha values from the previous node (for use in transition feature
 *            ExpF computation
 *         *prevAlphaAligned - vector of alpha values for 'soft' numerator
 *         firstFrame - boolean indicator of whether this is the first sequence or not
 * Returns:
 *
 * Stub function.
 * Should compute gradient and expected values for features in this node and store them in *grad and
 *   *ExpF vectors respectively.  State features and transition features are computed in the same function.
 *   "Soft" means that we use a soft alignment of the label in the numerator of the gradient calculation
 *   instead of the hard assignment.
 */
double CRF_StateNode::computeSoftExpF(double* ExpF, double* grad, double Zx, double soft_Zx, double* prev_alpha, vector<double>* prevAlphaAligned, bool firstFrame)
{
	return 0;
}

/*
 * CRF_StateNode::computeAlphaSum
 *
 * Returns: Sum of the values in the alpha vector of this node
 *
 * Stub function.
 * Used to compute the normalization constant for the CRF.
 */
double CRF_StateNode::computeAlphaSum()
{
	return 0;
}

/*
 * CRF_StateNode::computeAlphaAlignedSum
 *
 * Returns: Sum of the values in the alphaAligned vector of this node
 *
 * Stub function.
 * Used to compute the "soft" normalization constant.
 */
double CRF_StateNode::computeAlphaAlignedSum()
{
	return 0;
}

/*
 * CRF_StateNode::reset
 *
 * Input: see constructor
 *
 * Clears out the contents of the node and replaces it with new contents for recalculation.
 */

void CRF_StateNode::reset(float *fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf_in)
{
	//memcpy(fb,this->ftrBuf,sizeof_fb);
	if (this->ftrBuf != NULL) { delete[] this->ftrBuf; }
	this->ftrBuf=fb;
	this->ftrBuf_size=sizeof_fb;
	this->label=lab;
	this->crf_ptr=crf_in;
	this->alphaScale=0.0;
}

/*
 *  CRF_StateNode::getAlpha
 *
 *  Returns: pointer to alpha vector
 *
 *  Accessor function for alpha vector
 */
double* CRF_StateNode::getAlpha()
{
	return this->alphaArray;
}

/*
 * CRF_StateNode::getBeta
 *
 * Returns: pointer to beta vector
 *
 * Accessor function for beta vector
 */
double* CRF_StateNode::getBeta()
{
	return this->betaArray;
}

/*
 * CRF_StateNode::getAlphaBeta
 *
 * Returns: pointer to alphabeta vector
 *
 * Accessor function for alphabeta vector
 */
double* CRF_StateNode::getAlphaBeta()
{
	return this->alphaBetaArray;
}

/*
 * CRF_StateNode::getLabel
 *
 * Returns: label associated with node
 */
QNUInt32 CRF_StateNode::getLabel()
{
	return this->label;
}

/*
 * CRF_StateNode::getAlphaScale
 *
 * Returns: scaling value for alpha vector
 */
double CRF_StateNode::getAlphaScale()
{
	return this->alphaScale;
}

/*
 * CRF_StateNode::getAlphaAligned
 *
 * Returns: aligned alpha vector for "soft" numerator processing
 */
vector<double>* CRF_StateNode::getAlphaAligned()
{
	return &(this->alphaArrayAligned);
}

/*
 * CRF_StateNode::getBetaAligned
 *
 * Returns: aligned beta vector for "soft" numerator processing
 */
vector<double>* CRF_StateNode::getBetaAligned()
{
	return &(this->betaArrayAligned);
}

/*
 * CRF_StateNode::getPrevAlpha
 *
 * Returns: alpha vector of the previous node (used for debugging)
 */
double* CRF_StateNode::getPrevAlpha()
{
	return this->prevAlpha;
}

/*
 * CRF_StateNode::getAlphaAlignedBase
 *
 * Returns: base vector for "soft" alpha computation
 *
 */
vector<double>* CRF_StateNode::getAlphaAlignedBase()
{
	return &(this->alphaArrayAlignedBase);
}

/*
 * CRF_StateNode::getBetaAlignedBase
 *
 * Returns: base vector for "soft" beta computation
 */
vector<double>* CRF_StateNode::getBetaAlignedBase()
{
	return &(this->betaArrayAlignedBase);
}

/*
 * CRF_StateNode::getTransValue
 *
 * Returns: Transition matrix value for transition prev_lab->cur_lab
 */
double CRF_StateNode::getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	return 0.0;
}

/*
 * CRF_StateNode::getStateValue
 *
 * Returns: State vector value for state label cur_lab
 */
double CRF_StateNode::getStateValue(QNUInt32 cur_lab)
{
	return 0.0;
}


/*
 * CRF_StateNode::getFullTransValue
 *
 * Returns: Full value for transition prev_lab->cur_lab
 *   (i.e. transition value prev_lab->cur_lab plus the state value for cur_lab)
 */

double CRF_StateNode::getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab)
{
	return 0.0;
}

/*
 *  CRF_StateNode::createStateNode
 *
 *  Input: see constructor
 *
 *  Factory class: This depends on the fact that the StateVector saves its nodes between
 *   each pass.  The calls on this function are bounded by the size of the longest sequence
 *   being examined.
 *
 */

CRF_StateNode* CRF_StateNode::createStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf) {

	if (crf->getFeatureMap()->getNumStates()>1) {
		return new CRF_StdNStateNode(fb, sizeof_fb, lab, crf);
	}
	else {
		return new CRF_StdStateNode(fb, sizeof_fb, lab, crf);
	}
}

/*
 * CRF_StateNode::getFtrBuffer
 *
 * Returns: Feature buffer for this node
 */

float *CRF_StateNode::getFtrBuffer() {
	return this->ftrBuf;
}

/*
 * CRF_StateNode::getFtrBufferSize
 *
 * Returns: size of feature buffer for this node
 */

QNUInt32 CRF_StateNode::getFtrBufferSize() {
	return this->ftrBuf_size;
}
