/*
 * CRF_StdSegStateNode.h
 *
 *  Created on: Sep 13, 2011
 *      Author: hey
 */

#ifndef CRF_STDSEGSTATENODE_H_
#define CRF_STDSEGSTATENODE_H_

#include "CRF_StateNode.h"

class CRF_StdSegStateNode : public CRF_StdStateNode
{
protected:
	QNUInt32 labMaxDur;
	QNUInt32 nodeLabMaxDur;
	QNUInt32 nFullLabs;
	QNUInt32 nFtrsPerSeg;

public:
	CRF_StdSegStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 curMaxDur);
	virtual ~CRF_StdSegStateNode();
	virtual double computeTransMatrix();
	virtual double computeAlpha(double** prev_alpha);
	virtual double computeFirstAlpha(double** prev_alpha);
	virtual double computeBeta(double** result_beta, double scale=1.0);
//	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double** prev_alpha, QNUInt32 prev_lab);
	virtual double* computeAlphaBeta(double Zx);
//	virtual double computeAlphaSum();
//	virtual double computeAlphaAlignedSum();
//	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
//	virtual double getStateValue(QNUInt32 cur_lab);
//	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);

	// Override these functions which come from CRF_StdStateNode but do not work for CRF_StdSegStateNode.
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
};

#endif /* CRF_STDSEGSTATENODE_H_ */
