/*
 * CRF_StdSegStateNode.h
 *
 *  Created on: Sep 13, 2011
 *      Author: hey
 */

#ifndef CRF_STDSEGSTATENODE_H_
#define CRF_STDSEGSTATENODE_H_

#include "CRF_StateNode.h"
#include "CRF_StdStateNode.h"

class CRF_StdSegStateNode : public CRF_StdStateNode
{
protected:
	QNUInt32 labMaxDur;
	QNUInt32 nodeLabMaxDur;
	QNUInt32 nActualLabs;
	QNUInt32 nFtrsPerSeg;
	QNUInt32 prevNodeNLabs;
	QNUInt32 nextNodeNActualLabs;

public:
	CRF_StdSegStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs);
	virtual ~CRF_StdSegStateNode();
	virtual double computeTransMatrix();
	virtual double computeAlpha();
	virtual double computeFirstAlpha();
	virtual double computeBeta(double scale=1.0);
//	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 prev_lab);
	virtual double* computeAlphaBeta(double Zx);
	virtual double computeAlphaSum();
	virtual double computeAlphaAlignedSum();
//	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
//	virtual double getStateValue(QNUInt32 cur_lab);
//	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);

	virtual QNUInt32 getNActualLabs();
	virtual bool checkNumPrevNodes();
	virtual bool checkNumNextNodes();

	// Override these functions from CRF_StdStateNode but do not work for CRF_StdSegStateNode.
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
};

#endif /* CRF_STDSEGSTATENODE_H_ */
