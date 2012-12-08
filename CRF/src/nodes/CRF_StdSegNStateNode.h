/*
 * CRF_StdSegNStateNode.h
 *
 *  Created on: Apr 29, 2012
 *      Author: hey
 */

#ifndef CRF_STDSEGNSTATENODE_H_
#define CRF_STDSEGNSTATENODE_H_

#include "CRF_StateNode.h"
#include "CRF_StdNStateNode.h"

class CRF_StdSegNStateNode : public CRF_StdNStateNode
{
protected:
	QNUInt32 labMaxDur;
	QNUInt32 nodeLabMaxDur;
	QNUInt32 nActualLabs;  // nActualLabs = nLabs / labMaxDur. TODO: change the name to make it more clear and distinct from nFullLabs and nLabs
	QNUInt32 nFtrsPerSeg;
	QNUInt32 prevNodeNLabs;
	QNUInt32 nextNodeNActualLabs;

public:
	CRF_StdSegNStateNode(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs);
	virtual ~CRF_StdSegNStateNode();
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

	// Override these functions from CRF_StdNStateNode which do not work for CRF_StdSegNStateNode.
	virtual double computeAlpha(double* prev_alpha);
	virtual double computeFirstAlpha(double* prev_alpha);
	virtual double computeBeta(double* result_beta, double scale=1.0);
	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);
};

#endif /* CRF_STDSEGNSTATENODE_H_ */
