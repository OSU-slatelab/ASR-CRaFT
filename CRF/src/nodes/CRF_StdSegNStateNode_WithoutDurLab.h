/*
 * CRF_StdSegNStateNode_WithoutDurLab.h
 *
 *  Created on: May 2, 2012
 *      Author: Yanzhang (Ryan) He
 */

#include "CRF_StateNode.h"
#include "CRF_StdSegNStateNode.h"
#include "../io/CRF_FeatureStream.h"  // for declaration of CRF_LAB_BAD

#ifndef CRF_STDSEGNSTATENODE_WITHOUTDURLAB_H_
#define CRF_STDSEGNSTATENODE_WITHOUTDURLAB_H_

class CRF_StdSegNStateNode_WithoutDurLab : public CRF_StdSegNStateNode
{
protected:
	double* alphaArray_WithDur;

public:
	CRF_StdSegNStateNode_WithoutDurLab(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs);
	virtual ~CRF_StdSegNStateNode_WithoutDurLab();

	virtual double computeTransMatrix();
	virtual double computeAlpha();
	virtual double computeFirstAlpha();
	virtual double computeBeta(double scale=1.0);
////	virtual void setTailBeta();
	virtual double computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 prev_lab);
//	virtual double* computeAlphaBeta(double Zx);
//	virtual double computeAlphaSum();
//	virtual double computeAlphaAlignedSum();

//	virtual QNUInt32 getNActualLabs();
//	virtual bool checkNumPrevNodes();
//	virtual bool checkNumNextNodes();
//
//	// Override these functions from CRF_StdNStateNode which do not work for CRF_StdSegNStateNode.
//	virtual double computeAlpha(double* prev_alpha);
//	virtual double computeFirstAlpha(double* prev_alpha);
//	virtual double computeBeta(double* result_beta, double scale=1.0);
//	virtual double computeExpF(double* ExpF, double* grad, double Zx, double* prev_alpha, QNUInt32 prev_lab);

	// These are the correct versions for CRF_StdSegNStateNode_WithoutDurLab.
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 cur_dur);
	virtual double getStateValue(QNUInt32 cur_lab, QNUInt32 cur_dur);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 cur_dur);

	// Disable all these functions by overriding them with exception handling. Use their modified versions above.
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);
	virtual double getStateValue(QNUInt32 cur_lab);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);

	virtual double getTempBeta(QNUInt32 next_lab, QNUInt32 next_dur);
};

#endif /* CRF_STDSEGNSTATENODE_WITHOUTDURLAB_H_ */
