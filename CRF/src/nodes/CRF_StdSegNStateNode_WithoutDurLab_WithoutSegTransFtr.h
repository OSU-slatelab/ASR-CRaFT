/*
 * CRF_StdSegNStateNode_WithoutDurLab_WithoutSegTransFtr.h
 *
 *  Created on: May 2, 2012
 *      Author: Yanzhang (Ryan) He
 */

#ifndef CRF_STDSEGNSTATENODE_WITHOUTDURLAB_WITHOUTSEGTRANSFTR_H_
#define CRF_STDSEGNSTATENODE_WITHOUTDURLAB_WITHOUTSEGTRANSFTR_H_

#include "CRF_StateNode.h"
#include "CRF_StdSegNStateNode_WithoutDurLab.h"
#include "../io/CRF_FeatureStream.h"  // for declaration of CRF_LAB_BAD

class CRF_StdSegNStateNode_WithoutDurLab_WithoutSegTransFtr : public CRF_StdSegNStateNode_WithoutDurLab {

protected:
	// the array of the summation of [alpha(cur_lab) + trans_value(cur_lab, next_lab)]
	// over all cur_lab for a given next_lab,
	// for future calculation of alphas in next nodes.
	double* alphaPlusTrans;
	// the array of the summation of [next_node_beta(next_lab, dur) + next_node_state_value(next_lab, dur)]
	// over all duration for a given next_lab,
	// for calculation of betas and ExpF(trans) in the current node.
	double* tmpBetaArray_nextBetasPlusNextStateValue_sumOverDur;

public:
	CRF_StdSegNStateNode_WithoutDurLab_WithoutSegTransFtr(float* fb, QNUInt32 sizeof_fb, QNUInt32 lab, CRF_Model* crf, QNUInt32 nodeMaxDur, QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs);
	virtual ~CRF_StdSegNStateNode_WithoutDurLab_WithoutSegTransFtr();

	virtual double computeTransMatrix();
	virtual double computeAlpha();
	virtual double computeFirstAlpha();
	virtual double computeBeta(double scale=1.0);
////	virtual void setTailBeta();

	// [important note: this is different from other models which pass prev_lab instead of t_next_lab]
	virtual double computeExpF(double* ExpF, double* grad, double Zx, QNUInt32 t_next_lab);

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

	virtual double computeAlphaPlusTrans(double* nextNodeDenseTransMatrix, double* nextNodeDiagTransMatrix, double* nextNodeOffDiagTransMatrix);

	// These are the correct versions for CRF_StdSegNStateNode_WithoutDurLab_WithoutSegTransFtr.
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);    // without dur
//	virtual double getStateValue(QNUInt32 cur_lab, QNUInt32 cur_dur);
	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 cur_dur);

	// Disable all these functions by overriding them with exception handling. Use their modified versions above.
	virtual double getTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab, QNUInt32 cur_dur);   // with dur
//	virtual double getStateValue(QNUInt32 cur_lab);
//	virtual double getFullTransValue(QNUInt32 prev_lab, QNUInt32 cur_lab);

	// Disabled in CRF_StdSegNStateNode_WithoutDurLab_WithoutTransFtr
	virtual double getTempBeta(QNUInt32 next_lab, QNUInt32 next_dur);

	// only specific to CRF_StdSegNStateNode_WithoutDurLab_WithoutTransFtr
	virtual double getAlphaPlusTrans(QNUInt32 next_lab);
};

#endif /* CRF_STDSEGNSTATENODE_WITHOUTDURLAB_WITHOUTSEGTRANSFTR_H_ */
