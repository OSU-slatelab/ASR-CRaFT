#ifndef CRF_STATEVECTOR_H_
#define CRF_STATEVECTOR_H_
/*
 * CRF_StateVector.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_StateVector
 */
#include <vector>

#include "../CRF.h"
#include "../CRF_Model.h"
#include "CRF_StateNode.h"

/*
 * class CRF_StateVector
 *
 * Used in training and decoding processing.  Holds the nodes across the entire linear chain crf.
 *
 * CRF_StateVector uses the factor class createStateNode in CRF_StateNode to generate the right
 * kind of state based on the crf model used in the function "set".
 *
 */
class CRF_StateVector : public vector <CRF_StateNode*>
{
private:
	QNUInt32 nodeCount;
public:
	CRF_StateVector();
	virtual ~CRF_StateVector();
	virtual void set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in);
	virtual void setNodeCount(QNUInt32 cnt);
	virtual QNUInt32 getNodeCount();
	virtual void deleteAll();

	// Added by Ryan
	virtual void set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs,
			QNUInt32 lab_buf, CRF_Model* crf_in, QNUInt32 nodeMaxDur,
			QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs);
};

#endif /*CRF_STATEVECTOR_H_*/
