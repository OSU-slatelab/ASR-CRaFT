/*
 * CRF_StateVector.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */

#include "CRF_StateVector.h"

#include "CRF_StateNode.h"

/*
 * CRF_StateVector constructor
 */
CRF_StateVector::CRF_StateVector()
{
}

/*
 * CRF_StateVector destructor
 */
CRF_StateVector::~CRF_StateVector()
{
	this->deleteAll();
	//cerr << "in statevector destructor" << endl;
	// vector destroys objects, but not pointers to objects
}

/*
 * CRF_StateVector::set
 *
 * Input: idx - index for current state node being created/set
 *        new_buf, num_ftrs, lab_buf, crf_in - see constructor for CRF_StateNode
 */

void CRF_StateVector::set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in)
{
	if (idx >= this->size() ) {
		this->push_back( CRF_StateNode::createStateNode(new_buf,num_ftrs,lab_buf,crf_in));
	}
	else {
		this->at(idx)->reset(new_buf,num_ftrs,lab_buf,crf_in);
	}
}

// Added by Ryan
/*
 * CRF_StateVector::set
 *
 * Input: idx - index for current state node being created/set
 *        new_buf, num_ftrs, lab_buf, crf_in,
 *        nodeMaxDur, prevNode_nLabs, nextNode_nActualLabs
 *        - see constructor for CRF_StateNode
 */
void CRF_StateVector::set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs,
		QNUInt32 lab_buf, CRF_Model* crf_in, QNUInt32 nodeMaxDur,
		QNUInt32 prevNode_nLabs, QNUInt32 nextNode_nActualLabs)
{
	if (idx >= this->size() ) {
		this->push_back( CRF_StateNode::createStateNode(new_buf,num_ftrs,lab_buf,crf_in,nodeMaxDur,prevNode_nLabs,nextNode_nActualLabs));
	}
	else {
		// By Ryan: For CRF_StdSegStateNode, the parameters nodeMaxDur, prevNode_nLabs and nextNode_nActualLabs should not be changed for the same node index.
		this->at(idx)->reset(new_buf,num_ftrs,lab_buf,crf_in);
	}
}

/*
 * CRF_StateVector::setNodeCount
 *
 * Input: cnt - value to set the nodeCount to
 *
 * nodeCount is different from size because the StateVector can be reused.  nodecount should be
 * set to the length of the current sequence.
 */
void CRF_StateVector::setNodeCount(QNUInt32 cnt)
{
	this->nodeCount=cnt;
}

/*
 * CRF_StateVector::getNodeCount
 *
 * returns: index of final node in this sequence
 */
QNUInt32 CRF_StateVector::getNodeCount()
{
	return this->nodeCount;
}

/*
 * CRF_StateVector::deleteAll
 *
 * deletes the objects in the statevector
 */
void CRF_StateVector::deleteAll()
{
	//cerr << "in statevector::deleteall" << endl;
	for (size_t n=0;n<this->size();n++) {
		delete this->at(n);
	}
	//cerr << "leaving statevector::deleteall" << endl;
}
