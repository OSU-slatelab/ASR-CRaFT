#include "CRF_StateVector.h"

#include "CRF_StateNode.h"

CRF_StateVector::CRF_StateVector()
{
}

CRF_StateVector::~CRF_StateVector()
{
	this->deleteAll();
	//cerr << "in statevector destructor" << endl;
	// vector destroys objects, but not pointers to objects
}

void CRF_StateVector::set(QNUInt32 idx, float* new_buf, QNUInt32 num_ftrs, QNUInt32 lab_buf, CRF_Model* crf_in)
{
	if (idx >= this->size() ) {
		this->push_back( CRF_StateNode::createStateNode(new_buf,num_ftrs,lab_buf,crf_in));
	}
	else {
		this->at(idx)->reset(new_buf,num_ftrs,lab_buf,crf_in);
	}
}

void CRF_StateVector::setNodeCount(QNUInt32 cnt)
{
	this->nodeCount=cnt;
}

QNUInt32 CRF_StateVector::getNodeCount()
{
	return this->nodeCount;
}

void CRF_StateVector::deleteAll()
{
	//cerr << "in statevector::deleteall" << endl;
	for (size_t n=0;n<this->size();n++) {
		delete this->at(n);
	}
	//cerr << "leaving statevector::deleteall" << endl;
}
