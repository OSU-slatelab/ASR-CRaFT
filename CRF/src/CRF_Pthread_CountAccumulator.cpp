/*
 * CRF_Pthread_CountAccumulator.cpp
 *
 *  Created on: Aug 24, 2009
 *      Author: fosler
 */

#include "CRF_Pthread_CountAccumulator.h"

CRF_Pthread_CountAccumulator_Thread::CRF_Pthread_CountAccumulator_Thread(CRF_Model *crf_ptr_in,CRF_FeatureStream *ftr_strm_in, int uttRpt_in)
: CRF_CountAccumulator(crf_ptr_in,ftr_strm_in,uttRpt_in) {

}

CRF_Pthread_CountAccumulator_Thread::~CRF_Pthread_CountAccumulator_Thread() {

}

int CRF_Pthread_CountAccumulator_Thread::join() {
	return pthread_join(this->threadId,NULL);
}

int CRF_Pthread_CountAccumulator_Thread::start() {
	int code = pthread_create(&threadId,
			NULL,
			CRF_Pthread_CountAccumulator_Thread::threadEntry,
			(void *)this);
	return code;
}

/*static */
void * CRF_Pthread_CountAccumulator_Thread::threadEntry(void * pthis)
{
	CRF_Pthread_CountAccumulator_Thread * pt = (CRF_Pthread_CountAccumulator_Thread*)pthis;
	cout << "Entering thread " << pt->threadId << endl;
	pt->reset();
	pt->accumulate();
	return pthis;
}

CRF_Pthread_CountAccumulator::CRF_Pthread_CountAccumulator(CRF_Model *crf_ptr_in,CRF_FeatureStreamManager *ftr_strm_in, int nthreads_in, int uttRpt_in)
:CRF_CountAccumulator(crf_ptr_in,ftr_strm_in->trn_stream,uttRpt_in) {

	this->ftr_str_mgr=ftr_strm_in;
	this->nthreads=nthreads_in;
	this->children = new CRF_Pthread_CountAccumulator_Thread *[nthreads];
	for (int i=0;i<this->nthreads;i++) {
		this->children[i]=new CRF_Pthread_CountAccumulator_Thread(this->crf_ptr,ftr_str_mgr->getChild(i)->trn_stream,this->uttRpt);
	}

}

CRF_Pthread_CountAccumulator::~CRF_Pthread_CountAccumulator() {
	for (int i=0;i<this->nthreads;i++) {
		delete this->children[i];
	}
	delete [] this->children;
}

void CRF_Pthread_CountAccumulator::reset() {
	for (int i=0;i<this->nthreads;i++) {
		this->children[i]->reset();
	}
	for (int i=0;i<this->nlambdas;i++) {
		this->numeratorcounts[i]=0.0;
		this->denominatorcounts[i]=0.0;
	}
}

void CRF_Pthread_CountAccumulator::accumulate() {

	for (int i=0;i<this->nthreads;i++) {
		this->children[i]->start();
	}
	for (int i=0;i<this->nthreads;i++) {
		this->children[i]->join();
		this->children[i]->add_results(this->numeratorcounts,this->denominatorcounts);
	}

	if (0) {
	cout << "pthreadnum:";
	for (int ii=0;ii<this->nlambdas;ii++) {
		cout << " " << this->numeratorcounts[ii];
	}
	cout << endl;

	cout << "pthreadden:";
	for (int ii=0;ii<this->nlambdas;ii++) {
		cout << " " << this->denominatorcounts[ii];
	}
	cout << endl;
	}
}



