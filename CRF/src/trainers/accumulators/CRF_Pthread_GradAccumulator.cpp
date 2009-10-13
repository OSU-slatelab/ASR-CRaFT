/*
 * CRF_Pthread_GradAccumulator.cpp
 *
 *  Created on: Jul 7, 2009
 *      Author: fosler
 */

#include "CRF_Pthread_GradAccumulator.h"


CRF_Pthread_GradAccumulator_Thread::CRF_Pthread_GradAccumulator_Thread(CRF_GradBuilder *bldr,double *gradient)
: builder(bldr), grad(gradient), tmpZx(0.), logLi(0.)
{
}

CRF_Pthread_GradAccumulator_Thread::~CRF_Pthread_GradAccumulator_Thread() {

}

int CRF_Pthread_GradAccumulator_Thread::start(CRF_FeatureStream *ftr)
{
	//cout << "Starting new thread" << endl;
	this->ftr_str=ftr;
	int code = pthread_create(&threadId,
			                  NULL,
							  CRF_Pthread_GradAccumulator_Thread::threadEntry,
							  (void *)this);
	return code;
}

int CRF_Pthread_GradAccumulator_Thread::run()
{
	//cout << "running thread " << this->threadId << endl;
	this->logLi=0.;
	// make this local?
	this->tmpZx=0.;

	this->uttCount=0;
	QN_SegID segid=QN_SEGID_BAD;

	// now loop over all utterances
	do {
		//cout << "thread " << this->threadId << " utt " << uCounter << " start" << endl;
		//cout << this->builder << endl;
		//cout << this->ftr_str << endl;
		//cout << this->grad << endl;
		double tmpLogLi=this->builder->buildGradient(this->ftr_str,
														this->grad,
														&this->tmpZx);
		//cout << "thread " << this->threadId << " utt " << uCounter << " end" << endl;
		double uttLogLi=tmpLogLi - tmpZx;
		this->logLi += uttLogLi;
		this->uttCount++;


		segid = ftr_str->nextseg();
	} while (segid != QN_SEGID_BAD);



	// clean up stream
	ftr_str->rewind();
	ftr_str->nextseg();
	return 0;
}

/*static */
void * CRF_Pthread_GradAccumulator_Thread::threadEntry(void * pthis)
{
	CRF_Pthread_GradAccumulator_Thread * pt = (CRF_Pthread_GradAccumulator_Thread*)pthis;
	cout << "Entering thread " << pt->threadId << endl;
	pt->run( );
	return pthis;
}

int CRF_Pthread_GradAccumulator_Thread::join() {
	return pthread_join(this->threadId,NULL);
}

CRF_Pthread_GradAccumulator::CRF_Pthread_GradAccumulator(CRF_Model *myCrf,bool myLogspace,int myNStates)
 : CRF_GradAccumulator(myCrf,myLogspace,myNStates)
 {

}


CRF_Pthread_GradAccumulator::~CRF_Pthread_GradAccumulator() {

}

double CRF_Pthread_GradAccumulator::accumulateGradient(CRF_FeatureStreamManager* ftr_str_mgr,
												int nStreams,
												double* grad,
												QNUInt32 *uttCount) {

	CRF_Pthread_GradAccumulator_Thread **threads;
	int stream;

	time_t rawtime;
	double totalLogLi=0.0;
	QNUInt32 uCounter=0;
	double tmp_Zx;
	QNUInt32 nlambda=this->crf->getLambdaLen();

	// set up stream gradients
	double **sgrad=new double*[nStreams];
	double *sgrad_data=sgrad[0]=new double[nStreams*nlambda];
	memset(sgrad[0],0,nStreams*nlambda*sizeof(double));

	for(QNUInt32 s=1;s<nStreams;s++) {
		sgrad[s]=sgrad[s-1]+nlambda;
	}

	// initialize global gradient
	for (QNUInt32 i=0;i<nlambda;i++) {
		grad[i]=0.0;
	}

	// set up threads
	threads=new CRF_Pthread_GradAccumulator_Thread *[nStreams];
	CRF_GradBuilder **builders=new CRF_GradBuilder *[nStreams];

	for (stream=0;stream<nStreams;stream++) {
		//cout << "Using stream " << stream << endl;
		CRF_FeatureStream *ftr_str=ftr_str_mgr->getChild(stream)->trn_stream;

		// Rewind the stream
		ftr_str->rewind();
		QN_SegID segid = ftr_str->nextseg();

		// SHOULD CODE BETTER TO THROW EXCEPTION
		if (segid == QN_SEGID_BAD) {
			cerr << "Feature stream contains no utterances!" << endl;
			exit(1);
		}

		builders[stream]=CRF_GradBuilder::create(crf,useLogspace,nStates);
		threads[stream]=
			new CRF_Pthread_GradAccumulator_Thread(builders[stream],sgrad[stream]);

		threads[stream]->start(ftr_str);
	}

	double totLogLi=0.;

	// now get them to join
	for (stream=0;stream<nStreams;stream++) {
		threads[stream]->join();
		totLogLi+=threads[stream]->getLogLi();
		uCounter+=threads[stream]->getUttCount();

		for(QNUInt32 i=0;i<nlambda;i++)
			grad[i]+=sgrad[stream][i];

		delete builders[stream];
		delete threads[stream];
	}

	delete[] builders;
	delete[] threads;
	delete[] sgrad_data;
	delete[] sgrad;

	*uttCount=uCounter;
	return totLogLi;
}


