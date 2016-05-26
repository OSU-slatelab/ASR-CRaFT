/*
 * CRF_Pthread_GradAccumulator.cpp
 *
 *  Created on: Jul 7, 2009
 *      Author: fosler
 */

#include "CRF_Pthread_GradAccumulator.h"


CRF_Pthread_GradAccumulator_Thread::CRF_Pthread_GradAccumulator_Thread(CRF_GradBuilder *bldr,double *gradient)
: builder(bldr), grad(gradient), tmpZx(0.), logLi(0.), tmpLogLi(0.),
  uttCount(0), segid(QN_SEGID_BAD), uttReport(100), minibatchPerThread(CRF_UINT32_MAX)
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

	// Added by Ryan
	this->tmpLogLi = 0.;

	this->uttCount=0;
	//this->segid=QN_SEGID_BAD;
	time_t rawtime;

	// Added by Ryan
	cout << " minibatch size per thread for thread " << this->threadId << ": " << this->minibatchPerThread << endl;

	// now loop over all utterances
	do {
		//cout << "thread " << this->threadId << " utt " << uCounter << " start" << endl;
		//cout << this->builder << endl;
		//cout << this->ftr_str << endl;
		//cout << this->grad << endl;
		time(&rawtime);
		char* ntime = ctime(&rawtime);
		if (this->uttReport>0 && (uttCount % this->uttReport == 0)) {
		//if (uttCount % 100 == 0) {
		cout << ntime << " Thread: " << this->threadId;
		cout << " Starting Utt: " << uttCount << endl;
		}

		// Modified by Ryan, to save tmpLogLi
		double uttTmpLogLi=this->builder->buildGradient(this->ftr_str,
														this->grad,
														&this->tmpZx);
		this->tmpLogLi += uttTmpLogLi;
		//cout << "thread " << this->threadId << " utt " << uCounter << " end" << endl;
		double uttLogLi = uttTmpLogLi - tmpZx;
		this->logLi += uttLogLi;

		if (this->uttReport>0 && (uttCount % this->uttReport == 0)) {
		//if (uttCount % 100 == 0) {
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time << " Thread: " << this->threadId;
			cout << " Finished Utt: " << uttCount << " logLi: " << logLi;
			cout << " Avg LogLi: " << logLi/(uttCount+1) << " Zx: " << tmpZx;
			cout << " Numerator: " << tmpLogLi << endl;
		}

		this->uttCount++;


		this->segid = ftr_str->nextseg();

		// Added by Ryan
		if (this->uttCount >= this->minibatchPerThread) {
			break;
		}

	} while (this->segid != QN_SEGID_BAD);



	// Commented out by Ryan, don't clean up stream due to minibatch
	// clean up stream
	//ftr_str->rewind();
	//this->segid = ftr_str->nextseg();

	cout << "exiting thread  " << this->threadId << endl;
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
		CRF_FeatureStream *ftr_str;

		if (nStreams == 1) {
			ftr_str = ftr_str_mgr->trn_stream;
		} else {
			ftr_str = ftr_str_mgr->getChild(stream)->trn_stream;
		}

		// Rewind the stream
		ftr_str->rewind();
		QN_SegID segid = ftr_str->nextseg();

		// SHOULD CODE BETTER TO THROW EXCEPTION
		if (segid == QN_SEGID_BAD) {
			cerr << "Feature stream contains no utterances!" << endl;
			exit(1);
		}

		//builders[stream]=CRF_GradBuilder::create(crf,useLogspace,nStates);
		//cerr << "Creating thread " << stream << endl;
		builders[stream]=CRF_GradBuilder::create(crf,this->objective);

		// Added by Ryan
		builders[stream]->setNodeList(new CRF_StateVector());

		threads[stream]=
			new CRF_Pthread_GradAccumulator_Thread(builders[stream],sgrad[stream]);

		// Added by Ryan
		threads[stream]->setSegID(segid);
		threads[stream]->setUttReport(this->uttReport);
		//threads[stream]->minibatchPerThread is QN_ALL by default, i.e. infinity

		//cerr << "Starting thread " << stream << endl;
		threads[stream]->start(ftr_str);
		//cerr << "thread " << stream << " started" << endl;
	}

	double totLogLi=0.;

	// now get them to join
	for (stream=0;stream<nStreams;stream++) {
		//cerr << "Joining thread " << stream << endl;
		threads[stream]->join();
		//cerr << "Thread joined " << endl;
		totLogLi+=threads[stream]->getLogLi();
		uCounter+=threads[stream]->getUttCount();

		for(QNUInt32 i=0;i<nlambda;i++)
			grad[i]+=sgrad[stream][i];

		//cerr << "Deleting builders " << stream << endl;
		delete builders[stream];
		//cerr << "Deleting thread " << stream << endl;
		delete threads[stream];
		//cerr << "Thread deleted" << endl;
	}
	//cerr << "Deleting array builders" << endl;
	delete[] builders;
	//cerr << "Deleting array threads" << endl;
	delete[] threads;
	//cerr << "Deleting array sgrad_data" << endl;
	delete[] sgrad_data;
	//cerr << "Deleting array sgrad" << endl;
	delete[] sgrad;

	*uttCount=uCounter;
	return totLogLi;
}

// Added by Ryan
/*
 *
 * CRF_Pthread_GradAccumulator::accumulateGradientMinibatch()
 *
 * The return value is tmpLogLi, i.e. only the numerator.
 * The denominator is return through the input argument *Zx_out.
 *
 * whereas CRF_Pthread_GradAccumulator::accumulateGradient() return LogLi directly,
 * which is = tmpLogLi - Zx_out.
 *
 */
double CRF_Pthread_GradAccumulator::accumulateGradientMinibatch(CRF_FeatureStreamManager* ftr_str_mgr,
												int nStreams,
												double* grad,
												double* Zx_out,
												QNUInt32 *uttCount,
												bool *endOfIter) {

	int stream;

	QNUInt32 uCounter = 0;
	QNUInt32 nlambda = this->crf->getLambdaLen();

	// set up stream gradients
	double **sgrad = new double*[nStreams];
	double *sgrad_data = sgrad[0] = new double[nStreams*nlambda];
	memset(sgrad[0],0,nStreams*nlambda*sizeof(double));

	for(QNUInt32 s = 1; s < nStreams; s++) {
		sgrad[s] = sgrad[s-1] + nlambda;
	}

	// initialize global gradient
	for (QNUInt32 i = 0; i < nlambda; i++) {
		grad[i] = 0.0;
	}

	// set up threads
	CRF_Pthread_GradAccumulator_Thread **threads = new CRF_Pthread_GradAccumulator_Thread *[nStreams];
	CRF_GradBuilder **builders = new CRF_GradBuilder *[nStreams];
	for (stream = 0; stream < nStreams; ++stream) {
		threads[stream] = NULL;
		builders[stream] = NULL;
	}

	// divide up the minibatch to nStreams.
	// if not divisible, the remainder is evenly spreaded to threads.
	QNUInt32 minibatch_per_thread, remainder;
	if (this->minibatch < nStreams) {
		cerr << "Error in CRF_Pthread_GradAccumulator::accumulateGradientMinibatch(): "
				<< "minibatch size (" << this->minibatch << ") is less than the number of threads (" << nStreams << ")." << endl;
		exit(1);
	}
	if (this->minibatch == CRF_UINT32_MAX) {
		minibatch_per_thread = CRF_UINT32_MAX;
		remainder = 0;
	} else {
		minibatch_per_thread = this->minibatch / nStreams;
		remainder = this->minibatch % nStreams;
	}

	int nEndOfStreams = 0;
	for (stream = 0; stream < nStreams; ++stream) {
		//cout << "Using stream " << stream << endl;
		CRF_FeatureStreamManager *child_ftr_str_mgr;

		if (nStreams == 1) {
			child_ftr_str_mgr = ftr_str_mgr;
		} else {
			child_ftr_str_mgr = ftr_str_mgr->getChild(stream);
		}

		if (!child_ftr_str_mgr) {
			cerr << "CRF_Pthread_GradAccumulator::accumulateGradientMinibatch() Error: Cannot get the child feature stream manager for "
					<< stream << endl;
			exit(1);
		}
		CRF_FeatureStream *ftr_str = child_ftr_str_mgr->trn_stream;

		// Rewind the stream
		//ftr_str->rewind();
		QN_SegID segid = ftr_str->nextseg();

		// SHOULD CODE BETTER TO THROW EXCEPTION
		if (segid == QN_SEGID_BAD) {
			++nEndOfStreams;
			//cerr << "Feature stream contains no utterances!" << endl;
			//exit(1);
		} else {
			//builders[stream]=CRF_GradBuilder::create(crf,useLogspace,nStates);
			//cerr << "Creating thread " << stream << endl;
			builders[stream] = CRF_GradBuilder::create(crf,this->objective);
			builders[stream]->setNodeList(new CRF_StateVector());
			threads[stream] =
				new CRF_Pthread_GradAccumulator_Thread(builders[stream],sgrad[stream]);
			threads[stream]->setSegID(segid);
			threads[stream]->setUttReport(this->uttReport);
			threads[stream]->setMinibatchPerThread(minibatch_per_thread + stream < remainder ? 1 : 0);

			//cerr << "Starting thread " << stream << endl;
			threads[stream]->start(ftr_str);
			//cerr << "thread " << stream << " started" << endl;
		}
	}

	double totTmpLogLi = 0.;
	*Zx_out = 0.;

	// now get them to join
	for (stream = 0; stream < nStreams; ++stream) {
		if (!threads[stream]) {
			continue;
		}
		//cerr << "Joining thread " << stream << endl;
		threads[stream]->join();
		//cerr << "Thread joined " << endl;

		totTmpLogLi += threads[stream]->getTmpLogLi();
		*Zx_out += threads[stream]->getTmpZx();
		uCounter += threads[stream]->getUttCount();

		if (threads[stream]->getSegID() == QN_SEGID_BAD) {
			++nEndOfStreams;
		}

		for(QNUInt32 i = 0; i < nlambda; ++i)
			grad[i] += sgrad[stream][i];

		//cerr << "Deleting builders " << stream << endl;
		delete builders[stream];
		//cerr << "Deleting thread " << stream << endl;
		delete threads[stream];
		//cerr << "Thread deleted" << endl;
	}

	if (nEndOfStreams == nStreams) {
		// if all streams are at the end, then all training segments have been read
		// so it's the end of iteration. rewind all the streams.
		*endOfIter = true;
		for (stream = 0; stream < nStreams; ++stream) {
			if (nStreams == 1) {
				ftr_str_mgr->trn_stream->rewind();
			} else {
				ftr_str_mgr->getChild(stream)->trn_stream->rewind();
			}
		}
	} else {
		*endOfIter = false;
	}

	//cerr << "Deleting array builders" << endl;
	delete[] builders;
	//cerr << "Deleting array threads" << endl;
	delete[] threads;
	//cerr << "Deleting array sgrad_data" << endl;
	delete[] sgrad_data;
	//cerr << "Deleting array sgrad" << endl;
	delete[] sgrad;

	*uttCount=uCounter;
	return totTmpLogLi;
}
