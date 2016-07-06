/*
 * CRF_Minibatch_GradAccumulator.cpp
 *
 *  Created on: Nov 22, 2014
 *      Author: Yanzhang (Ryan) He
 */

#include "CRF_Minibatch_GradAccumulator.h"

CRF_Minibatch_GradAccumulator_Thread::CRF_Minibatch_GradAccumulator_Thread(CRF_GradBuilder *bldr, double *gradient)
: gBuilder(bldr), grad(gradient), bunchLogLiDenom(0.), bunchLogLiNumer(0.), bunchLogLi(0.),
  uttCount(0), segid(QN_SEGID_BAD), uttReport(100), minibatchPerThread(CRF_UINT32_MAX)
{
}

CRF_Minibatch_GradAccumulator_Thread::~CRF_Minibatch_GradAccumulator_Thread() {

}

int CRF_Minibatch_GradAccumulator_Thread::start(CRF_FeatureStream *ftr)
{
	//cout << "Starting new thread" << endl;
	this->ftr_str=ftr;
	int code = pthread_create(&threadId,
			                  NULL,
							  CRF_Minibatch_GradAccumulator_Thread::threadEntry,
							  (void *)this);
	return code;
}

int CRF_Minibatch_GradAccumulator_Thread::run()
{
	//cout << "running thread " << this->threadId << endl;
	this->bunchLogLi = 0.;
	this->bunchLogLiDenom = 0.;
	this->bunchLogLiNumer = 0.;
	this->uttCount=0;
	time_t rawtime;

	double uttLogLiDenom = 0.;

	// Added by Ryan
	//cout << " minibatch size per thread for thread " << this->threadId << ": " << this->minibatchPerThread << endl;

	// now loop over all utterances
	do {
		//cout << "thread " << this->threadId << " utt " << uCounter << " start" << endl;
		//cout << this->builder << endl;
		//cout << this->ftr_str << endl;
		//cout << this->grad << endl;
		time(&rawtime);
		char* ntime = ctime(&rawtime);
		if (this->uttReport>0 && this->uttCount > 0 && (this->uttCount % this->uttReport == 0)) {
			//if (uttCount % 100 == 0) {
			cout << ntime << " Thread: " << this->threadId;
			cout << " Starting Utt: " << uttCount << endl;
		}

		// Modified by Ryan, to save tmpLogLi
		double uttLogLiNumer = this->gBuilder->buildGradient(this->ftr_str,
														this->grad,
														&uttLogLiDenom);
		this->bunchLogLiNumer += uttLogLiNumer;
		this->bunchLogLiDenom += uttLogLiDenom;
		//cout << "thread " << this->threadId << " utt " << uCounter << " end" << endl;
		double uttLogLi = uttLogLiNumer - uttLogLiDenom;
		this->bunchLogLi += uttLogLi;

		if (this->uttReport>0 && this->uttCount > 0 && (this->uttCount % this->uttReport == 0)) {
		//if (uttCount % 100 == 0) {
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time << " Thread: " << this->threadId;
			cout << " Finished Utt: " << uttCount;
			cout << " Utt Numerator: " << uttLogLiNumer;
			cout << " Utt Zx: " << uttLogLiDenom;
			cout << " Utt LogLi: " << uttLogLi;
			cout << " Thread-Avg LogLi: " << this->bunchLogLi / (uttCount+1);
			cout << endl;
		}

		this->uttCount++;

		this->segid = ftr_str->nextseg();

		if (this->uttCount >= this->minibatchPerThread) {
			break;
		}

	} while (this->segid != QN_SEGID_BAD);

	// don't rewind here since next minibatch will continue from this segment

	//cout << "exiting thread  " << this->threadId << endl;
	return 0;
}

/*static */
void * CRF_Minibatch_GradAccumulator_Thread::threadEntry(void * pthis)
{
	CRF_Minibatch_GradAccumulator_Thread * pt = (CRF_Minibatch_GradAccumulator_Thread*)pthis;
	//cout << "Entering thread " << pt->threadId << endl;
	pt->run( );
	return pthis;
}

int CRF_Minibatch_GradAccumulator_Thread::join() {
	return pthread_join(this->threadId,NULL);
}



CRF_Minibatch_GradAccumulator::CRF_Minibatch_GradAccumulator(
		CRF_Model *myCrf, CRF_FeatureStreamManager* myFtrStrmMgr, QNUInt32 myNStreams)
	: crf(myCrf), nStreams(myNStreams) {

	this->uttReport = 0;
	this->objective = EXPF;
	this->minibatch = CRF_UINT32_MAX;

	this->ftrStrms = new CRF_FeatureStream*[nStreams];
	if (!myFtrStrmMgr) {
		cerr << "CRF_Minibatch_GradAccumulator::CRF_Minibatch_GradAccumulator() Error: "
				<< "Cannot open the feature stream manager myFtrStrmMgr" << endl;
		exit(1);
	}
	assert(myNStreams == myFtrStrmMgr->getNThreads());

	if (nStreams == 1) {
		this->ftrStrms[0] = myFtrStrmMgr->trn_stream;
	} else {
		for (int stream = 0; stream < nStreams; ++stream) {
			CRF_FeatureStreamManager *child_ftr_str_mgr;
			child_ftr_str_mgr = myFtrStrmMgr->getChild(stream);
			if (!child_ftr_str_mgr) {
				cerr << "CRF_Minibatch_GradAccumulator::CRF_Minibatch_GradAccumulator() Error: "
						<< "Cannot get the child feature stream manager for " << stream << endl;
				exit(1);
			}
			this->ftrStrms[stream] = child_ftr_str_mgr->trn_stream;
		}
	}

	this->strmsSegids = new QN_SegID[nStreams];
	this->gBuilders = new CRF_GradBuilder *[nStreams];
	for (int stream = 0; stream < nStreams; ++stream) {
		this->strmsSegids[stream] = QN_SEGID_BAD;
		this->gBuilders[stream] = CRF_GradBuilder::create(crf, this->objective);
		this->gBuilders[stream]->setNodeList(new CRF_StateVector());
	}
}

CRF_Minibatch_GradAccumulator::~CRF_Minibatch_GradAccumulator() {

	delete[] this->ftrStrms;
	delete[] this->strmsSegids;

	for (int stream = 0; stream < nStreams; ++stream) {
		delete this->gBuilders[stream];
	}
	delete[] this->gBuilders;
}

void CRF_Minibatch_GradAccumulator::setMinibatch(QNUInt32 minibatch) {

	if (minibatch < this->nStreams) {
		cerr << "CRF_Minibatch_GradAccumulator::setMinibatch() Error: "
				<< "minibatch size (" << minibatch << ") is less than "
				<< "the number of threads (" << this->nStreams << ")." << endl;
		exit(1);
	}
	if (minibatch == 0) {
		// totally batch
		this->minibatch = CRF_UINT32_MAX;
	} else {
		this->minibatch = minibatch;
	}
}

void CRF_Minibatch_GradAccumulator::rewindAllAndNextSegs() {

	for (int stream = 0; stream < this->nStreams; ++stream) {
		this->ftrStrms[stream]->rewind();
		this->strmsSegids[stream] = this->ftrStrms[stream]->nextseg();
	}
}

// Added by Ryan
/*
 *
 * CRF_Minibatch_GradAccumulator::accumulateGradient()
 *
 * The returned value is the numerator of the log-likelihood of the minibatch.
 * The denominator is returned through the input argument *Zx_out.
 *
 * whereas CRF_Pthread_GradAccumulator::accumulateGradient() returns LogLi directly,
 * which is = numerator - denominator.
 *
 */
double CRF_Minibatch_GradAccumulator::accumulateGradient(
		double* grad, double* Zx_out, QNUInt32 *uttCount, bool *isEndOfIter) {

	*uttCount = 0;
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
	CRF_Minibatch_GradAccumulator_Thread **threads = new CRF_Minibatch_GradAccumulator_Thread *[nStreams];
	for (int stream = 0; stream < nStreams; ++stream) {
		threads[stream] = NULL;
	}

	// divide up the minibatch to nStreams.
	// if not divisible, the remainder is evenly spread to threads.
	QNUInt32 minibatch_per_thread, remainder;
	if (minibatch < nStreams) {
		cerr << "CRF_Minibatch_GradAccumulator::accumulateGradient() Error: "
				<< "minibatch size (" << minibatch << ") is less than the number of threads (" << nStreams << ")." << endl;
		exit(1);
	}
	if (minibatch == CRF_UINT32_MAX) {
		minibatch_per_thread = CRF_UINT32_MAX;
		remainder = 0;
	} else {
		minibatch_per_thread = minibatch / nStreams;
		remainder = minibatch % nStreams;
	}

	int nStreamsEnd = 0;
	for (int stream = 0; stream < nStreams; ++stream) {
		//cout << "Using stream " << stream << endl;
		CRF_FeatureStream *ftr_str = this->ftrStrms[stream];

		if (this->strmsSegids[stream] == QN_SEGID_BAD) {
			++nStreamsEnd;
		} else {
			//builders[stream]=CRF_GradBuilder::create(crf,useLogspace,nStates);
			//cerr << "Creating thread " << stream << endl;
			threads[stream] =
				new CRF_Minibatch_GradAccumulator_Thread(this->gBuilders[stream],sgrad[stream]);
			threads[stream]->setSegID(this->strmsSegids[stream]);
			threads[stream]->setUttReport(this->uttReport);
			threads[stream]->setMinibatchPerThread(minibatch_per_thread + (stream < remainder ? 1 : 0));

			//cerr << "Starting thread " << stream << endl;
			threads[stream]->start(ftr_str);
			//cerr << "thread " << stream << " started" << endl;
		}
	}

	if (nStreamsEnd == nStreams) {
		cerr << "All feature streams are at the end! "
				<< "You don't have any utterances or you forget to rewind all the streams." << endl
				<< "For the latter case, run rewindAllAndNextSegs()." << endl;
		exit(1);
	}

	double totLogLiNumer = 0.;
	*Zx_out = 0.;

	// now get them to join
	int nStreams_active = 0;
	for (int stream = 0; stream < nStreams; ++stream) {
		if (!threads[stream]) {
			continue;
		}
		//cerr << "Joining thread " << stream << endl;
		threads[stream]->join();
		//cerr << "Thread joined " << endl;

		++nStreams_active;

		totLogLiNumer += threads[stream]->getBunchLogLiNumer();
		*Zx_out += threads[stream]->getBunchLogLiDenom();
		*uttCount += threads[stream]->getUttCount();

		this->strmsSegids[stream] = threads[stream]->getSegID();
		if (this->strmsSegids[stream] == QN_SEGID_BAD) {
			++nStreamsEnd;
		}

		for(QNUInt32 i = 0; i < nlambda; ++i) {
			grad[i] += sgrad[stream][i];
		}

		//cerr << "Deleting thread " << stream << endl;
		delete threads[stream];
		//cerr << "Thread deleted" << endl;
	}

	// averaging the gradient over streams
	for(QNUInt32 i = 0; i < nlambda; ++i) {
		grad[i] /= nStreams_active;
	}

	// if all streams are at the end, then all training segments have been read
	// i.e. it's the end of iteration.
	*isEndOfIter = (nStreamsEnd == nStreams);

	//cerr << "Deleting array threads" << endl;
	delete[] threads;
	//cerr << "Deleting array sgrad_data" << endl;
	delete[] sgrad_data;
	//cerr << "Deleting array sgrad" << endl;
	delete[] sgrad;

	return totLogLiNumer;
}
