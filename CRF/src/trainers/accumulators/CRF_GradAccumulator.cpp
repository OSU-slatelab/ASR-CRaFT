/*
 * CRF_GradAccumulator.cpp
 *
 *  Created on: Jul 7, 2009
 *      Author: fosler
 */

#include "CRF_GradAccumulator.h"
#include "../gradbuilders/CRF_GradBuilder.h"

CRF_GradAccumulator::CRF_GradAccumulator(CRF_Model *myCrf,
										bool myLogspace,
										int myNStates)
	: crf(myCrf), useLogspace(myLogspace), nStates(myNStates) {

	uttReport=0;
}

CRF_GradAccumulator::~CRF_GradAccumulator() {

}

double CRF_GradAccumulator::accumulateGradient(CRF_FeatureStreamManager* ftr_str_mgr,
												int nStreams,
												double* grad,
												QNUInt32 *uttCount) {

	CRF_GradBuilder *gbuild=CRF_GradBuilder::create(crf,useLogspace,nStates);

	if (nStreams!=1) {
		cerr << "CRF_GradAccumulator can only handle one stream" << endl;
		exit(1);
	}

	CRF_FeatureStream *ftr_str=ftr_str_mgr->trn_stream;

	// Rewind the stream
	ftr_str->rewind();
	QN_SegID segid = ftr_str->nextseg();

	// SHOULD CODE BETTER TO THROW EXCEPTION
	if (segid == QN_SEGID_BAD) {
		cerr << "Feature stream contains no utterances!" << endl;
		exit(1);
	}

	time_t rawtime;
	double totLogLi=0.0;
	unsigned int uCounter=0;
	double tmp_Zx;
	QNUInt32 nlambda=this->crf->getLambdaLen();

	for (QNUInt32 i=0;i<nlambda;i++) {
		grad[i]=0.0;
	}

	// now loop over all utterances
	do {
		// each call to buildGradient updates grad
		double tmpLogLi=gbuild->buildGradient(ftr_str,grad,&tmp_Zx);
		double logLi=tmpLogLi - tmp_Zx;
		totLogLi += logLi;

		if (this->uttReport>0 && (uCounter % this->uttReport == 0)) {
			time(&rawtime);
			char* time = ctime(&rawtime);
			time[strlen(time)-1]='\0';
			cout << time << " Finished Utt: " << uCounter << " logLi: " << logLi;
			cout << " Avg LogLi: " << totLogLi/(uCounter+1) << " Zx: " << tmp_Zx;
			cout << " Numerator: " << tmpLogLi << endl;
		}
		uCounter++;

		segid = ftr_str->nextseg();
	} while (segid != QN_SEGID_BAD);


	delete gbuild;
	    // clean up stream
	ftr_str->rewind();
	ftr_str->nextseg();
	*uttCount=uCounter;
	return totLogLi;
}
