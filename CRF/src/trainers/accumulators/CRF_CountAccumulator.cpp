/*
 * CRF_CountAccumulator.cpp
 *
 *  Created on: Aug 24, 2009
 *      Author: fosler
 */

#include "CRF_CountAccumulator.h"

CRF_CountAccumulator::CRF_CountAccumulator(CRF_Model *crf_ptr_in,
										   CRF_FeatureStream *ftr_strm_in,
										   int uttRpt_in) :
	crf_ptr(crf_ptr_in), ftr_strm(ftr_strm_in), uttRpt(uttRpt_in) {
	this->nlambdas=this->crf_ptr->getLambdaLen();
	this->lb=new CRF_LatticeBuilder(this->ftr_strm,this->crf_ptr);
	this->numeratorcounts=new double[this->nlambdas];
	this->denominatorcounts=new double[this->nlambdas];

	this->bunch_size=1024;
	this->ftr_buf=new float[ftr_strm->num_ftrs()*bunch_size];
	this->lab_buf=new QNUInt32[bunch_size];
	this->uttOffset=0;
	this->threadno=-1;
}

CRF_CountAccumulator::~CRF_CountAccumulator() {
	delete this->lb;
	delete [] this->ftr_buf;
	delete [] this->lab_buf;
	delete [] this->numeratorcounts;
	delete [] this->denominatorcounts;
}

void CRF_CountAccumulator::reset() {
	for (int ii=0;ii<this->nlambdas;ii++) {
		numeratorcounts[ii]=0.0;
		denominatorcounts[ii]=0.0;
	}
}

void CRF_CountAccumulator::accumulate() {
	int count=0;
	bool train_transitions=false;
	int nframes;
	int nlabs=this->crf_ptr->getNLabs();

	this->ftr_strm->rewind();
	QN_SegID segid=this->ftr_strm->nextseg();

	while (segid != QN_SEGID_BAD) {
		if (this->uttRpt > 0 && ((count+this->uttOffset) % this->uttRpt == 0)) {
			cout << "Processing segment " << (count+this->uttOffset);
			if (this->threadno>=0) {
				cout << " (thread " << this->threadno << ")";
			}
			cout << endl;
		}

		try {
			if (train_transitions)
				nframes=lb->getAlignmentGammas(&modelgamma,&aligngamma,&modeltransgamma,&aligntransgamma);
			else {
				nframes=lb->getAlignmentGammas(&modelgamma,&aligngamma,NULL,NULL);
			}
		}
		catch (exception &e) {
			cerr << "Exception: " << e.what() << endl;
			exit(-1);
		}


		CRF_StateVector *nodelist=lb->getNodeList();

		try {
			for (int j=0; j<nframes; j++) {
				float *ftrs=nodelist->at(j)->getFtrBuffer();
				QNUInt32 lab=nodelist->at(j)->getLabel();

				for (int i=0; i<nlabs; i++) {

					double alignalphabeta=expE(aligngamma.at(j*nlabs+i));
					double modelalphabeta=expE(modelgamma.at(j*nlabs+i));

					//cout << "frame " << j << " lab " << i << " align: " << alignalphabeta << " model: " << modelalphabeta << endl;
					alignalphabeta=(lab==i)?1.0:0.0;

					if (alignalphabeta>0.0)
						this->crf_ptr->getFeatureMap()->computeStateExpF(ftrs,NULL,numeratorcounts,NULL,alignalphabeta,i,i,0);
					if (modelalphabeta>0.0)
						this->crf_ptr->getFeatureMap()->computeStateExpF(ftrs,NULL,denominatorcounts,NULL,modelalphabeta,i,i,0);

					if (train_transitions && j>0) {
						for (int k=0; k<nlabs; k++) {
							alignalphabeta=expE(aligntransgamma.at((j-1)*nlabs*nlabs+i*nlabs+k));
							modelalphabeta=expE(modeltransgamma.at((j-1)*nlabs*nlabs+i*nlabs+k));
							if (alignalphabeta>0.0)
								this->crf_ptr->getFeatureMap()->computeTransExpF(ftrs,NULL,numeratorcounts,NULL,alignalphabeta,i,k,i,k,0);
							if (modelalphabeta>0.0)
								this->crf_ptr->getFeatureMap()->computeTransExpF(ftrs,NULL,denominatorcounts,NULL,modelalphabeta,i,k,i,k,0);

						}
					}
				}

			}
			segid=this->ftr_strm->nextseg();
			count++;

			modelgamma.clear();
			modeltransgamma.clear();
			aligngamma.clear();
			aligntransgamma.clear();

		} catch (exception& e) {
			cerr << "Exception: " << e.what() << endl;
			exit(-1);
		}
	}
}

void CRF_CountAccumulator::add_results(double *num_in, double *den_in) {
	for (int ii=0;ii<this->nlambdas;ii++) {
		num_in[ii]+=numeratorcounts[ii];
		den_in[ii]+=denominatorcounts[ii];
	}
}
