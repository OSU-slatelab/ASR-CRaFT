#include "CRF_AISTrainer.h"

#define QN_UINT32_MAX (0xffffffff)

CRF_AISTrainer::CRF_AISTrainer(CRF_Model* crf_in, CRF_FeatureStreamManager* ftr_str_mgr, char* wt_fname)
	: CRF_Trainer(crf_in, ftr_str_mgr, wt_fname)
{
	l1alpha=0;
}

#ifdef OLDCODE
void CRF_AISTrainer::train()
{
	bool train_transitions=false;
	bool compute_transitions=true;
	ofstream ofile;
	QNUInt32 nStates = this->crf_ptr->getFeatureMap()->getNumStates();

	CRF_FeatureStream *crf_ftr_str=this->ftr_strm_mgr->trn_stream;
	CRF_LatticeBuilder *lb=new CRF_LatticeBuilder(crf_ftr_str,this->crf_ptr);
	CRF_NewLocalPosteriorBuilder *lpb=NULL;
	if (!train_transitions) {
		lpb=new CRF_NewLocalPosteriorBuilder(this->crf_ptr,true);
	}

	crf_ftr_str->rewind();

	QN_SegID segid = crf_ftr_str->nextseg();
	int count=0;

	int nframes;
	vector<double> modelgamma;
	vector<double> aligngamma;
	vector<double> modeltransgamma;
	vector<double> aligntransgamma;
	int lambdalen=this->crf_ptr->getLambdaLen();
	int nlabs=this->crf_ptr->getNLabs();
	double numeratorcounts[lambdalen];
	double denominatorcounts[lambdalen];
	for (int ii=0;ii<lambdalen;ii++) {
		numeratorcounts[ii]=0.0;
		denominatorcounts[ii]=0.0;
	}

		cout << "COMPUTING TRANSITIONS for " << nlabs << endl;

		size_t bunch_size=1024;
		float *ftr_buf=new float[crf_ftr_str->num_ftrs()*bunch_size];
		QNUInt32 *lab_buf=new QNUInt32[bunch_size];
		double *transitions=new double[nlabs*nlabs];
		double *counts=new double[nlabs];
		double epsilon=0.0000001;
		for (int i=0;i<nlabs;i++) {
			counts[i]=epsilon*nlabs;
			for (int j=0;j<nlabs;j++) {
				transitions[i*nlabs+j]=epsilon;
			}
		}


		QNUInt32 ftr_count;
		if (compute_transitions) {
			while (segid != QN_SEGID_BAD) {
			QNUInt32 last_lab;
			bool sent_start=true;

			do {
				ftr_count=crf_ftr_str->read(bunch_size,ftr_buf,lab_buf);
			//cout << "read " << ftr_count << " frames" << endl;

			//cout << lab_buf[0] << " " << endl;
				if (!sent_start && ftr_count>0) {
					transitions[last_lab*nlabs+lab_buf[0]]++;
					counts[last_lab]++;
				//cout << last_lab << "-" << lab_buf[0] << " " << endl;
				}

				last_lab=lab_buf[0];
				sent_start=false;
				for (QNUInt32 i=1; i<ftr_count; i++) {
					QNUInt32 lab=lab_buf[i];
				//cout << lab << " " << endl;
					transitions[last_lab*nlabs+lab]++;
					counts[last_lab]++;
				//cout << last_lab << "-" << lab << " " << endl;
					last_lab=lab;
				}
			} while(ftr_count>=bunch_size);

			segid = crf_ftr_str->nextseg();
		}
	}
	// reset stream
	crf_ftr_str->rewind();
	segid = crf_ftr_str->nextseg();

	cout << "initializing lambdas with transition biases" << endl;

	double *lambda=crf_ptr->getLambda();

	for (int i=0;i<nlabs;i++) {
		//cout << "ct(" << i << ")= " << counts[i];
		for (int j=0;j<nlabs;j++) {
			//cout << " " << j << ":" << transitions[i*nlabs+j];
			//if (counts[i]>0.0)
			if (compute_transitions) {
				transitions[i*nlabs+j]=log(transitions[i*nlabs+j]/counts[i]);
			} else {
				if (i==j) {
					transitions[i*nlabs+j]=log(0.5);
				} else {
					transitions[i*nlabs+j]=log(0.5/(nlabs-1));
				}
				//transitions[i*nlabs+j]=0.0;
			}
			QNUInt32 ix=this->crf_ptr->getFeatureMap()->getTransBiasIdx(i,j);
			// if QN_UINT32_MAX then this feature is not used
			// or the function is not defined by the feature map
			if (ix != QN_UINT32_MAX) {
				lambda[ix]=transitions[i*nlabs+j];
			}
			//cout << "/" << ix << "->" << transitions[i*nlabs+j];
		}
		//cout << endl;
	}

	//for (QNUInt32 ii=0;ii<lambdalen;ii++) {
	//	string s=this->crf_ptr->getFeatureMap()->getMapDescriptor(ii);
	//	cout << ii << ": " << s << " = " << lambda[ii] << endl;
	//}


	cout << "STARTING AIS TRAINING" << endl;
	cout << "l1alpha = " << this->l1alpha << endl;

	while (iCounter<this->maxIters) {
		for (int ii=0;ii<lambdalen;ii++) {
			numeratorcounts[ii]=0.0;
			denominatorcounts[ii]=0.0;
		}
		count=0;
		while (segid != QN_SEGID_BAD) {
			if (count % this->uttRpt == 0)
				cout << "Processing segment " << count << endl;
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
	//		cout << " ... gammas ... " << endl;
	//		cerr << "aligngamma.length=" << aligngamma.size() << endl;
	//		cerr << "modelgamma.length=" << modelgamma.size() << endl;

			CRF_StateVector *nodelist=lb->getNodeList();

			try {
				for (int j=0; j<nframes; j++) {
					for (int i=0; i<nlabs; i++) {
						//cerr << "frame " << j << " lab " << i;
						double alignalphabeta=expE(aligngamma.at(j*nlabs+i));
						double modelalphabeta=expE(modelgamma.at(j*nlabs+i));
												//cerr << " alignab " << alignalphabeta << " modelab " << modelalphabeta << endl;
						float *ftrs=nodelist->at(j)->getFtrBuffer();

#ifdef DEBUG
						if (alignalphabeta>0.001) {
							cerr << "frame " << j << " lab " << i << " gamma " << alignalphabeta << endl;
							for (int ii=0;ii<crf_ftr_str->num_ftrs();ii++) {
								cerr << " " << ftrs[ii];
							}
							cerr << endl;
							for (int ii=0;ii<crf_ftr_str->num_ftrs();ii++) {
								cerr << " " << numeratorcounts[this->crf_ptr->getFeatureMap()->getStateFeatureIdx(i,ii)];
							}
							cerr << endl;
						}
#endif
						if (alignalphabeta>0.0)
							this->crf_ptr->getFeatureMap()->computeStateExpF(ftrs,NULL,numeratorcounts,NULL,alignalphabeta,i,i,0);
						if (modelalphabeta>0.0)
							this->crf_ptr->getFeatureMap()->computeStateExpF(ftrs,NULL,denominatorcounts,NULL,modelalphabeta,i,i,0);
#ifdef DEBUG
						if (alignalphabeta>0.001) {
							for (int ii=0;ii<crf_ftr_str->num_ftrs();ii++) {
								cerr << " " << numeratorcounts[this->crf_ptr->getFeatureMap()->getStateFeatureIdx(i,ii)];
							}
							cerr << endl;
							cerr << endl;
						}
#endif

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
				segid=crf_ftr_str->nextseg();
				count++;


								//delete posteriorList;
				modelgamma.clear();
				modeltransgamma.clear();
				aligngamma.clear();
				aligntransgamma.clear();

			} catch (exception& e) {
				cerr << "Exception: " << e.what() << endl;
				exit(-1);
			}
		}
	#ifdef DEBUG
		cout << "numcounts: ";
		for (QNUInt32 ii=0;ii<lambdalen;ii++){
			cout << " " << ii << ":" << numeratorcounts[ii];
		}
		cout << endl;
		cout << "dencounts: ";
		for (QNUInt32 ii=0;ii<lambdalen;ii++){
			cout << " " << ii << ":" << denominatorcounts[ii];
		}
		cout << endl;
		//cout << " counts generated" << endl;
		cout << "constraint diffs > " << this->l1alpha << endl;
		for (QNUInt32 ii=0;ii<lambdalen;ii++) {
			if (abs(numeratorcounts[ii]-denominatorcounts[ii])>this->l1alpha) {
				string s=this->crf_ptr->getFeatureMap()->getMapDescriptor(ii);
				cout << ii << ": " << s << " = " << (numeratorcounts[ii]-denominatorcounts[ii]) << endl;
			}
		}
	#endif

		for (QNUInt32 ii=0;ii<lambdalen;ii++) {
			//string s=this->crf_ptr->getFeatureMap()->getMapDescriptor(ii);
			//cout << ii << ": " << s << " = " << lambda[ii] << " n=" << numeratorcounts[ii] << " d=" << denominatorcounts[ii] << endl;
			if (numeratorcounts[ii] > 0 && denominatorcounts[ii] > 0 &&
				((lambda[ii]!=0.0) ||
					(numeratorcounts[ii] > 0.1 && abs(numeratorcounts[ii]-denominatorcounts[ii])> this->l1alpha ))) {
				lambda[ii]=lambda[ii] + this->lr * log(numeratorcounts[ii]/denominatorcounts[ii]);
			}
			//cout << ii << ": " << s << " = " << lambda[ii] << " n=" << numeratorcounts[ii] << " d=" << denominatorcounts[ii] << endl;
		}

		string fname;
		stringstream ss;
		ss << this->weight_fname << ".i" << iCounter << ".out";
		ss >> fname;
		cout << "Writing Iteration " << iCounter << " weights to file " << fname << endl;
		bool chkwrite=this->crf_ptr->writeToFile(fname.c_str());
		if (!chkwrite) {
			cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
			exit(-1);
		}
		iCounter++;
		crf_ftr_str->rewind();
		segid=crf_ftr_str->nextseg();
	}
}

#else

void CRF_AISTrainer::train()
{
	bool train_transitions=false;
	bool compute_transitions=false;
	ofstream ofile;
	QNUInt32 nStates = this->crf_ptr->getFeatureMap()->getNumStates();

	CRF_FeatureStream *crf_ftr_str=this->ftr_strm_mgr->trn_stream;
	CRF_CountAccumulator *ca;
	int nthreads=this->ftr_strm_mgr->getNThreads();
	if (nthreads<=1) {
		ca=new CRF_CountAccumulator(this->crf_ptr,crf_ftr_str,this->uttRpt);
	} else {
		ca=new CRF_Pthread_CountAccumulator(this->crf_ptr,this->ftr_strm_mgr,nthreads,this->uttRpt);
	}
	int nframes;
	int lambdalen=this->crf_ptr->getLambdaLen();
	int nlabs=this->crf_ptr->getNLabs();

	double numeratorcounts[lambdalen];
	double denominatorcounts[lambdalen];


	cout << "COMPUTING TRANSITIONS for " << nlabs << endl;

	size_t bunch_size=1024;
	int nftrs=this->crf_ptr->getFeatureMap()->getNumStateFuncs(0);
	int nftrs_ftrstrm=this->ftr_strm_mgr->getNumFtrs();

	float *ftr_buf=new float[nftrs_ftrstrm*bunch_size];
	QNUInt32 *lab_buf=new QNUInt32[bunch_size];
	double *transitions=new double[nlabs*nlabs];
	double *counts=new double[nlabs];
	double *ftrcounts=new double[nlabs*nftrs];
	double *sumftrcounts=new double[nftrs];
	double *penaltyfunction= new double[lambdalen];
	double epsilon=0.0000001;
	for (int i=0;i<nlabs;i++) {
		counts[i]=epsilon*nlabs;
		for (int j=0;j<nlabs;j++) {
			transitions[i*nlabs+j]=epsilon;
		}
		for (int j=0;j<nftrs;j++) {
			ftrcounts[i*nftrs+j]=0.0;
		}
	}
	for (int i=0;i<nftrs;i++) {
		sumftrcounts[i]=0.0;
	}
	for (int i=0;i<lambdalen;i++) {
		penaltyfunction[i]=1.0;
	}

	QNUInt32 ftr_count;
	// this section of code looks at the label file and computes
	// transition probabilities to initialize the lambda vector
	if (compute_transitions) {
		crf_ftr_str->rewind();
		QN_SegID segid=crf_ftr_str->nextseg();

		while (segid != QN_SEGID_BAD) {
			QNUInt32 last_lab;
			bool sent_start=true;

			do {
				ftr_count=crf_ftr_str->read(bunch_size,ftr_buf,lab_buf);
				// catch overflow from >bunch_size sentences
				if (!sent_start && ftr_count>0) {
					transitions[last_lab*nlabs+lab_buf[0]]++;
					counts[last_lab]++;
				}

				last_lab=lab_buf[0];
				//crf_ptr->getFeatureMap()->accumulateFeatures(ftr_buf,ftrcounts,last_lab);
				//crf_ptr->getFeatureMap()->accumulateFeatures(ftr_buf,sumftrcounts,0);
				for (int ii=0;ii<nftrs_ftrstrm;ii+=2) {
					int ix=(int)ftr_buf[ii];
					//cout << ix << ":" << ftr_buf[ii+1] << " ";
					ftrcounts[last_lab*nftrs+ix]+=ftr_buf[ii+1];
					sumftrcounts[ix]+=ftr_buf[ii+1];
				}
				//cout << endl;
				sent_start=false;
				for (QNUInt32 i=1; i<ftr_count; i++) {
					QNUInt32 lab=lab_buf[i];
					transitions[last_lab*nlabs+lab]++;
					counts[last_lab]++;

					//for (int ii=0;ii<nftrs;ii++) {
					//	ftrcounts[lab*nftrs+ii]+=ftr_buf[i*nftrs+ii];
					//	sumftrcounts[ii]+=ftr_buf[i*nftrs+ii];
					//}
					for (int ii=0;ii<nftrs_ftrstrm;ii+=2) {
						int iii=ii+nftrs_ftrstrm*i;
						int ix=(int)ftr_buf[iii];
						//cout << ix << ":" << ftr_buf[iii+1] << " ";
						ftrcounts[lab*nftrs+ix]+=ftr_buf[iii+1];
						sumftrcounts[ix]+=ftr_buf[iii+1];
					}
					//cout << endl;
					//crf_ptr->getFeatureMap()->accumulateFeatures(ftr_buf,ftrcounts,lab);
					//crf_ptr->getFeatureMap()->accumulateFeatures(ftr_buf,sumftrcounts,0);

					last_lab=lab;
				}
			} while(ftr_count>=bunch_size);

			segid = crf_ftr_str->nextseg();
		}
	}
	// reset stream
	crf_ftr_str->rewind();

	cout << "initializing lambdas with transition biases" << endl;

	double *lambda=crf_ptr->getLambda();

	for (int i=0;i<nlabs;i++) {
		//cout << "ct(" << i << ")= " << counts[i];
		for (int j=0;j<nlabs;j++) {
			//cout << " " << j << ":" << transitions[i*nlabs+j];
			//if (counts[i]>0.0)
			if (compute_transitions) {
				transitions[i*nlabs+j]=log(transitions[i*nlabs+j]/counts[i]);
			} else {
				if (i==j) {
					transitions[i*nlabs+j]=log(0.5);
				} else {
					transitions[i*nlabs+j]=log(0.5/(nlabs-1));
				}
				//transitions[i*nlabs+j]=0.0;
			}
			QNUInt32 ix=this->crf_ptr->getFeatureMap()->getTransBiasIdx(i,j);
			// if QN_UINT32_MAX then this feature is not used
			// or the function is not defined by the feature map
			if (ix != QN_UINT32_MAX) {
				lambda[ix]=transitions[i*nlabs+j];
			}
			//cout << "/" << ix << "->" << transitions[i*nlabs+j];
		}
		//cout << endl;

		if (0) {
 		// now initialize ftr counts
		int maxix=0;
		double maxval=ftrcounts[i*nftrs];
		cout << "ftrct " << i << " 0:" << ftrcounts[i*nftrs];
		for (int j=1;j<nftrs;j++) {
			cout << " " << j << ":" << ftrcounts[i*nftrs+j];
			if (ftrcounts[i*nftrs+j]>maxval) {
				maxval=ftrcounts[i*nftrs+j];
				maxix=j;
			}
		}
		cout << endl;
		cout << endl;
		cout  << "Setting lab " << i << " ftr " << maxix << "to 1" << endl;
		QNUInt32 ix=this->crf_ptr->getFeatureMap()->getStateFeatureIdx(i,maxix);
		if (ix != QN_UINT32_MAX) {
			lambda[ix]=1.0;
		}
		}



	}

	if (1) {
 	cout << "sumftrct:";
	for(int ii=0;ii<nftrs;ii++) {
		cout << " " << ii << ":" << sumftrcounts[ii];
		for (int j=0;j<nlabs;j++) {
			penaltyfunction[this->crf_ptr->getFeatureMap()->getStateFeatureIdx(j,ii)]=log(sumftrcounts[ii]/nlabs);
		}
	}
	cout << endl;
	}
	cout << "STARTING AIS TRAINING" << endl;

	while (iCounter<this->maxIters) {
		cout << "initialize local vars" << endl;
		for (int ii=0;ii<lambdalen;ii++) {
			numeratorcounts[ii]=0.0;
			denominatorcounts[ii]=0.0;
		}

		// use the accumulator to get counts
		cout << "initialize CA vars" << endl;
		ca->reset();
		cout << "accumulate" << endl;
		ca->accumulate();
		cout << "results" << endl;
		ca->add_results(numeratorcounts,denominatorcounts);

		// now update the lambda vector
		for (QNUInt32 ii=0;ii<lambdalen;ii++) {
			cout << ii << " " << numeratorcounts[ii] << " " << denominatorcounts[ii] << endl;
			if (numeratorcounts[ii] > 0 && denominatorcounts[ii] > 0 &&
				((lambda[ii]!=0.0) ||
					(numeratorcounts[ii] > 0.1 && abs(numeratorcounts[ii]-denominatorcounts[ii])> this->l1alpha ))) {
				//lambda[ii]=lambda[ii] + this->lr * log(numeratorcounts[ii]/denominatorcounts[ii])/penaltyfunction[ii];
				double step=log(numeratorcounts[ii]/denominatorcounts[ii]);
				//step=this->lr *((step<-1.0)?-1.0:(step>1.0)?1.0:step);
				step= this->lr * step;
				lambda[ii]=lambda[ii] + step;

			}
		}

		// write out the weight file
		string fname;
		stringstream ss;
		ss << this->weight_fname << ".i" << iCounter << ".out";
		ss >> fname;
		cout << "Writing Iteration " << iCounter << " weights to file " << fname << endl;
		bool chkwrite=this->crf_ptr->writeToFile(fname.c_str());
		if (!chkwrite) {
			cerr << "ERROR! File " << fname << " unable to be opened for writing.  ABORT!" << endl;
			exit(-1);
		}
		iCounter++;
	}
}
#endif



