/*
 * CRF_CountAccumulator.h
 *
 *  Created on: Aug 24, 2009
 *      Author: fosler
 */

#ifndef CRF_COUNTACCUMULATOR_H_
#define CRF_COUNTACCUMULATOR_H_
#include "../../CRF.h"
#include "../../io/CRF_FeatureStream.h"
#include "../../nodes/CRF_StateVector.h"
#include "../../decoders/CRF_LatticeBuilder.h"
#include "../../ftrmaps/CRF_FeatureMap.h" // for QN_UINT32_MAX -- move elsewhere?

class CRF_CountAccumulator {
protected:
	bool count_transitions;
	CRF_Model *crf_ptr;
	CRF_FeatureStream *ftr_strm;
	CRF_LatticeBuilder *lb;
	int nlambdas;
	double *numeratorcounts;
	double *denominatorcounts;
	vector<double> modelgamma;
	vector<double> aligngamma;
	vector<double> modeltransgamma;
	vector<double> aligntransgamma;
	size_t bunch_size;
	float *ftr_buf;
	QNUInt32 *lab_buf;
	int uttRpt;
	int uttOffset;
	int threadno;
public:
	CRF_CountAccumulator(CRF_Model *crf_ptr_in,CRF_FeatureStream *ftr_strm_in, int uttRpt_in);
	virtual ~CRF_CountAccumulator();
	virtual void reset();
	virtual void accumulate();
	virtual void add_results(double *global_numcounts, double *global_dencounts);
};

#endif /* CRF_COUNTACCUMULATOR_H_ */
