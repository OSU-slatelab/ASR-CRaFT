#ifndef CRF_FEATURESTREAM_H_
#define CRF_FEATURESTREAM_H_
/*
 * CRF_FeatureStream.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_FeatureStream
 * Uses the feature stream model/classes defined for the ICSI Quicknet neural networks
 * package for compatibility with ICSI Quicknet.
 */

// Added by Ryan
#ifndef CRF_LAB_BAD
#define CRF_LAB_BAD (0xffffffff)
#endif

#include "QuickNet.h"
#include "QN_streams.h"
#include "QN_seqgen.h"

#include "../CRF.h"

/*
 * class CRF_FeatureStream
 *
 * Controls the presentation of input features and labels to the CRF training and
 * decoding code.  Features and labels are assumed to be held in two different files
 * according to this model.  The CRF_FeatureStream object is based on the ICSI Quicknet
 * file stream objects and implements similar functionality for compatibility between
 * the two packages.
 *
 */
class CRF_FeatureStream
{
	friend class CRF_FeatureStreamManager;
protected:
	QN_InFtrStream* ftr_stream;
	QN_InLabStream* lab_stream;
	int debug;
	QNUInt32 start_offset;
	QNUInt32 numsegs;
	//unsigned long numsegs;
	QN_SegID segid;

public:
	// Modified by Ryan
	// allowing for initialization of start_offset, numsegs, which is needed by the multi-threaded version
	CRF_FeatureStream(QN_InFtrStream* ftr, QN_InLabStream* lab, int dbg=0, QNUInt32 offset=0, QNUInt32 segs=QN_ALL);
	CRF_FeatureStream(QN_InFtrStream* ftr, int dbg=0, QNUInt32 offset=0, QNUInt32 segs=QN_ALL);

	virtual ~CRF_FeatureStream();
	virtual QN_SegID nextseg();
	virtual size_t read(size_t bs, float* fb, QNUInt32* lb);
	virtual int get_pos(size_t* segno, size_t* frameno);
	virtual CRF_FeatureStream* join(CRF_FeatureStream* in_stream);
	virtual void rewind();
	virtual QNUInt32 num_ftrs();
	virtual void set_pos(QNUInt32 segno, QNUInt32 frmno);
	virtual QNUInt32 num_segs();
	virtual void display();
	virtual void view(QNUInt32 startseg,QNUInt32 nsegs);

	// Added by Ryan
	virtual QNUInt32 num_labs();
};

#endif /*CRF_FEATURESTREAM2_H_*/
