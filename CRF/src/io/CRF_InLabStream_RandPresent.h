#ifndef CRF_INLABSTREAM_RANDPRESENT_H_
#define CRF_INLABSTREAM_RANDPRESENT_H_
/*
 * CRF_InLabStream_RandPresent.h
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 * Contains the class definitions for CRF_InLabStream_RandPresent
 * Uses the feature stream model/classes defined for the ICSI Quicknet neural networks
 * package for compatibility with ICSI Quicknet.
 */

#include <QN_streams.h>
#include <QN_seqgen.h>
#include "CRF_FeatureStream.h"

/*
 * class CRF_InLabStream_RandPresent
 *
 * Controls the presentation of input labels to the CRF training and
 * decoding code.  Features and labels are assumed to be held in two different files
 * according to this model.  The CRF_InLabStream_RandPresent object is based on the ICSI Quicknet
 * feature stream objects and implements similar functionality for compatibility between
 * the two packages.
 *
 * This feature stream allows access to segments in a random order.  Two possible randomization
 * methods (no replacement and with replacement) are supported.  To use segments in sequential
 * order, use the CRF_FeatureStream class instead.
 *
 */
class CRF_InLabStream_RandPresent : public QN_InLabStream
{
protected:
    QN_ClassLogger log;         // Logging object.
    QN_InLabStream& str;        // The stream we filtering.
	QN_SeqGen* seqgen;			// Randomizer
	QNUInt32 seed;				// Random number seed
	size_t cur_seg;
	size_t max_segs;
	size_t epoch;
	QNUInt32 real_seg;
public:
	CRF_InLabStream_RandPresent(int a_debug, const char* a_classname, const char* a_dbgname,
									QN_InLabStream& a_str, seqtype gentype, QNUInt32 seed=0);
	virtual ~CRF_InLabStream_RandPresent();
	size_t num_labs();
	QN_SegID nextseg();
	size_t read_labs(size_t cnt, QNUInt32* labs);
	size_t read_labs(size_t cnt, QNUInt32* labs, size_t stride);
	int rewind();
	size_t num_segs();
	size_t num_frames(size_t segno = QN_ALL);
	int get_pos(size_t* segno, size_t* frameno);
	QN_SegID set_pos(size_t segno, size_t frameno);
	seqtype seqgen_type;


};

#endif /*CRF_INLABSTREAM_RANDPRESENT_H_*/
