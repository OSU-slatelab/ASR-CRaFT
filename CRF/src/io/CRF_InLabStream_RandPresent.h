#ifndef CRF_INLABSTREAM_RANDPRESENT_H_
#define CRF_INLABSTREAM_RANDPRESENT_H_

#include <QN_streams.h>
#include <QN_seqgen.h>
#include "CRF_FeatureStream.h"

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
