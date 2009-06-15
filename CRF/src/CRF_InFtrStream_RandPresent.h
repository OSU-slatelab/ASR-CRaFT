#ifndef CRF_INFTRSTREAM_RANDPRESENT_H_
#define CRF_INFTRSTREAM_RANDPRESENT_H_

#include "CRF.h"
#include "CRF_FeatureStream.h"

class CRF_InFtrStream_RandPresent : public QN_InFtrStream
{
protected:
    QN_ClassLogger log;         // Logging object.
    QN_InFtrStream& str;        // The stream we filtering.
    seqtype seqgen_type;
   	QN_SeqGen* seqgen;			// Randomizer
	QNUInt32 seed;				// Random number seed
	size_t cur_seg;
	QNUInt32 real_seg;
	size_t max_segs;
	size_t epoch;

public:
	CRF_InFtrStream_RandPresent(int a_debug, const char* a_classname, const char* a_dbgname,
									QN_InFtrStream& a_str, seqtype gen_type, QNUInt32 seed=0);
	virtual ~CRF_InFtrStream_RandPresent();	
    size_t num_ftrs();
    QN_SegID nextseg();
    size_t read_ftrs(size_t cnt, float* ftrs);
    int rewind();
    size_t num_segs();
    size_t num_frames(size_t segno = QN_ALL);
    int get_pos(size_t* segno, size_t* frameno);
    QN_SegID set_pos(size_t segno, size_t frameno);



};

#endif /*CRF_INFTRSTREAM_RANDPRESENT_H_*/
