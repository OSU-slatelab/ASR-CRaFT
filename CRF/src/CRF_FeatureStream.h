#ifndef CRF_FEATURESTREAM_H_
#define CRF_FEATURESTREAM_H_

#include "QuickNet.h"
#include "QN_streams.h"
#include "QN_seqgen.h"

#include "CRF.h"

class CRF_FeatureStream
{
	friend class CRF_FeatureStreamManager;
protected:
	QN_InFtrStream* ftr_stream;
	QN_InLabStream* lab_stream;
	int debug;
	QNUInt32 start_offset;
	//QNUInt32 numsegs;
	unsigned long numsegs;
	QN_SegID segid;

public:
	CRF_FeatureStream(QN_InFtrStream* ftr, QN_InLabStream* lab, int dbg=0);
	CRF_FeatureStream(QN_InFtrStream* ftr, int dbg=0);
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
};

#endif /*CRF_FEATURESTREAM2_H_*/
