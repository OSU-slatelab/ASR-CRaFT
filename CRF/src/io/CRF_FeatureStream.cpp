/*
 * CRF_FeatureStream.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_FeatureStream.h"

// The following is to correct for the fact that QN_ALL ends up being defined as
// larger than a 32bit integer can handle on 64bit machines.  CRF_ALL is always the
// maximum a 32bit integer can hold.  This should be revisited and examined.
#define CRF_ALL 4294967295

// Modified by Ryan
// allowing for initialization of start_offset, numsegs, which is needed by the multi-threaded version
/*
 * CRF_FeatureStream constructor
 *
 * Input: *ftr - feature stream pre-created as a QN_InFtrStream object
 *        *lab - label stream pre-created as a QN_InLabStream object
 *        dbg - debugging flag
 *
 */
CRF_FeatureStream::CRF_FeatureStream(QN_InFtrStream* ftr, QN_InLabStream* lab, int dbg, QNUInt32 offset, QNUInt32 segs)
	: ftr_stream(ftr),
	  lab_stream(lab),
	  debug(dbg),
	  start_offset(offset),
	  numsegs(segs),
	  segid(QN_SEGID_BAD)
{
}

// Modified by Ryan
// allowing for initialization of start_offset, numsegs, which is needed by the multi-threaded version
/*
 * CRF_FeatureStream constructor
 *
 * Input: *ftr - feature stream pre-created as a QN_InFtrStream object
 *        dbg - debugging flag
 *
 * This constructor is used when a feature file has no associated label file
 * (i.e. during decoding).
 */
CRF_FeatureStream::CRF_FeatureStream(QN_InFtrStream* ftr, int dbg, QNUInt32 offset, QNUInt32 segs)
	: ftr_stream(ftr),
	  lab_stream(NULL),
	  debug(dbg),
	  start_offset(offset),
	  numsegs(segs),
	  segid(QN_SEGID_BAD)
{
}

/*
 * CRF_FeatureStream destructor
 */
CRF_FeatureStream::~CRF_FeatureStream()
{
}

/*
 * CRF_FeatureStream::nextseg
 *
 * Returns: Segment id of the current segment being addressed by the feature stream
 *
 * Advances the feature stream to the next segment.  If a label stream is defined for
 * the object, advances it as well.
 */
QN_SegID CRF_FeatureStream::nextseg()
{
	QN_SegID ftr_segid, lab_segid;

	if (this->segid == QN_SEGID_BAD) {
		this->segid=0;
	} else {
		this->segid++;
		//if ((this->numsegs != QN_ALL) &&
		if ((this->numsegs != ULONG_MAX) &&
			(this->segid >= this->numsegs)) {
			this->segid=QN_SEGID_BAD;
			return QN_SEGID_BAD;
		}
	}


	ftr_segid = this->ftr_stream->nextseg();
	if (this->lab_stream == NULL) {
		lab_segid=ftr_segid;
	} else {
		lab_segid = this->lab_stream->nextseg();
	}

	//cout << "Sizeof(size_t)" << sizeof(size_t) << "QN_ALL" << QN_ALL << endl;
	//cout << "Feature ID: " << ftr_segid << " Label ID: " << lab_segid << endl;
	assert (ftr_segid == lab_segid);
	//if (this->numsegs==QN_ALL && ftr_segid==QN_SEGID_BAD) {
	if (this->numsegs==CRF_ALL && ftr_segid==QN_SEGID_BAD) {
		this->segid=QN_SEGID_BAD;
	}
	return this->segid;
}

/*
 * CRF_FeatureStream::read
 *
 * Input:  bs - bunch size - number of frames of features and labels to read from the file
 *         fb - frame buffer to store read features
 *         lb - frame buffer to store read lable
 * Returns: Segment id of the current segment being addressed by the feature stream
 *
 * Advances the feature stream to the next segment.  If a label stream is defined for
 * the object, advances it as well.
 */
size_t CRF_FeatureStream::read(size_t bs, float* fb, QNUInt32* lb)
{
	size_t ftr_cnt, lab_cnt;

	ftr_cnt=this->ftr_stream->read_ftrs(bs,fb);
	if (this->lab_stream == NULL) {
		lab_cnt=ftr_cnt;
	}
	else {
		lab_cnt=this->lab_stream->read_labs(bs,lb);
	}
	//cout << "Feature Count: " << ftr_cnt << "Label Count: " << lab_cnt << endl;

	//commented out by Ryan, since ftr_cnt might be larger than lab_cnt for segmental features extraction
	//assert (ftr_cnt == lab_cnt);

	/*if (ftr_cnt != lab_cnt)
		ftr_cnt=-1;*/ // This doesn't work - size_t is unsigned - need a different indicator of failure here.
	return ftr_cnt;

}

/*
 * CRF_FeatureStream::get_pos
 *
 * Input: segno - parameter to store the current segment number
 *        frameno - parameter to store the current frame number
 *
 * Returns: Error code indicating success or failure in reading the current segment and
 *          frame number
 *
 */
int CRF_FeatureStream::get_pos(size_t* segno, size_t* frameno)
{
	size_t ftrsegno;
	*segno=this->segid;
	if (this->segid==QN_SEGID_BAD) {
		*frameno=0;
		return QN_SIZET_BAD;
	} else {
		return this->ftr_stream->get_pos(&ftrsegno, frameno);
	}
}

/*
 * CRF_FeatureStream::nextseg
 *
 * Input: in_stream - second independent feature stream object
 *
 * Returns: new feature stream object that is the concatenation of the current object
 *          and the object passed in by the parameter in_stream
 *
 * Used to build larger feature streams from multiple smaller files.  The streams must
 * match segment-wise and frame-wise for this function to work without error.
 *
 */
CRF_FeatureStream* CRF_FeatureStream::join(CRF_FeatureStream* in_stream)
{
	QN_InFtrStream* new_stream;
	new_stream = new QN_InFtrStream_JoinFtrs(this->debug, "join_ftrfile",
									*(this->ftr_stream), *(in_stream->ftr_stream));

	// Modified by Ryan
	// We need to keep the original start_offset, numsegs and segid so that
	// the multi-threaded version is not screwed.
	//
	//return new CRF_FeatureStream(new_stream, this->lab_stream);
	return new CRF_FeatureStream(new_stream, this->lab_stream, this->debug, this->start_offset, this->numsegs);
}

/*
 * CRF_FeatureStream::rewind
 *
 * Resets the stream back to the beginning.
 *
 * If start_offset is defined, "the beginning" is the segment determined by start_offset
 */
void CRF_FeatureStream::rewind()
{
#ifdef NON_VIEW_VERSION
	this->ftr_stream->rewind();
	if (this->lab_stream != NULL) {
		this->lab_stream->rewind();
	}
#else
	this->segid=QN_SEGID_BAD;
	if (this->start_offset==0) {
		this->ftr_stream->rewind();
		if (this->lab_stream != NULL) {
			this->lab_stream->rewind();
		}
	} else {
		QN_SegID ftrid= this->ftr_stream->set_pos(this->start_offset-1,0);
		if (ftrid==QN_SEGID_BAD) {
			// do it the hard way
			this->ftr_stream->rewind();
			for(QNUInt32 i=0;i<start_offset;i++)
				this->ftr_stream->nextseg();
		}
		if (this->lab_stream != NULL) {
			QN_SegID labid=this->lab_stream->set_pos(this->start_offset-1,0);
			if (labid==QN_SEGID_BAD) {
				// do it the hard way
				this->lab_stream->rewind();
				for(QNUInt32 i=0;i<start_offset;i++)
					this->lab_stream->nextseg();
			}
		}
	}
#endif

}

/*
 * CRF_FeatureStream::num_ftrs
 *
 * Returns: Number of features in the stream
 */

QNUInt32 CRF_FeatureStream::num_ftrs()
{
	return this->ftr_stream->num_ftrs();
}

/*
 * CRF_FeatureStream::set_pos
 *
 * Input: segno - segment number to jump to
 *        frmno - frame number to jump to
 *
 *
 * Advances the feature stream to the segment/frame defined by segno and frmno.
 */
void CRF_FeatureStream::set_pos(QNUInt32 segno, QNUInt32 frmno)
{
	//eeek!  better error handling needed
	//if ((this->numsegs != QN_ALL) &&
	if ((this->numsegs != CRF_ALL) &&
			(segno>this->start_offset+this->numsegs)) {
		cerr << "CRF_FeatureStream: set_pos(" << segno << "," << frmno << ") out of bounds in view" << endl;
		exit(1);
	}
	this->ftr_stream->set_pos(segno+this->start_offset,frmno);
	if (this->lab_stream != NULL) {
		this->lab_stream->set_pos(segno+this->start_offset,frmno);
	}
}

/*
 * CRF_FeatureStream::num_segs
 *
 * Returns: Total number of segments in the feature stream
 */
QNUInt32 CRF_FeatureStream::num_segs()
{
	//if (this->numsegs==QN_ALL) {
	if (this->numsegs==CRF_ALL) {
		return this->ftr_stream->num_segs();
	} else {
		return this->numsegs;
	}
}

/*
 * CRF_FeatureStream::display
 *
 * Dumps the content of the feature stream to the output stream.
 * Debugging function.
 *
 */
void CRF_FeatureStream::display() {

	QNUInt32 total_segs, ftr_count;
	QNUInt32 bunch_size=2;

	assert(bunch_size != 0);

	total_segs=this->num_segs();
	ftr_count=0;

	int seg_count=-1;
	int frm_count=-1;

	QNUInt32 num_ftrs = this->num_ftrs();
	cout << "Num Feas: " << num_ftrs << endl;

	float* ftr_buf = new float[num_ftrs*bunch_size];
	QNUInt32* lab_buf = new QNUInt32[bunch_size];


	while (1) {
		if (ftr_count < bunch_size) {
			QN_SegID ftr_segid;
			ftr_segid=this->nextseg();
//			ftr_segid=this->nextseg(instr,inlab);
			if (ftr_segid == QN_SEGID_BAD)
				break;
			seg_count++;
			frm_count=-1;
		}
		ftr_count=this->read(bunch_size,ftr_buf,lab_buf);
//		ftr_count=this->read(instr,inlab,bunch_size,ftr_buf,lab_buf);
		/*if (ftr_count != lab_count) {
			QN_OUTPUT("Feature Count %d / Lab Count %d mismatch",ftr_count,lab_count);
			break;
		}*/ // This check needs to be replaced with some kind of check on the read function - an error flag
		for (QNUInt32 i=0; i<ftr_count; i++) {
			frm_count++;
			cout << "Seg: " << seg_count << " Frame: " << frm_count;
			cout << " Label:" << lab_buf[i] << " Features: ";
			for (QNUInt32 j=0; j<num_ftrs; j++) {
				QNUInt32 idx=i*num_ftrs+j;
				cout << " " << ftr_buf[idx];
			}
			cout << endl;
		}
	}
}

/*
 * CRF_FeatureStream::view
 *
 * Input: startseg - segment to use as the initial segment for the view
 *        nsegs - total number of segments in the view
 *
 * Transforms the feature stream into a subset (or "view") of the overall feature stream.
 * A view contains a portion of the segments of the full stream in the same order as they
 * appear in the full stream.
 */
void CRF_FeatureStream::view(QNUInt32 startseg,QNUInt32 nsegs) {
	QNUInt32 ftrnseg=this->ftr_stream->num_segs();
	if (ftrnseg<startseg || startseg<0 ||
			(nsegs!=CRF_ALL && ftrnseg<startseg+nsegs)) {
		//(nsegs!=QN_ALL && ftrnseg<startseg+nsegs)) {
		cerr << "CRF_FeatureStream::view out of range " << startseg << " " << nsegs << endl;
		exit(1);
	}
	this->start_offset=startseg;
	this->numsegs=nsegs;
	this->rewind();
}


// Added by Ryan
/*
 * CRF_FeatureStream::num_labs
 *
 * Returns: the width of the label stream
 */
QNUInt32 CRF_FeatureStream::num_labs()
{
	if (this->lab_stream == NULL) {
		return 0;
	}
	else {
		return this->lab_stream->num_labs();
	}
}
