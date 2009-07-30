#include "CRF_FeatureStream.h"

CRF_FeatureStream::CRF_FeatureStream(QN_InFtrStream* ftr, QN_InLabStream* lab, int dbg)
	: ftr_stream(ftr),
	  lab_stream(lab),
	  debug(dbg),
	  start_offset(0),
	  numsegs(QN_ALL),
	  segid(QN_SEGID_BAD)
{
}

CRF_FeatureStream::CRF_FeatureStream(QN_InFtrStream* ftr, int dbg)
	: ftr_stream(ftr),
	  lab_stream(NULL),
	  debug(dbg),
	  start_offset(0),
	  numsegs(QN_ALL),
	  segid(QN_SEGID_BAD)
{
}

CRF_FeatureStream::~CRF_FeatureStream()
{
}

QN_SegID CRF_FeatureStream::nextseg()
{
	QN_SegID ftr_segid, lab_segid;

	if (this->segid == QN_SEGID_BAD) {
		this->segid=0;
	} else {
		this->segid++;
		if ((this->numsegs != QN_ALL) &&
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
	if (this->numsegs==QN_ALL && ftr_segid==QN_SEGID_BAD) {
		this->segid=QN_SEGID_BAD;
	}
	return this->segid;
}

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
	assert (ftr_cnt == lab_cnt);
	/*if (ftr_cnt != lab_cnt)
		ftr_cnt=-1;*/ // This doesn't work - size_t is unsigned - need a different indicator of failure here.
	return ftr_cnt;

}

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

CRF_FeatureStream* CRF_FeatureStream::join(CRF_FeatureStream* in_stream)
{
	QN_InFtrStream* new_stream;
	new_stream = new QN_InFtrStream_JoinFtrs(this->debug, "join_ftrfile",
									*(this->ftr_stream), *(in_stream->ftr_stream));
	return new CRF_FeatureStream(new_stream,this->lab_stream);
}

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


QNUInt32 CRF_FeatureStream::num_ftrs()
{
	return this->ftr_stream->num_ftrs();
}

void CRF_FeatureStream::set_pos(QNUInt32 segno, QNUInt32 frmno)
{
	//eeek!  better error handling needed
	if ((this->numsegs != QN_ALL) &&
			(segno>this->start_offset+this->numsegs)) {
		cerr << "CRF_FeatureStream: set_pos(" << segno << "," << frmno << ") out of bounds in view" << endl;
		exit(1);
	}
	this->ftr_stream->set_pos(segno+this->start_offset,frmno);
	if (this->lab_stream != NULL) {
		this->lab_stream->set_pos(segno+this->start_offset,frmno);
	}
}

QNUInt32 CRF_FeatureStream::num_segs()
{
	if (this->numsegs==QN_ALL) {
		return this->ftr_stream->num_segs();
	} else {
		return this->numsegs;
	}
}

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

void CRF_FeatureStream::view(QNUInt32 startseg,QNUInt32 nsegs) {
	QNUInt32 ftrnseg=this->ftr_stream->num_segs();
	if (ftrnseg<startseg || startseg<0 ||
		(nsegs!=QN_ALL && ftrnseg<startseg+nsegs)) {
		cerr << "CRF_FeatureStream::view out of range " << startseg << " " << nsegs << endl;
		exit(1);
	}
	this->start_offset=startseg;
	this->numsegs=nsegs;
	this->rewind();
}
