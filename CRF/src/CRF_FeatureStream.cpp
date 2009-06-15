#include "CRF_FeatureStream.h"

CRF_FeatureStream::CRF_FeatureStream(QN_InFtrStream* ftr, QN_InLabStream* lab, int dbg)
	: ftr_stream(ftr),
	  lab_stream(lab),
	  debug(dbg)
{
}

CRF_FeatureStream::CRF_FeatureStream(QN_InFtrStream* ftr, int dbg)
	: ftr_stream(ftr),
	  debug(dbg)
{
	this->lab_stream=NULL;
}

CRF_FeatureStream::~CRF_FeatureStream()
{
}

QN_SegID CRF_FeatureStream::nextseg()
{
	QN_SegID ftr_segid, lab_segid;
	
	ftr_segid = this->ftr_stream->nextseg();
	if (this->lab_stream == NULL) {
		lab_segid=ftr_segid;
	}
	else {
		lab_segid = this->lab_stream->nextseg();
	}
	//cout << "Feature ID: " << ftr_segid << " Label ID: " << lab_segid << endl;
	assert (ftr_segid == lab_segid);
	return ftr_segid;	
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
	return this->ftr_stream->get_pos(segno, frameno);
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
	this->ftr_stream->rewind();
	if (this->lab_stream != NULL) {
		this->lab_stream->rewind();
	}
}

QNUInt32 CRF_FeatureStream::num_ftrs()
{
	return this->ftr_stream->num_ftrs();
}

void CRF_FeatureStream::set_pos(QNUInt32 segno, QNUInt32 frmno)
{
	this->ftr_stream->set_pos(segno,frmno);
	if (this->lab_stream != NULL) {
		this->lab_stream->set_pos(segno,frmno);
	}
}

QNUInt32 CRF_FeatureStream::num_segs()
{
	return this->ftr_stream->num_segs();
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
