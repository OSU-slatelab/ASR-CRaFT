#include "CRF_InLabStream_RandPresent.h"

CRF_InLabStream_RandPresent::CRF_InLabStream_RandPresent(int a_debug, 
														 const char* a_classname, 
														 const char* a_dbgname,
														 QN_InLabStream& a_str,
														 seqtype gentype, QNUInt32 seed)
	:	log(a_debug, a_classname, a_dbgname),
		str(a_str),
		seqgen_type(gentype)
{
	this->seed=seed;
	//seed=0;
	//seqgen=new QN_SeqGen_Random48(a_str.num_segs(),seed);
	//seqgen=new QN_SeqGen_RandomNoReplace(a_str.num_segs(),seed);
	//cur_seg=0;
	epoch=0;
	max_segs=a_str.num_segs();
	this->rewind();
}

CRF_InLabStream_RandPresent::~CRF_InLabStream_RandPresent()
{
}

size_t CRF_InLabStream_RandPresent::num_labs()
{
	return str.num_labs();
}

QN_SegID CRF_InLabStream_RandPresent::nextseg()
{
	if (cur_seg>=max_segs){
		return QN_SEGID_BAD;
	}
	else {
		cur_seg++;
		this->real_seg=seqgen->next();
		QN_SegID retval = str.set_pos(this->real_seg,0);
		if (retval != QN_SEGID_BAD) {
			return this->real_seg;
		}
		else {
			return retval;
		}
	}
}

size_t CRF_InLabStream_RandPresent::read_labs(size_t cnt, QNUInt32* labs)
{
	return str.read_labs(cnt, labs);
}

size_t CRF_InLabStream_RandPresent::read_labs(size_t cnt, QNUInt32* labs, size_t stride)
{
	return str.read_labs(cnt, labs, stride);
}

int CRF_InLabStream_RandPresent::rewind()
{
	cur_seg=0;
	epoch++;
	//seqgen=new QN_SeqGen_Random48(str.num_segs(),12345*epoch+seed);
	//seqgen=new QN_SeqGen_RandomNoReplace(str.num_segs(),12345*epoch+seed);
	switch ( this->seqgen_type )
	{
		case RANDOM_NO_REPLACE:
			seqgen=new QN_SeqGen_RandomNoReplace(str.num_segs(),12345*epoch+seed);
			break;
		case RANDOM_REPLACE:
			seqgen=new QN_SeqGen_Random48(str.num_segs(),12345*epoch+seed);
			break;
		case SEQUENTIAL:
			exit(-1);
			break;
		default:
			exit(-1);
	}
	
	return str.rewind();
}

size_t CRF_InLabStream_RandPresent::num_segs()
{
	return str.num_segs();
}

size_t CRF_InLabStream_RandPresent::num_frames(size_t segno)
{
	return str.num_frames(segno);
}

int CRF_InLabStream_RandPresent::get_pos(size_t* segno, size_t* frameno)
{
	return str.get_pos(segno, frameno);
}

QN_SegID CRF_InLabStream_RandPresent::set_pos(size_t segno, size_t frameno)
{
	return str.set_pos(segno, frameno);
}
