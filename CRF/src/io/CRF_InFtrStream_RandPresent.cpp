#include "CRF_InFtrStream_RandPresent.h"

CRF_InFtrStream_RandPresent::CRF_InFtrStream_RandPresent(int a_debug, 
														const char* a_classname, 
														const char* a_dbgname,
														QN_InFtrStream& a_str,
														seqtype gen_type,
														QNUInt32 seed)
		: log(a_debug, a_classname, a_dbgname),
		  str(a_str),
		  seqgen_type(gen_type)
{
	this->seed=seed;
	//seed=0;
	//seqgen=new QN_SeqGen_Random48(a_str.num_segs(),seed);
	//seqgen=new QN_SeqGen_RandomNoReplace(a_str.num_segs(),seed);
	//cur_seg=0;
	epoch=0;
	max_segs=a_str.num_segs();
	seqgen=NULL;
	this->rewind();
}

CRF_InFtrStream_RandPresent::~CRF_InFtrStream_RandPresent()
{
}

size_t CRF_InFtrStream_RandPresent::num_ftrs()
{
	return str.num_ftrs();
}

QN_SegID CRF_InFtrStream_RandPresent::nextseg()
{
	if (cur_seg>=max_segs){
		return QN_SEGID_BAD;
	}
	else {
		cur_seg++;
		this->real_seg=seqgen->next();
		//std::cout << "Setting sequence to segment " << this->real_seg << std::endl;
		QN_SegID retval = str.set_pos(this->real_seg,0);
		if (retval != QN_SEGID_BAD) {
			return this->real_seg;
		}
		else {
			return retval;
		}
	}
	//return str.nextseg();
}

size_t CRF_InFtrStream_RandPresent::read_ftrs(size_t cnt, float* ftrs)
{
	return str.read_ftrs(cnt, ftrs);
}

int CRF_InFtrStream_RandPresent::rewind()
{
	cur_seg=0;
	epoch++;
	if (seqgen != NULL)
		delete seqgen;
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

size_t CRF_InFtrStream_RandPresent::num_segs()
{
	return str.num_segs();
}

size_t CRF_InFtrStream_RandPresent::num_frames(size_t segno)
{
	return str.num_frames(segno);
}


int CRF_InFtrStream_RandPresent::get_pos(size_t* segno, size_t* frameno)
{
	QN_SegID retval=str.get_pos(segno, frameno);
	*segno=real_seg;
	return this->real_seg;
}

QN_SegID CRF_InFtrStream_RandPresent::set_pos(size_t segno, size_t frameno)
{
	return str.set_pos(segno, frameno);
}
