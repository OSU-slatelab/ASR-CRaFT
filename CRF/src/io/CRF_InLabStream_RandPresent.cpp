/*
 * CRF_InLabStream_RandPresent.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_InLabStream_RandPresent.h"

/*
 * CRF_InLabStream_RandPresent constructor
 *
 * Input: a_debug - debug flag
 *        a_classname - debug class name
 *        a_dbgname - debug filename
 *        a_str - input feature stream
 *        gen_type - presentation order
 *        seed - random seed for random presentation
 *
 */
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

/*
 *  CRF_InLabStream_RandPresent destructor
 */
CRF_InLabStream_RandPresent::~CRF_InLabStream_RandPresent()
{
}

/*
 * CRF_InLabStream_RandPresent::num_labs
 *
 * Returns: number of labels in feature stream
 *
 * Accessor function to retreive number of labels
 *
 */
size_t CRF_InLabStream_RandPresent::num_labs()
{
	return str.num_labs();
}

/*
 * CRF_InLabStream_RandPresent::nextseg
 *
 * Returns: segment id of the current segment
 *
 * Moves the internal pointer in the label stream to the next segment.  returns the segment id or
 * an error if the last segment in the stream has been accessed.
 *
 */
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

/*
 * CRF_InLabStream_RandPresent::read_ftrs
 *
 * Input: cnt - number of frames of features to read
 *        *labs - buffer to store read labels in
 *
 * Returns: error code based on success or failure of feature read
 *
 * Reads cnt frames of input labels and stores them in *labs buffer
 *
 */
size_t CRF_InLabStream_RandPresent::read_labs(size_t cnt, QNUInt32* labs)
{
	return str.read_labs(cnt, labs);
}


size_t CRF_InLabStream_RandPresent::read_labs(size_t cnt, QNUInt32* labs, size_t stride)
{
	return str.read_labs(cnt, labs, stride);
}

/*
 * CRF_InLabStream_RandPresent::rewind
 *
 * Rewinds the file back to the first segment.
 *
 * for random presentation, the order of presentation changes with each epoch and the "first segment"
 * changes for each epoch.  Seed is used to ensure that even though the order is random, runs can be
 * reproduced.
 */
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

/*
 * CRF_InLabStream_RandPresent::num_segs
 *
 * Returns number of segments in the stream
 */
size_t CRF_InLabStream_RandPresent::num_segs()
{
	return str.num_segs();
}

/*
 * CRF_InLabStream_RandPresent::num_frames
 *
 * Input: segno - segment number to be quered
 *
 * Returns: number of frames of feature vectors in segment segno
 *
 */
size_t CRF_InLabStream_RandPresent::num_frames(size_t segno)
{
	return str.num_frames(segno);
}

/*
 * CRF_InLabStream_RandPresent::get_pos
 *
 * Input: *segno - paramter to store segment number
 *        *frameno - parameter to store frame number
 *
 * Returns: the read segment id number for this feature stream
 *
 */
int CRF_InLabStream_RandPresent::get_pos(size_t* segno, size_t* frameno)
{
	return str.get_pos(segno, frameno);
}

/*
 * CRF_InLabStream_RandPresent::set_pos
 *
 * Input: segno - segment number to jump to
 *        frameno - frame number to jump to
 *
 * Returns: Segment number jumped to or error code
 *
 * Sets the position in the feature stream to the segment and frame defined by the segno and frameno
 *  parameters.
 *
 */
QN_SegID CRF_InLabStream_RandPresent::set_pos(size_t segno, size_t frameno)
{
	return str.set_pos(segno, frameno);
}
