/*
 * CRF_InFtrStream_RandPresent.cpp
 *
 * Copyright (c) 2010
 * Author: Jeremy Morris
 *
 */
#include "CRF_InFtrStream_RandPresent.h"

/*
 * CRF_InFtrStream_RandPresent constructor
 *
 * Input: a_debug - debug flag
 *        a_classname - debug class name
 *        a_dbgname - debug filename
 *        a_str - input feature stream
 *        gen_type - presentation order
 *        seed - random seed for random presentation
 *
 */
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

/*
 *  CRF_InFtrStream_RandPresent destructor
 */
CRF_InFtrStream_RandPresent::~CRF_InFtrStream_RandPresent()
{
}

/*
 * CRF_InFtrStream_RandPresent::num_ftrs
 *
 * Returns: number of features in feature stream
 *
 * Accessor function to retreive number of features
 *
 */
size_t CRF_InFtrStream_RandPresent::num_ftrs()
{
	return str.num_ftrs();
}

/*
 * CRF_InFtrStream_RandPresent::nextseg
 *
 * Returns: segment id of the current segment
 *
 * Moves the internal pointer in the feature stream to the next segment.  returns the segment id or
 * an error if the last segment in the stream has been accessed.
 *
 */
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

/*
 * CRF_InFtrStream_RandPresent::read_ftrs
 *
 * Input: cnt - number of frames of features to read
 *        *ftrs - buffer to store read features in
 *
 * Returns: error code based on success or failure of feature read
 *
 * Reads cnt frames of input features and stores them in *ftrs buffer
 *
 */
size_t CRF_InFtrStream_RandPresent::read_ftrs(size_t cnt, float* ftrs)
{
	return str.read_ftrs(cnt, ftrs);
}

/*
 * CRF_InFtrStream_RandPresent::rewind
 *
 * Rewinds the file back to the first segment.
 *
 * for random presentation, the order of presentation changes with each epoch and the "first segment"
 * changes for each epoch.  Seed is used to ensure that even though the order is random, runs can be
 * reproduced.
 */
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

/*
 * CRF_InFtrStream_RandPresent::num_segs
 *
 * Returns number of segments in the stream
 */
size_t CRF_InFtrStream_RandPresent::num_segs()
{
	return str.num_segs();
}

/*
 * CRF_InFtrStream_RandPresent::num_frames
 *
 * Input: segno - segment number to be quered
 *
 * Returns: number of frames of feature vectors in segment segno
 *
 */
size_t CRF_InFtrStream_RandPresent::num_frames(size_t segno)
{
	return str.num_frames(segno);
}

/*
 * CRF_InFtrStream_RandPresent::get_pos
 *
 * Input: *segno - paramter to store segment number
 *        *frameno - parameter to store frame number
 *
 * Returns: the read segment id number for this feature stream
 *
 */
int CRF_InFtrStream_RandPresent::get_pos(size_t* segno, size_t* frameno)
{
	QN_SegID retval=str.get_pos(segno, frameno);
	*segno=real_seg;
	return this->real_seg;
}

/*
 * CRF_InFtrStream_RandPresent::set_pos
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

QN_SegID CRF_InFtrStream_RandPresent::set_pos(size_t segno, size_t frameno)
{
	return str.set_pos(segno, frameno);
}
