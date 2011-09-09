/*
 * CRF_InFtrStream_SeqMultiWindow.cpp
 *
 *  Created on: Aug 19, 2011
 *      Author: hey
 */

#include "CRF_InFtrStream_SeqMultiWindow.h"
#include <math.h>

CRF_InFtrStream_SeqMultiWindow::CRF_InFtrStream_SeqMultiWindow(int a_debug, const char* a_dbgname,
        QN_InFtrStream& a_str,
        size_t a_max_win_len, size_t a_top_margin,
        size_t a_bot_margin,
        size_t a_bunch_frames)
	: log(a_debug, "CRF_InFtrStream_SeqMultiWindow", a_dbgname),
	  in_str(a_str),
	  max_win_len(a_max_win_len),
	  top_margin(a_top_margin),
	  bot_margin(a_bot_margin),
	  in_width(a_str.num_ftrs()),
	  max_win_in_width(max_win_len * in_width),
	  bunch_frames(a_bunch_frames==QN_SIZET_BAD
	             ? top_margin + max_win_len + bot_margin
	             : a_bunch_frames),
	  max_buf_lines(bunch_frames + max_win_len + bot_margin-1),
	  max_win_buf(new float[max_buf_lines*in_width]),
	  first_line_ptr(&max_win_buf[top_margin*in_width]),
	  segno(-1)
{
	out_width = 488 + max_win_len;  //TODO: Currently hard-coded. Need to be parameterized.

	if (bunch_frames<(top_margin+bot_margin+max_win_len))
		log.error("Insufficient bunch length for size of window.");
	log.log(QN_LOG_PER_RUN, "Windowing, window=%lu*%lu values, "
	        "top margin=%lu frames, bottom margin=%lu frames, "
	        "buffer size=%lu frames.", (unsigned long) in_width,
	        (unsigned long) max_win_len, (unsigned long) top_margin,
	        (unsigned long) bot_margin, (unsigned long) max_buf_lines);
}

CRF_InFtrStream_SeqMultiWindow::~CRF_InFtrStream_SeqMultiWindow() {
	delete [] max_win_buf;
}

size_t CRF_InFtrStream_SeqMultiWindow::num_ftrs()
{
	return out_width;
}

QN_SegID CRF_InFtrStream_SeqMultiWindow::nextseg()
{
	QN_SegID segid;             // Returned segment ID.

	segid = in_str.nextseg();
	if (segid!=QN_SEGID_BAD)
	{
	    // Try and read in a buffer full of frames.
	    buf_lines = in_str.read_ftrs(bunch_frames, max_win_buf);
	    cur_line = top_margin;
	    cur_line_ptr = first_line_ptr;

	    cur_avail_max_win_len = 1;

	    //cout << "top_margin:" << top_margin << ", bot_margin:" << bot_margin << ", buf_lines:" << buf_lines << endl;

	    // Count how many frames are available in the longest window for this sentence.
	    if (top_margin + bot_margin >= buf_lines)   // No frames available.
	    {
	    	seg_max_win_len = 0;
	    }
	    else if (top_margin + max_win_len + bot_margin <= buf_lines)  // More than max_win_len frames available.
	    {
	    	seg_max_win_len = max_win_len;
	    }
	    else                          // Less than max_win_len frames available.
	    {
	    	seg_max_win_len = buf_lines - top_margin - bot_margin;
	    }

	    segno++;
	}
	log.log(QN_LOG_PER_SENT, "Moved on to segment %li.", segno);
	return segid;
}

size_t CRF_InFtrStream_SeqMultiWindow::read_ftrs(size_t cnt, float* ftrs)
{
	//return read_ftrs(cnt, ftrs, QN_ALL);

	if (cnt != 1)
		log.error("Batch size has to be 1: every time read_ftrs() can just move forward by 1 frame.");

	if (segno==-1)
		log.error("Trying to read before start of first sentence.");

	if (cur_avail_max_win_len > seg_max_win_len)
		return 0;

	if (cur_avail_max_win_len == max_win_len)
	{
		size_t frame;               // Current frame being transferred.
		const size_t l_in_width = in_width; // Local copy of class variable.
		for (frame=0; frame<cnt; frame++)
		{
			// If at end of buffer, read in more lines.
			if (cur_line + max_win_len + bot_margin > buf_lines)
			{
				// Move the remaining lines to the beginning of the buffer.
				size_t line;        // Iterator over lines being moved.
				size_t old_lines = max_win_len + bot_margin - 1; // Lines remaining.
				float *from = cur_line_ptr; // Where we get the line from.
				float *to = &max_win_buf[0]; // Where we put the line.
				for (line=0; line<old_lines; line++)
				{
					qn_copy_vf_vf(l_in_width, from, to);
					from += l_in_width;
					to += l_in_width;
				}
				cur_line = 0;
				cur_line_ptr = &max_win_buf[0];

				// Try and fill up the rest of the buffer.
				size_t count = in_str.read_ftrs(bunch_frames,
											  &max_win_buf[old_lines*l_in_width]);
				buf_lines = old_lines + count;
				// If nothing to read in, terminate read.
				if (count==0)
				{
					return 0;
				}
			}
		}
	}

	size_t multi_win_count = cur_avail_max_win_len;

	float* vals = ftrs;       // Where to put windowed data for each individual features.

	vals = sample_ftrs(cur_line_ptr, vals, multi_win_count);
	vals = avg_ftrs(cur_line_ptr, vals, multi_win_count);
	vals = max_ftrs(cur_line_ptr, vals, multi_win_count);
	vals = min_ftrs(cur_line_ptr, vals, multi_win_count);
	vals = dur_ftrs(vals, multi_win_count);

	if (cur_avail_max_win_len < max_win_len)   //for leading frames
	{
		cur_avail_max_win_len++;
	}
	else                                       //for all following frames
	{
		cur_line += 1;
		cur_line_ptr += in_width;
	}

	log.log(QN_LOG_PER_BUNCH, "Read %lu windows for 1 frame.", (unsigned long) multi_win_count);

	return multi_win_count;
}

size_t CRF_InFtrStream_SeqMultiWindow::read_ftrs(size_t cnt, float* ftrs, size_t stride)
{
	size_t ret = QN_SIZET_BAD;
	if (stride==num_ftrs() || stride==QN_ALL)
	{
		ret = read_ftrs(cnt, ftrs);
	}
	else
	{
		log.error("The strided version of the read_ftrs function currently does not work for CRF_InFtrStream_SeqMultiWindow.");
	}
	return QN_SIZET_BAD;
}

// Rewind works
int CRF_InFtrStream_SeqMultiWindow::rewind()
{
	int ec;

	ec = in_str.rewind();
	if (ec==QN_OK)
		segno = -1;
	return ec;
}

// Some streams functions that do not work for SeqMultiWindows.

size_t CRF_InFtrStream_SeqMultiWindow::num_segs()
{
	return QN_SIZET_BAD;
}

size_t CRF_InFtrStream_SeqMultiWindow::num_frames(size_t segno)
{
	return QN_SIZET_BAD;
}

int CRF_InFtrStream_SeqMultiWindow::get_pos(size_t* segno, size_t* frameno)
{
	return QN_BAD;
}

QN_SegID CRF_InFtrStream_SeqMultiWindow::set_pos(size_t segno, size_t frameno)
{
	return QN_SEGID_BAD;
}


// Each individual window-based feature function

float* CRF_InFtrStream_SeqMultiWindow::sample_ftrs(float* in_win_buf, float* ftrs, size_t avail_max_win_len)
{
	//for win_len=1
	float* cur_win_buf_ptr = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_ftrs_ptr = ftrs;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		cout << "cur_win_len=" << cur_win_len << " frameSteps:";
		float* cur_sample_frame_buf_ptr = cur_win_buf_ptr;
		float* cur_sample_frame_ftrs_ptr = cur_win_ftrs_ptr;
		float one_tenth_win_len = cur_win_len * 0.1;
		for (int i = 1; i < 10; i += 2)  //sample frames at 10%, 30%, 50%, 70% and 90% points
		{
			QNUInt32 frameStep = (QNUInt32)ceil(one_tenth_win_len * i) - 1;   //it might be not accurate for some integral value since it is floating point calculation
			//QNUInt32 frameStep = (QNUInt32)(one_tenth_win_len * i);   //it might be not accurate for some integral value since it is floating point calculation
			cur_sample_frame_buf_ptr += in_width * frameStep;
			if (cur_sample_frame_ftrs_ptr != NULL)
			{
				qn_copy_vf_vf(in_width, cur_sample_frame_buf_ptr, cur_sample_frame_ftrs_ptr);
				cur_sample_frame_ftrs_ptr += in_width;
			}
			cur_sample_frame_buf_ptr = cur_win_buf_ptr;
			cout << " " << frameStep;
		}
		cur_win_buf_ptr -= in_width;
		cur_win_ftrs_ptr += out_width;
		cout << endl;
	}

	return ftrs + 5 * in_width;	    //TODO: Currently hard-coded. Need to be parameterized.
}

float* CRF_InFtrStream_SeqMultiWindow::avg_ftrs(float* in_win_buf, float* ftrs, size_t avail_max_win_len)
{
	//for win_len=1
	float* cur_win_buf_ptr = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_ftrs_ptr = ftrs;

	float* acc_sum_ftrs(new float[in_width]);
	for (size_t acc_i = 0; acc_i < in_width; acc_i++)
	{
		acc_sum_ftrs[acc_i] = 0.0f;
	}

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
		{
			acc_sum_ftrs[ftr_i] += cur_win_buf_ptr[ftr_i];
			cur_win_ftrs_ptr[ftr_i] = acc_sum_ftrs[ftr_i] / cur_win_len;
		}
		cur_win_buf_ptr -= in_width;
		cur_win_ftrs_ptr += out_width;
	}

	delete [] acc_sum_ftrs;
	return ftrs + in_width;        //TODO: Currently hard-coded. Need to be parameterized.
}

float* CRF_InFtrStream_SeqMultiWindow::max_ftrs(float* in_win_buf, float* ftrs, size_t avail_max_win_len)
{
	//for win_len=1
	float* cur_win_buf_ptr = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_ftrs_ptr = ftrs;

	float* acc_max_ftrs(new float[in_width]);
	for (size_t acc_i = 0; acc_i < in_width; acc_i++)
	{
		acc_max_ftrs[acc_i] = cur_win_buf_ptr[acc_i];
	}

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
		{
			if (cur_win_buf_ptr[ftr_i] > acc_max_ftrs[ftr_i])
			{
				acc_max_ftrs[ftr_i] = cur_win_buf_ptr[ftr_i];
			}
			cur_win_ftrs_ptr[ftr_i] = acc_max_ftrs[ftr_i];
		}
		cur_win_buf_ptr -= in_width;
		cur_win_ftrs_ptr += out_width;
	}

	delete [] acc_max_ftrs;
	return ftrs + in_width;       //TODO: Currently hard-coded. Need to be parameterized.
}

float* CRF_InFtrStream_SeqMultiWindow::min_ftrs(float* in_win_buf, float* ftrs, size_t avail_max_win_len)
{
	//for win_len=1
	float* cur_win_buf_ptr = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_ftrs_ptr = ftrs;

	float* acc_min_ftrs(new float[in_width]);
	for (size_t acc_i = 0; acc_i < in_width; acc_i++)
	{
		acc_min_ftrs[acc_i] = cur_win_buf_ptr[acc_i];
	}

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
		{
			if (cur_win_buf_ptr[ftr_i] < acc_min_ftrs[ftr_i])
			{
				acc_min_ftrs[ftr_i] = cur_win_buf_ptr[ftr_i];
			}
			cur_win_ftrs_ptr[ftr_i] = acc_min_ftrs[ftr_i];
		}
		cur_win_buf_ptr -= in_width;
		cur_win_ftrs_ptr += out_width;
	}

	delete [] acc_min_ftrs;
	return ftrs + in_width;      //TODO: Currently hard-coded. Need to be parameterized.
}

float* CRF_InFtrStream_SeqMultiWindow::dur_ftrs(float* ftrs, size_t avail_max_win_len)
{
	//for win_len=1
	float* cur_win_ftrs_ptr = ftrs;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (QNUInt32 tmp_win_len=1; tmp_win_len <= max_win_len; tmp_win_len++)
		{
			if (tmp_win_len == cur_win_len)
			{
				cur_win_ftrs_ptr[tmp_win_len-1] = 1;
			}
			else
			{
				cur_win_ftrs_ptr[tmp_win_len-1] = 0;
			}
		}
		cur_win_ftrs_ptr += out_width;
	}

	return ftrs + max_win_len;   //TODO: Currently hard-coded. Need to be parameterized.
}
