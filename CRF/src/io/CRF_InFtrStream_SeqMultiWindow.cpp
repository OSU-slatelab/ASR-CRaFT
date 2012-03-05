/*
 * CRF_InFtrStream_SeqMultiWindow.cpp
 *
 *  Created on: Aug 19, 2011
 *      Author: hey
 */

#include "CRF_InFtrStream_SeqMultiWindow.h"

// modified for context features
//CRF_InFtrStream_SeqMultiWindow::CRF_InFtrStream_SeqMultiWindow(int a_debug, const char* a_dbgname,
//        QN_InFtrStream& a_str,
//        size_t a_max_win_len, size_t a_top_margin,
//        size_t a_bot_margin,
//        size_t a_bunch_frames)
CRF_InFtrStream_SeqMultiWindow::CRF_InFtrStream_SeqMultiWindow(int a_debug, const char* a_dbgname,
		QN_InFtrStream& a_str,
		size_t a_max_win_len, size_t a_top_margin,
		size_t a_bot_margin, size_t a_left_ctx_len,
		size_t a_right_ctx_len, bool a_seg_ftr,
		bool a_use_boundary_delta_ftr,
		size_t a_bunch_frames)
	: log(a_debug, "CRF_InFtrStream_SeqMultiWindow", a_dbgname),
	  in_str(a_str),
	  max_win_len(a_max_win_len),
	  top_margin(a_top_margin),
	  bot_margin(a_bot_margin),

	  // added for context features
	  left_context_len(a_left_ctx_len),
	  right_context_len(a_right_ctx_len),

	  in_width(a_str.num_ftrs()),
	  max_win_in_width(max_win_len * in_width),

	  // modified for context features
//	  bunch_frames(a_bunch_frames==QN_SIZET_BAD
//	             ? top_margin + max_win_len + bot_margin
//	             : a_bunch_frames),
//	  max_buf_lines(bunch_frames + max_win_len + bot_margin-1),
	  bunch_frames(a_bunch_frames==QN_SIZET_BAD
				 ? top_margin + left_context_len + max_win_len + right_context_len + bot_margin
				 : a_bunch_frames),
	  max_buf_lines(bunch_frames + left_context_len + max_win_len +
			  right_context_len + bot_margin-1),
	  extract_segment_features(a_seg_ftr),
	  use_boundary_delta_ftrs(a_use_boundary_delta_ftr),

	  max_win_buf(new float[max_buf_lines*in_width]),

	  // modified for context features
//	  first_line_ptr(&max_win_buf[top_margin*in_width]),
	  first_line_ptr(&max_win_buf[(top_margin+left_context_len)*in_width]),

	  segno(-1)
{
	if (max_win_len == 0)
		log.error("Window size must be larger than 0.");

	if (max_win_len == 1)
	{
		// modified for context features
		//out_width = in_width;
		out_width = (left_context_len + 1 + right_context_len) * in_width;
	} else {

		// modified for context features
		if (extract_segment_features)
		{
			//out_width = 8 * in_width + max_win_len;   //TODO: Currently hard-coded. Need to be parameterized.
			////out_width = 488 + max_win_len;
			out_width = 8 * in_width + max_win_len +
					(left_context_len + right_context_len) * in_width;  //TODO: Currently hard-coded. Need to be parameterized.
//			out_width = 5 * in_width + max_win_len +
//								(left_context_len + right_context_len) * in_width;  //TODO: Currently hard-coded. Need to be parameterized.
//			out_width = in_width + max_win_len +
//								(left_context_len + right_context_len) * in_width;  //TODO: Currently hard-coded. Need to be parameterized.
		} else {
			if (use_boundary_delta_ftrs)
			{
				QNUInt32 both_context_len = left_context_len;
				if (both_context_len > right_context_len + 1)
				{
					both_context_len = right_context_len + 1;
				}
				out_width = both_context_len * in_width;
			}
			else
			{
				out_width = (left_context_len + 1 + right_context_len) * in_width;
			}
		}
	}

	if (this->extract_segment_features && this->use_boundary_delta_ftrs)
	{
		log.error("extract_segment_features must be false to use boundary_delta_ftrs.");
	}

	// modified for context features
//	if (bunch_frames<(top_margin+bot_margin+max_win_len))
//		log.error("Insufficient bunch length for size of window.");
//	log.log(QN_LOG_PER_RUN, "Windowing, window=%lu*%lu values, "
//	        "top margin=%lu frames, bottom margin=%lu frames, "
//	        "buffer size=%lu frames.", (unsigned long) out_width,
//	        (unsigned long) max_win_len, (unsigned long) top_margin,
//	        (unsigned long) bot_margin, (unsigned long) max_buf_lines);
	if (bunch_frames<(top_margin + left_context_len + max_win_len +
			 right_context_len + bot_margin))
			log.error("bunch_frame=%lu < %lu, Insufficient bunch length for size of window.",
					(unsigned long) bunch_frames,
					(unsigned long) top_margin + left_context_len + max_win_len + right_context_len + bot_margin);
	log.log(QN_LOG_PER_RUN, "Windowing, window=%lu*%lu values, "
			"top margin=%lu frames, bottom margin=%lu frames, "
			"left context=%lu frames, right context=%lu frames, "
			"buffer size=%lu frames.", (unsigned long) out_width,
			(unsigned long) max_win_len, (unsigned long) top_margin,
			(unsigned long) bot_margin, (unsigned long) left_context_len,
			(unsigned long) right_context_len, (unsigned long) max_buf_lines);
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

	    //modified for context features
//	    cur_line = top_margin;
	    cur_line = top_margin + left_context_len;

		cur_line_ptr = first_line_ptr;

	    cur_avail_max_win_len = 1;

	    // just for debugging
//	    cout << "top_margin:" << top_margin << ", bot_margin:" << bot_margin << ", buf_lines:" << buf_lines << endl;

	    // Count how many frames are available in the longest window for this sentence.
	    // modified for context features
//	    if (top_margin + bot_margin >= buf_lines)   // No frames available.
//	    {
//	    	seg_max_win_len = 0;
//	    }
//	    else if (top_margin + max_win_len + bot_margin <= buf_lines)  // More than max_win_len frames available.
//	    {
//	    	seg_max_win_len = max_win_len;
//	    }
//	    else                          // Less than max_win_len frames available.
//	    {
//	    	seg_max_win_len = buf_lines - top_margin - bot_margin;
//	    }
	    if (top_margin + left_context_len + right_context_len + bot_margin >= buf_lines)   // No frames available.
		{
			seg_max_win_len = 0;
		}
		else if (top_margin + left_context_len + max_win_len +
				right_context_len + bot_margin <= buf_lines)  // More than max_win_len frames available.
		{
			seg_max_win_len = max_win_len;
		}
		else                          // Less than max_win_len frames available.
		{
			seg_max_win_len = buf_lines - top_margin - bot_margin - left_context_len - right_context_len;
		}

	    segno++;
	}
	log.log(QN_LOG_PER_SENT, "Moved on to segment %li.", segno);
	return segid;
}

size_t CRF_InFtrStream_SeqMultiWindow::read_ftrs(size_t cnt, float* ftrs)
{
	return read_ftrs(cnt, ftrs, QN_ALL);

////	if (cnt != 1)
////		log.error("Batch size has to be 1: every time read_ftrs() can just move forward by 1 frame.");
//	if (cnt > cur_avail_max_win_len)
//		log.error("The number of windows(" << cnt << ") is set too large, which has to be " << cur_avail_max_win_len << " at current frame.");
//	if (cnt < cur_avail_max_win_len)
//		log.error("The number of windows(" << cnt << ") is set too small, which has to be " << cur_avail_max_win_len << " at current frame.");
//
//	if (segno==-1)
//		log.error("Trying to read before start of first sentence.");
//
//	if (cur_avail_max_win_len > seg_max_win_len)
//		return 0;
//
//	if (cur_avail_max_win_len == max_win_len)
//	{
//		const size_t l_in_width = in_width; // Local copy of class variable.
//
//		// the max window moves forward by one frame
//
//		// If at end of buffer, read in more lines.
//		// modified for context features
////			if (cur_line + max_win_len + bot_margin > buf_lines)
//		if (cur_line + max_win_len + right_context_len + bot_margin > buf_lines)
//		{
//			// Move the remaining lines to the beginning of the buffer.
//			size_t line;        // Iterator over lines being moved.
//
//			// modified for context features
////				size_t old_lines = max_win_len + bot_margin - 1; // Lines remaining.
//			size_t old_lines = left_context_len + max_win_len +
//					  right_context_len + bot_margin - 1; // Lines remaining.
//
//			// modified for context features
////				float *from = cur_line_ptr; // Where we get the line from.
//			float *from = cur_line_ptr - left_context_len * in_width; // Where we get the line from.
//
//			float *to = &max_win_buf[0]; // Where we put the line.
//			for (line=0; line<old_lines; line++)
//			{
//				qn_copy_vf_vf(l_in_width, from, to);
//				from += l_in_width;
//				to += l_in_width;
//			}
//
//			//modified for context features
////				cur_line = 0;
////				cur_line_ptr = &max_win_buf[0];
//			cur_line = left_context_len;
//			cur_line_ptr = &max_win_buf[left_context_len * in_width];
//
//			// Try and fill up the rest of the buffer.
//			size_t count = in_str.read_ftrs(bunch_frames,
//										  &max_win_buf[old_lines*l_in_width]);
//			buf_lines = old_lines + count;
//			// If nothing to read in, terminate read.
//			if (count==0)
//			{
//				return 0;
//			}
//		}
//	}
//
////	size_t multi_win_count = cur_avail_max_win_len;
//	size_t multi_win_count = cnt;
//
//	size_t numWrittenFtrPerWin = 0;   // the total number of features that have been written to each output window.
//
//	// added for context features
//	// left context features
//	numWrittenFtrPerWin += first_frame_left_ctx_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//
//	if (max_win_len == 1)
//	{
//		// modified for context features
////		numWrittenFtrPerWin += avg_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count);
//		numWrittenFtrPerWin += first_frame_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//	} else {
//		// modified for context features
////		numWrittenFtrPerWin += sample_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count);
////		numWrittenFtrPerWin += avg_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count);
////		numWrittenFtrPerWin += max_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count);
////		numWrittenFtrPerWin += min_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count);
////		numWrittenFtrPerWin += dur_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count);
//		if (extract_segment_features)
//		{
//			numWrittenFtrPerWin += sample_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//			numWrittenFtrPerWin += avg_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//			numWrittenFtrPerWin += max_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//			numWrittenFtrPerWin += min_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//			numWrittenFtrPerWin += dur_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//		} else {
//			numWrittenFtrPerWin += first_frame_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//		}
//	}
//
//	// added for context features
//	// right context features
//	if (extract_segment_features)
//	{
//		numWrittenFtrPerWin += last_frame_right_ctx_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//	} else {
//		numWrittenFtrPerWin += first_frame_right_ctx_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, stride);
//	}
//
//	if (cur_avail_max_win_len < max_win_len)   //for leading frames
//	{
//		cur_avail_max_win_len++;
//	}
//	else                                       //for all following frames
//	{
//		cur_line += 1;
//		cur_line_ptr += in_width;
//	}
//
//	log.log(QN_LOG_PER_BUNCH, "Read %lu windows for 1 frame.", (unsigned long) multi_win_count);
//
//	return multi_win_count;
}

size_t CRF_InFtrStream_SeqMultiWindow::read_ftrs(size_t cnt, float* ftrs, size_t stride)
{
//	size_t ret = QN_SIZET_BAD;
//	if (stride==num_ftrs() || stride==QN_ALL)
//	{
//		ret = read_ftrs(cnt, ftrs);
//	}
//	else
//	{
//		log.error("The strided version of the read_ftrs function currently does not work for CRF_InFtrStream_SeqMultiWindow.");
//	}
//	return QN_SIZET_BAD;

	const size_t l_out_width = out_width; // Local copy of variable.
	size_t real_stride;         // Stride we use for output.

	if (stride==QN_ALL)
		real_stride = l_out_width;
	else
		real_stride = stride;

//	if (cnt != 1)
//		log.error("Batch size has to be 1: every time read_ftrs() can just move forward by 1 frame.");
	if (cnt > cur_avail_max_win_len)
		log.error("The number of windows(%lu) is set too large, which has to be %lu at current frame.", cnt, cur_avail_max_win_len);
	if (cnt < cur_avail_max_win_len)
		log.error("The number of windows(%lu) is set too small, which has to be %lu at current frame.", cnt, cur_avail_max_win_len);

	if (segno==-1)
		log.error("Trying to read before start of first sentence.");

	if (cur_avail_max_win_len > seg_max_win_len)
		return 0;

	if (cur_avail_max_win_len == max_win_len)
	{
		const size_t l_in_width = in_width; // Local copy of class variable.

		// the max window moves forward by one frame

		// If at end of buffer, read in more lines.
		// modified for context features
//			if (cur_line + max_win_len + bot_margin > buf_lines)
		if (cur_line + max_win_len + right_context_len + bot_margin > buf_lines)
		{
			// Move the remaining lines to the beginning of the buffer.
			size_t line;        // Iterator over lines being moved.

			// modified for context features
//				size_t old_lines = max_win_len + bot_margin - 1; // Lines remaining.
			size_t old_lines = left_context_len + max_win_len +
					  right_context_len + bot_margin - 1; // Lines remaining.

			// modified for context features
//				float *from = cur_line_ptr; // Where we get the line from.
			float *from = cur_line_ptr - left_context_len * in_width; // Where we get the line from.

			float *to = &max_win_buf[0]; // Where we put the line.
			for (line=0; line<old_lines; line++)
			{
				qn_copy_vf_vf(l_in_width, from, to);
				from += l_in_width;
				to += l_in_width;
			}

			//modified for context features
//				cur_line = 0;
//				cur_line_ptr = &max_win_buf[0];
			cur_line = left_context_len;
			cur_line_ptr = &max_win_buf[left_context_len * in_width];

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

//	size_t multi_win_count = cur_avail_max_win_len;
	size_t multi_win_count = cnt;

	size_t numWrittenFtrPerWin = 0;   // the total number of features that have been written to each output window.

	if (this->use_boundary_delta_ftrs)
	{
		numWrittenFtrPerWin += boundary_delta_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
	}
	else
	{
		// added for context features
		// left context features
		numWrittenFtrPerWin += first_frame_left_ctx_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);

		if (max_win_len == 1)
		{
			// modified for context features
	//		numWrittenFtrPerWin += avg_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count);
			numWrittenFtrPerWin += first_frame_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
		} else {
			// modified for context features
	//		numWrittenFtrPerWin += sample_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
	//		numWrittenFtrPerWin += avg_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
	//		numWrittenFtrPerWin += max_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
	//		numWrittenFtrPerWin += min_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
	//		numWrittenFtrPerWin += dur_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
			if (extract_segment_features)
			{
				numWrittenFtrPerWin += sample_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
				numWrittenFtrPerWin += avg_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
				numWrittenFtrPerWin += max_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
				numWrittenFtrPerWin += min_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
				numWrittenFtrPerWin += dur_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
			} else {
				numWrittenFtrPerWin += first_frame_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
			}
		}

		// added for context features
		// right context features
		if (extract_segment_features)
		{
			numWrittenFtrPerWin += last_frame_right_ctx_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
		} else {
			numWrittenFtrPerWin += first_frame_right_ctx_ftrs(ftrs + numWrittenFtrPerWin, multi_win_count, real_stride);
		}
	}

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

/*
 *  CRF_InFtrStream_SeqMultiWindow::sample_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::sample_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		//cout << "cur_win_len=" << cur_win_len << " frameSteps:";
		float* cur_sample_frame_in_buf = cur_win_in_buf;
		float* cur_sample_frame_out_buf = cur_win_out_buf;
		float one_tenth_win_len = cur_win_len * 0.1;
		for (int i = 1; i < 10; i += 2)  //sample frames at 10%, 30%, 50%, 70% and 90% points
		{
			QNUInt32 frameStep = (QNUInt32)ceil(one_tenth_win_len * i) - 1;   //it might be not accurate for some integral value since it is floating point calculation
			//QNUInt32 frameStep = (QNUInt32)(one_tenth_win_len * i);   //it might be not accurate for some integral value since it is floating point calculation

			// just for debugging
//			cout << "frameStep=" << frameStep << endl;

			cur_sample_frame_in_buf += in_width * frameStep;

			// just for debugging
//			cout << "cur_sample_frame_in_buf+=" << in_width << "*" << frameStep << endl;
//			cout << "cur_sample_frame_in_buf[0]=" << cur_sample_frame_in_buf[0];
//			cout << " cur_sample_frame_in_buf[1]=" << cur_sample_frame_in_buf[1];
//			cout << " cur_sample_frame_in_buf[2]=" << cur_sample_frame_in_buf[2] << endl;

			if (cur_sample_frame_out_buf != NULL)
			{
				qn_copy_vf_vf(in_width, cur_sample_frame_in_buf, cur_sample_frame_out_buf);

				// just for debugging
//				for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
//				{
//					cout << "ftr_buf(" << &cur_sample_frame_out_buf[ftr_i] << ")=" << cur_sample_frame_out_buf[ftr_i] << ", sample_ftrs" << endl;
//				}

				cur_sample_frame_out_buf += in_width;
			}
			cur_sample_frame_in_buf = cur_win_in_buf;
			//cout << " " << frameStep;
		}
		cur_win_in_buf -= in_width;
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
		//cout << endl;
	}

	return 5 * in_width;	    //TODO: Currently hard-coded. Need to be parameterized.
}

/*
 *  CRF_InFtrStream_SeqMultiWindow::avg_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::avg_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	float* acc_sum_ftrs(new float[in_width]);
	for (size_t acc_i = 0; acc_i < in_width; acc_i++)
	{
		acc_sum_ftrs[acc_i] = 0.0f;
	}

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
		{
			acc_sum_ftrs[ftr_i] += cur_win_in_buf[ftr_i];
			cur_win_out_buf[ftr_i] = acc_sum_ftrs[ftr_i] / cur_win_len;

			// just for debugging
//			cout << "ftr_buf(" << &cur_win_out_buf[ftr_i] << ")=" << cur_win_out_buf[ftr_i] << ", avg_ftrs" << endl;
		}
		cur_win_in_buf -= in_width;
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	delete [] acc_sum_ftrs;
	return in_width;        //TODO: Currently hard-coded. Need to be parameterized.
}

/*
 *  CRF_InFtrStream_SeqMultiWindow::max_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::max_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	float* acc_max_ftrs(new float[in_width]);
	for (size_t acc_i = 0; acc_i < in_width; acc_i++)
	{
		acc_max_ftrs[acc_i] = cur_win_in_buf[acc_i];
	}

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
		{
			if (cur_win_in_buf[ftr_i] > acc_max_ftrs[ftr_i])
			{
				acc_max_ftrs[ftr_i] = cur_win_in_buf[ftr_i];
			}
			cur_win_out_buf[ftr_i] = acc_max_ftrs[ftr_i];

			// just for debugging
//			cout << "ftr_buf(" << &cur_win_out_buf[ftr_i] << ")=" << cur_win_out_buf[ftr_i] << ", max_ftrs" << endl;
		}
		cur_win_in_buf -= in_width;
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	delete [] acc_max_ftrs;
	return in_width;       //TODO: Currently hard-coded. Need to be parameterized.
}

/*
 *  CRF_InFtrStream_SeqMultiWindow::min_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::min_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	float* acc_min_ftrs(new float[in_width]);
	for (size_t acc_i = 0; acc_i < in_width; acc_i++)
	{
		acc_min_ftrs[acc_i] = cur_win_in_buf[acc_i];
	}

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
		{
			if (cur_win_in_buf[ftr_i] < acc_min_ftrs[ftr_i])
			{
				acc_min_ftrs[ftr_i] = cur_win_in_buf[ftr_i];
			}
			cur_win_out_buf[ftr_i] = acc_min_ftrs[ftr_i];

			// just for debugging
//			cout << "ftr_buf(" << &cur_win_out_buf[ftr_i] << ")=" << cur_win_out_buf[ftr_i] << ", min_ftrs" << endl;
		}
		cur_win_in_buf -= in_width;
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	delete [] acc_min_ftrs;
	return in_width;      //TODO: Currently hard-coded. Need to be parameterized.
}

/*
 *  CRF_InFtrStream_SeqMultiWindow::dur_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::dur_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_out_buf = out_ftr_buf;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (QNUInt32 tmp_win_len=1; tmp_win_len <= max_win_len; tmp_win_len++)
		{
			if (tmp_win_len == cur_win_len)
			{
				cur_win_out_buf[tmp_win_len-1] = 1;

				// just for debugging
//				cout << "ftr_buf(" << &cur_win_out_buf[tmp_win_len-1] << ")=" << cur_win_out_buf[tmp_win_len-1] << ", dur_ftrs" << endl;
			}
			else
			{
				cur_win_out_buf[tmp_win_len-1] = 0;

				// just for debugging
//				cout << "ftr_buf(" << &cur_win_out_buf[tmp_win_len-1] << ")=" << cur_win_out_buf[tmp_win_len-1] << ", dur_ftrs" << endl;
			}
		}
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	return max_win_len;   //TODO: Currently hard-coded. Need to be parameterized.
}


// added for context features
/*
 *  CRF_InFtrStream_SeqMultiWindow::first_frame_left_ctx_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::first_frame_left_ctx_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_first_frame_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		float* cur_left_ctx_frame_in_buf = cur_win_first_frame_in_buf - left_context_len * in_width;
		float* cur_left_ctx_frame_out_buf = cur_win_out_buf;
		for (size_t l_ctx_i = 0; l_ctx_i < left_context_len; l_ctx_i++)
		{
			for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
			{
				cur_left_ctx_frame_out_buf[ftr_i] = cur_left_ctx_frame_in_buf[ftr_i];

				// just for debugging
//				cout << "ftr_buf(" << &cur_left_ctx_frame_out_buf[ftr_i] << ")=" << cur_left_ctx_frame_out_buf[ftr_i] << ", first_frame_left_ctx_ftrs" << endl;
			}
			cur_left_ctx_frame_in_buf += in_width;
			cur_left_ctx_frame_out_buf += in_width;
		}
		cur_win_first_frame_in_buf -= in_width;
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	return in_width * left_context_len;        //TODO: Currently hard-coded. Need to be parameterized.
}

// added for context features
/*
 *  CRF_InFtrStream_SeqMultiWindow::first_frame_right_ctx_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::first_frame_right_ctx_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_first_frame_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		float* cur_right_ctx_frame_in_buf = cur_win_first_frame_in_buf + in_width;
		float* cur_right_ctx_frame_out_buf = cur_win_out_buf;
		for (size_t r_ctx_i = 0; r_ctx_i < right_context_len; r_ctx_i++)
		{
			for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
			{
				cur_right_ctx_frame_out_buf[ftr_i] = cur_right_ctx_frame_in_buf[ftr_i];

				// just for debugging
//				cout << "ftr_buf(" << &cur_right_ctx_frame_out_buf[ftr_i] << ")=" << cur_right_ctx_frame_out_buf[ftr_i] << ", first_frame_right_ctx_ftrs" << endl;
			}
			cur_right_ctx_frame_in_buf += in_width;
			cur_right_ctx_frame_out_buf += in_width;
		}
		cur_win_first_frame_in_buf -= in_width;
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	return in_width * right_context_len;        //TODO: Currently hard-coded. Need to be parameterized.
}

// added for context features
/*
 *  CRF_InFtrStream_SeqMultiWindow::last_frame_right_ctx_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::last_frame_right_ctx_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_last_frame_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		float* cur_right_ctx_frame_in_buf = cur_win_last_frame_in_buf + in_width;
		float* cur_right_ctx_frame_out_buf = cur_win_out_buf;
		for (size_t r_ctx_i = 0; r_ctx_i < right_context_len; r_ctx_i++)
		{
			for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
			{
				cur_right_ctx_frame_out_buf[ftr_i] = cur_right_ctx_frame_in_buf[ftr_i];

				// just for debugging
//				cout << "ftr_buf(" << &cur_right_ctx_frame_out_buf[ftr_i] << ")=" << cur_right_ctx_frame_out_buf[ftr_i] << ", last_frame_right_ctx_ftrs" << endl;
			}
			cur_right_ctx_frame_in_buf += in_width;
			cur_right_ctx_frame_out_buf += in_width;
		}
		// the last frame is the same for all windows, so cur_win_last_frame_in_buf keeps the same.
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	return in_width * right_context_len;        //TODO: Currently hard-coded. Need to be parameterized.
}

// added for context features
/*
 *  CRF_InFtrStream_SeqMultiWindow::first_frame_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::first_frame_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	//for win_len=1
	float* cur_win_first_frame_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
		{
			cur_win_out_buf[ftr_i] = cur_win_first_frame_in_buf[ftr_i];

			// just for debugging
//			cout << "ftr_buf(" << &cur_win_out_buf[ftr_i] << ")=" << cur_win_out_buf[ftr_i] << ", first_frame_ftrs" << endl;
		}
		cur_win_first_frame_in_buf -= in_width;
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	return in_width;        //TODO: Currently hard-coded. Need to be parameterized.
}

// added for context features
/*
 *  CRF_InFtrStream_SeqMultiWindow::boundary_delta_ftrs
 *
 *  Input: out_multi_win_buf: output multiple windows feature buffer
 *         avail_max_win_len: available maximum window length
 *
 *  Return: the number of features that have been written in each output window
 *
 */
size_t CRF_InFtrStream_SeqMultiWindow::boundary_delta_ftrs(float* out_ftr_buf, size_t avail_max_win_len, size_t stride)
{
	if (this->extract_segment_features)
	{
		log.error("extract_segment_features must be false to use boundary_delta_ftrs.");
	}

	QNUInt32 both_context_len = left_context_len;
	if (both_context_len > right_context_len + 1)
	{
		both_context_len = right_context_len + 1;
	}

	//for win_len=1
	float* cur_win_first_frame_in_buf = cur_line_ptr + (avail_max_win_len - 1) * in_width;
	float* cur_win_out_buf = out_ftr_buf;

	for (QNUInt32 cur_win_len=1; cur_win_len <= avail_max_win_len; cur_win_len++)
	{
		float* cur_left_side_frame_in_buf = cur_win_first_frame_in_buf - in_width;
		float* cur_right_side_frame_in_buf = cur_win_first_frame_in_buf;
		float* cur_delta_frame_out_buf = cur_win_out_buf;
		for (size_t side_frame_idx = 0; side_frame_idx < both_context_len; side_frame_idx++)
		{
			// just for debugging
//			cout << "Context " << side_frame_idx << endl;

			for (size_t ftr_i = 0; ftr_i < in_width; ftr_i++)
			{
				float deltaFtr;
				if (cur_left_side_frame_in_buf[ftr_i] >= cur_right_side_frame_in_buf[ftr_i])
				{
					deltaFtr = cur_left_side_frame_in_buf[ftr_i] - cur_right_side_frame_in_buf[ftr_i];
				}
				else
				{
					deltaFtr = cur_right_side_frame_in_buf[ftr_i] - cur_left_side_frame_in_buf[ftr_i];
				}
				cur_delta_frame_out_buf[ftr_i] = deltaFtr;

				// just for debugging
//				cout << "ftr_buf(" << &cur_delta_frame_out_buf[ftr_i] << ")=" << deltaFtr << ", left frame=" << cur_left_side_frame_in_buf[ftr_i] << ", right frame=" << cur_right_side_frame_in_buf[ftr_i] << ", boundary_delta_ftrs" << endl;
			}
			cur_left_side_frame_in_buf -= in_width;
			cur_right_side_frame_in_buf += in_width;
			cur_delta_frame_out_buf += in_width;
		}
		cur_win_first_frame_in_buf -= in_width;
//		cur_win_out_buf += out_width;
		cur_win_out_buf += stride;
	}

	return in_width * both_context_len;        //TODO: Currently hard-coded. Need to be parameterized.
}
