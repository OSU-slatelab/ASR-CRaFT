/*
 * CRF_InLabStream_SeqMultiWindow.h
 *
 *  Created on: Sep 10, 2011
 *      Author: Yanzhang (Ryan) He
 */

#ifndef CRF_INLABSTREAM_SEQMULTIWINDOW_H_
#define CRF_INLABSTREAM_SEQMULTIWINDOW_H_

#include "../CRF.h"
#include "CRF_FeatureStream.h"   //for the definition of CRF_LAB_BAD

class CRF_InLabStream_SeqMultiWindow : public QN_InLabStream
{

protected:
	// Logging object.
    QN_ClassLogger log;
    // The input feature stream we are windowing.
    QN_InLabStream& in_str;
    // The number of frames in the longest window.
    const size_t max_win_len;
	// The number of frames we discard at the start of each segment.
	const size_t top_margin;
	// The number of frames we discard at the end of each segment.
	const size_t bot_margin;
	// The maximum number of frames in the input buffer in_buf.
	const size_t max_in_buf_size;
	// The maximum number of windows in win_lab and win_t_end.
	const size_t max_out_buf_size;
	// The buffer used for storing the input frame-based label.
	QNUInt32* in_buf;
	// The buffer to store output window-based labels for the current sentence.
	QNUInt32* out_lab_buf;
	// The buffer to store the start time-stamp of each output window.
	QNUInt32* out_t_start_buf;
	// The buffer to store the end time-stamp of each output window.
	QNUInt32* out_t_end_buf;
	// The buffer to store the flag of broken labels to indicate if each label is broken or complete.
	QNUInt32* out_broken_flag_buf;
	// The number of windows in win_lab and win_t_end.
	size_t num_out_wins;
	// The current window in win_lab and win_t_end.
	size_t cur_out_win;
	// The total number of frames in current sentence.
	size_t num_in_frames;
	// The current frame that read_labs() is going to read from.
	size_t cur_in_frame;
	// The exclusive end frame of the sentence excluding bottom margin.
	size_t end_in_frame;

	// The number of the current segment (for debugging).
	long segno;

	size_t groupLabels(size_t lastWinStartFrame, size_t nextWinStartFrame, QNUInt32 lastWinLab);

public:
	CRF_InLabStream_SeqMultiWindow(int a_debug, const char* a_dbgname,
			QN_InLabStream& a_str,
            size_t a_max_win_len, size_t a_top_margin,
            size_t a_bot_margin,
            size_t a_max_in_buf_size = 1024,
            size_t a_max_out_buf_size = 65536);
	virtual ~CRF_InLabStream_SeqMultiWindow();
	size_t num_labs();
	QN_SegID nextseg();
	size_t read_labs(size_t cnt, QNUInt32* labs);
	size_t read_labs(size_t cnt, QNUInt32* labs, size_t stride);
	int rewind();
	size_t num_segs();
	size_t num_frames(size_t segno = QN_ALL);
	int get_pos(size_t* segno, size_t* frameno);
	QN_SegID set_pos(size_t segno, size_t frameno);
};

#endif /* CRF_INLABSTREAM_SEQMULTIWINDOW_H_ */
