/*
 * CRF_InFtrStream_SeqMultiWindow.h
 *
 *  Created on: Aug 19, 2011
 *      Author: hey
 */

#ifndef CRF_INFTRSTREAM_SEQMULTIWINDOW_H_
#define CRF_INFTRSTREAM_SEQMULTIWINDOW_H_

#include "../CRF.h"

class CRF_InFtrStream_SeqMultiWindow : public QN_InFtrStream {

protected:
	// Logging object.
    QN_ClassLogger log;
    // The input feature stream we are windowing.
    QN_InFtrStream& in_str;
    // The number of frames in the longest window.
    const size_t max_win_len;
    // The number of frames in the longest window for the current sentence.
    // seg_max_win_len < max_win_len if the sentence is shorter than top_margin + max_win_len + bot_margin;
    // Otherwise, seg_max_win_len = max_win_len
    size_t seg_max_win_len;
    // The number of frames in the longest currently available window(cur_avail_max_win_len <= max_win_len).
    // cur_avail_max_win_len increases from 1 to max_win_len for the leading max_win_len frames;
    // cur_avail_max_win_len = max_win_len for all other following frames.
    size_t cur_avail_max_win_len;
    // The number of frames we discard at the start of each segment.
    const size_t top_margin;
    // The number of frames we discard at the end of each segment.
    const size_t bot_margin;
    // The number of features in one frame of the input stream.
    const size_t in_width;
    // The number of features in the longest input window
    const size_t max_win_in_width;
    // The number of features in one output window full.
    size_t out_width;
    // The number of frames we read from the input buffer in one bunch.
    const size_t bunch_frames;
    // The size of the longest input window buffer(max_win_buf) in frames.
    const size_t max_buf_lines;
    // The buffer used for storing the input features of the longest input window.
    float* max_win_buf;
    // The number of lines in max_win_buf.
    size_t buf_lines;
    // The current line (frame) we read from in max_win_buf;
    size_t cur_line;
    // Pointer to start of current line in max_win_buf.
    float* cur_line_ptr;
    // Pointer to first line in long_win_buf after top margin.
    float* first_line_ptr;

    // The number of the current segment (for debugging).
    long segno;

    size_t sample_ftrs(float* out_ftr_buf, size_t avail_max_win_len);
    size_t avg_ftrs(float* out_ftr_buf, size_t avail_max_win_len);
    size_t max_ftrs(float* out_ftr_buf, size_t avail_max_win_len);
    size_t min_ftrs(float* out_ftr_buf, size_t avail_max_win_len);
    size_t dur_ftrs(float* out_ftr_buf, size_t avail_max_win_len);

public:
	CRF_InFtrStream_SeqMultiWindow(int a_debug, const char* a_dbgname,
            QN_InFtrStream& a_str,
            size_t a_max_win_len, size_t a_top_margin,
            size_t a_bot_margin,
            size_t a_bunch_frames = QN_SIZET_BAD);
	virtual ~CRF_InFtrStream_SeqMultiWindow();
    size_t num_ftrs();
    QN_SegID nextseg();
    size_t read_ftrs(size_t cnt, float* ftrs);
    size_t read_ftrs(size_t cnt, float* ftrs, size_t stride);
    int rewind();
    size_t num_segs();
    size_t num_frames(size_t segno = QN_ALL);
    int get_pos(size_t* segno, size_t* frameno);
    QN_SegID set_pos(size_t segno, size_t frameno);

};

#endif /* CRF_INFTRSTREAM_SEQMULTIWINDOW_H_ */
