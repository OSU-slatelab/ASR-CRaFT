/*
 * CRF_InLabStream_SeqMultiWindow.cpp
 *
 *  Created on: Sep 10, 2011
 *      Author: hey
 */

#include "CRF_InLabStream_SeqMultiWindow.h"

CRF_InLabStream_SeqMultiWindow::CRF_InLabStream_SeqMultiWindow(int a_debug,
		const char* a_dbgname,
		QN_InLabStream& a_str,
        size_t a_max_win_len, size_t a_top_margin,
        size_t a_bot_margin,
        size_t a_max_in_buf_size,
        size_t a_max_out_buf_size)
	: log(a_debug, "CRF_InFtrStream_SeqMultiWindow", a_dbgname),
	  	  in_str(a_str),
	  	  max_win_len(a_max_win_len),
	  	  top_margin(a_top_margin),
	  	  bot_margin(a_bot_margin),
	  	  max_in_buf_size(a_max_in_buf_size),
	  	  max_out_buf_size(a_max_out_buf_size),
	      in_buf(new QNUInt32[max_in_buf_size]),
	      out_lab_buf(new QNUInt32[max_out_buf_size]),
	      out_t_start_buf(new QNUInt32[max_out_buf_size]),
	      out_t_end_buf(new QNUInt32[max_out_buf_size]),
	      out_broken_flag_buf(new QNUInt32[max_out_buf_size]),
	      segno(-1)
{
	//TODO: Currently CRF_InLabStream_SeqMultiWindow only works for the single label case. It needs to be modified to work for the multi-label case.
	if (a_str.num_labs() != 1)
		log.error("Currently CRF_InLabStream_SeqMultiWindow only works when the width of the input label stream = 1.");
	log.log(QN_LOG_PER_RUN, "Extracting window-based labels, max window length=%lu, "
		        "top margin=%lu frames, bottom margin=%lu frames, "
		        "input buffer size=%lu frames, output buffer size=%lu windows",
		        (unsigned long) max_win_len, (unsigned long) top_margin,
		        (unsigned long) bot_margin, (unsigned long) max_in_buf_size,
		        (unsigned long) max_out_buf_size);
}

CRF_InLabStream_SeqMultiWindow::~CRF_InLabStream_SeqMultiWindow()
{
	delete [] in_buf;
	delete [] out_lab_buf;
	delete [] out_t_start_buf;
	delete [] out_t_end_buf;
	delete [] out_broken_flag_buf;
}

size_t CRF_InLabStream_SeqMultiWindow::groupLabels(
		size_t lastWinStartFrame,
		size_t nextWinStartFrame,
		QNUInt32 lastWinLab)
{
	size_t ret_numWin = 0;
	size_t dur = nextWinStartFrame - lastWinStartFrame;
	if (dur <= max_win_len)
	{
		out_lab_buf[num_out_wins] = lastWinLab;
		out_t_start_buf[num_out_wins] = lastWinStartFrame;
		out_t_end_buf[num_out_wins] = nextWinStartFrame - 1;
		out_broken_flag_buf[num_out_wins] = 0;
		num_out_wins++;
		ret_numWin++;
	}
	else    //split the long window
	{
		size_t numPieces;
		if (dur % max_win_len == 0) {
			numPieces = dur / max_win_len;
		} else {
			numPieces = dur / max_win_len + 1;
		}
		assert(numPieces <= dur);
		size_t pieceDur = dur / numPieces;
		size_t remainder = dur % numPieces;
		size_t pieceStart = lastWinStartFrame;
		size_t pieceEnd;
		pieceDur++;
		for (size_t r = 0; r < remainder; r++)
		{
			pieceEnd = pieceStart + pieceDur - 1;
			assert(pieceEnd >= pieceStart);
			out_lab_buf[num_out_wins] = lastWinLab;
			out_t_start_buf[num_out_wins] = pieceStart;
			out_t_end_buf[num_out_wins] = pieceEnd;
			out_broken_flag_buf[num_out_wins] = 1;
			pieceStart = pieceEnd + 1;
			num_out_wins++;
			ret_numWin++;
		}
		pieceDur--;
		for (size_t r = 0; r < numPieces - remainder; r++)
		{
			pieceEnd = pieceStart + pieceDur - 1;
			assert(pieceEnd >= pieceStart);
			out_lab_buf[num_out_wins] = lastWinLab;
			out_t_start_buf[num_out_wins] = pieceStart;
			out_t_end_buf[num_out_wins] = pieceEnd;
			out_broken_flag_buf[num_out_wins] = 1;
			pieceStart = pieceEnd + 1;
			num_out_wins++;
			ret_numWin++;
		}
		assert (pieceStart == nextWinStartFrame);
	}

	return ret_numWin;
}

size_t CRF_InLabStream_SeqMultiWindow::num_labs()
{
	//TODO: Currently hard-coded four kinds of window-based labs: phone, start frame, end frame, broken/complete phone
	return in_str.num_labs() * 4;
}

//Read all frame-level labels in the next sentence and group them into window-level labels, and save them into output buffer for read_labs() to read.
QN_SegID CRF_InLabStream_SeqMultiWindow::nextseg()
{
	QN_SegID segid;             // Returned segment ID.

	segid = in_str.nextseg();

	num_out_wins = 0;
	cur_out_win = 0;
	num_in_frames = 0;
	cur_in_frame = top_margin;

	if (segid!=QN_SEGID_BAD)
	{
		QNUInt32 lastFrameLab = CRF_LAB_BAD;  // QN_UINT32_MAX doesn't work because it is wrongly defined as "0xffffffffu;" with an extra semi-colon.
		size_t winStartFrame = 0;
		size_t totalBufFrames = in_str.read_labs(max_in_buf_size, in_buf);
		num_in_frames += totalBufFrames;
		size_t curBufFrame = top_margin;
		size_t curInStrFrame = top_margin;

		//cout << "top_margin:" << top_margin << ", bot_margin:" << bot_margin << endl;

		while (totalBufFrames > bot_margin)
		{
			size_t availBufEndFrame = totalBufFrames - bot_margin;  // exclusive end
			for (; curBufFrame < availBufEndFrame; curBufFrame++)
			{
				QNUInt32 curFrameLab = in_buf[curBufFrame];
				if (curFrameLab != lastFrameLab)
				{
					if (lastFrameLab != CRF_LAB_BAD)
					{
						groupLabels(winStartFrame, curInStrFrame, lastFrameLab);
					}
					lastFrameLab = curFrameLab;
					winStartFrame = curInStrFrame;
				}
				curInStrFrame++;
			}

			// Move the remaining frames to the beginning of the buffer.
			size_t frame;        // Iterator over frames being moved.
			size_t old_frames = bot_margin; // Frames remaining.
			QNUInt32 *from = &in_buf[curBufFrame]; // Where we get the frames from.
			QNUInt32 *to = &in_buf[0]; // Where we put the frames.
			for (frame = 0; frame < old_frames; frame++)
			{
				*to++ = *from++;
			}
			size_t blankFrames = max_in_buf_size - old_frames;
			size_t new_frames = in_str.read_labs(blankFrames,
					  &in_buf[old_frames]);
			num_in_frames += new_frames;
			totalBufFrames = old_frames + new_frames;

			curBufFrame = 0;
		}

		if (lastFrameLab != CRF_LAB_BAD)
		{
			groupLabels(winStartFrame, curInStrFrame, lastFrameLab);
		}

		segno++;
	}

	if (num_in_frames <= bot_margin)
	{
		end_in_frame = 0;
	}
	else
	{
		end_in_frame = num_in_frames - bot_margin;
	}

	log.log(QN_LOG_PER_SENT, "Moved on to segment %li.", segno);
	return segid;
}

size_t CRF_InLabStream_SeqMultiWindow::read_labs(size_t cnt, QNUInt32* labs)
{
	if (cnt != 1)
		log.error("Batch size has to be 1: every time read_labs() can just move forward by 1 frame.");

	if (segno==-1)
		log.error("Trying to read before start of first sentence.");

	size_t ret_num_lab = 0;

	if (cur_out_win < num_out_wins && cur_in_frame < end_in_frame)
	{
		//just for debugging
//		cout << cur_in_frame << "\t"
//				<< out_lab_buf[cur_out_win] << "\t"
//				<< out_t_start_buf[cur_out_win] << "\t"
//				<< out_t_end_buf[cur_out_win] << "\t"
//				<< out_broken_flag_buf[cur_out_win] << endl;

		size_t ret_lab_idx = 0;
		//for (size_t frameMove = 0; frameMove < cnt; frameMove++)
		//{
			QNUInt32 win_t_end = out_t_end_buf[cur_out_win];
			if (cur_in_frame == win_t_end)
			{
				labs[ret_lab_idx++] = out_lab_buf[cur_out_win];
				labs[ret_lab_idx++] = out_t_start_buf[cur_out_win];
				labs[ret_lab_idx++] = out_t_end_buf[cur_out_win];
				labs[ret_lab_idx++] = out_broken_flag_buf[cur_out_win];
				cur_out_win++;
			}
			else
			{
				labs[ret_lab_idx++] = CRF_LAB_BAD;
				labs[ret_lab_idx++] = CRF_LAB_BAD;
				labs[ret_lab_idx++] = CRF_LAB_BAD;
				labs[ret_lab_idx++] = CRF_LAB_BAD;
			}
			ret_num_lab++;
			cur_in_frame++;
		//}
	}

	return ret_num_lab;
}

size_t CRF_InLabStream_SeqMultiWindow::read_labs(size_t cnt, QNUInt32* labs, size_t stride)
{
	size_t ret = QN_SIZET_BAD;
	if (stride==num_labs() || stride==QN_ALL)
	{
		ret = read_labs(cnt, labs);
	}
	else
	{
		log.error("The strided version of the read_labs function currently does not work for CRF_InLabStream_SeqMultiWindow.");
	}
	return QN_SIZET_BAD;
}

int CRF_InLabStream_SeqMultiWindow::rewind()
{
	int ec;

	ec = in_str.rewind();
	if (ec==QN_OK)
		segno = -1;
	return ec;
}

// Some streams functions that do not work for SeqMultiWindows.

size_t CRF_InLabStream_SeqMultiWindow::num_segs()
{
	return QN_SIZET_BAD;
}

size_t CRF_InLabStream_SeqMultiWindow::num_frames(size_t segno)
{
	return QN_SIZET_BAD;
}

int CRF_InLabStream_SeqMultiWindow::get_pos(size_t* segno, size_t* frameno)
{
	return QN_BAD;
}

QN_SegID CRF_InLabStream_SeqMultiWindow::set_pos(size_t segno, size_t frameno)
{
	return QN_SEGID_BAD;
}
