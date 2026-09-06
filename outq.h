/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OUTQ_H_
#define OUTQ_H_

#include "assert_helpers.h"
#include "ds.h"
#include "sstring.h"
#include "read.h"
#include "threading.h"
#include "mem_ids.h"
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "concurrentqueue.h"

/**
 * Encapsulates a list of lines of output.  If the earliest as-yet-unreported
 * read has id N and Bowtie 2 wants to write a record for read with id N+1, we
 * resize the lines_ and committed_ lists to have at least 2 elements (1 for N,
 * 1 for N+1) and return the BTString * associated with the 2nd element.  When
 * the user calls commit() for the read with id N, 
 */
class OutputQueue {

	static const size_t NFLUSH_THRESH = 8;

public:

	OutputQueue(
		OutFileBuf& obuf,
		bool reorder,
		size_t nthreads,
		bool threadSafe,
		TReadId rdid = 0) :
		obuf_(obuf),
		cur_(rdid),
		nstarted_(0),
		nfinished_(0),
		nflushed_(0),
		lines_(RES_CAT),
		started_(RES_CAT),
		finished_(RES_CAT),
		reorder_(reorder),
		threadSafe_(threadSafe),
        mutex_m()
	{
		assert(nthreads <= 1 || threadSafe);
		nstarted_2.store(0, std::memory_order_release);
		nfinished_2.store(0, std::memory_order_release);
		nflushed_2.store(0, std::memory_order_release);
		output_work.store(true, std::memory_order_release);
	}

	/**
	 * Signal that no more output will be produced and ask the output thread
	 * to drain any remaining records and exit.  The main thread has already
	 * flushed all pending records (oq.flush(true)) before calling this, so by
	 * the time output_thread.join() returns the output thread has drained
	 * the queue and written everything to obuf_.
	 */
	void endoutput()
	{
		output_work.store(false, std::memory_order_release);
		{
			std::lock_guard<std::mutex> lk(cv_m);
			cv.notify_one();
		}
	}

	/**
	 * Output thread entry point.  Drains the lock-free concurrent queue and
	 * writes each record to obuf_ without holding cv_m.  When the queue is
	 * empty and output is still active it parks on the condition variable;
	 * the producer notifies (holding cv_m) after each enqueue so no wakeup is
	 * lost.  When output_work goes false and the queue is drained, it exits.
	 */
	void get_output_from_queue_2()
	{
		BTString temp;
		while(true) {
			// Drain whatever is available (lock-free queue needs no lock, and
			// obuf_ writes must not hold cv_m).
			while(output_queue_2.try_dequeue(temp)) {
				obuf_.writeString(temp);
			}
			std::unique_lock<std::mutex> lk(cv_m);
			// Park until either more records arrive or output is being shut down.
			// The predicate is re-checked under cv_m, and the producer notifies
			// while holding cv_m, so this cannot miss a wakeup.
			cv.wait(lk, [this]() {
				return !output_work.load(std::memory_order_acquire) ||
				       output_queue_2.size_approx() != 0;
			});
			if(!output_work.load(std::memory_order_acquire) &&
			   output_queue_2.size_approx() == 0) {
				break;
			}
		}
	}

	/**
	 * Caller is telling us that they're about to write output record(s) for
	 * the read with the given id.
	 */
	void beginRead(TReadId rdid, size_t threadId);
	
	/**
	 * Writer is finished writing to 
	 */
	void finishRead(const BTString& rec, TReadId rdid, size_t threadId);
	
	/**
	 * Return the number of records currently being buffered.
	 */
	size_t size() const {
		return lines_.size();
	}
	
	/**
	 * Return the number of records that have been flushed so far.
	 */
	TReadId numFlushed() const {
		return nflushed_2.load(std::memory_order_acquire);
	}

	/**
	 * Return the number of records that have been started so far.
	 */
	TReadId numStarted() const {
		return nstarted_2.load(std::memory_order_acquire);
	}

	/**
	 * Return the number of records that have been finished so far.
	 */
	TReadId numFinished() const {
		return nfinished_2.load(std::memory_order_acquire);
	}

	/**
	 * Write already-committed lines starting from cur_.
	 */
	void flush(bool force = false, bool getLock = true);

protected:

	OutFileBuf&     obuf_;
	TReadId         cur_;
	TReadId         nstarted_;
	TReadId         nfinished_;
	TReadId         nflushed_;
	EList<BTString> lines_;
	EList<bool>     started_;
	EList<bool>     finished_;
	bool            reorder_;
	bool            threadSafe_;
	MUTEX_T         mutex_m;
	moodycamel::ConcurrentQueue<BTString> output_queue_2;
	std::atomic<TReadId> nstarted_2;
	std::atomic<TReadId> nfinished_2;
	std::atomic<TReadId> nflushed_2;
	std::atomic<bool> output_work;
	// Protects only the condition-variable wait/signal; the data itself is
	// carried by the lock-free output_queue_2.  cv_m is always held by the
	// producer at the moment it notifies, so the consumer's predicate re-check
	// cannot miss a wakeup.
	std::mutex cv_m;
	std::condition_variable cv;
};

class OutputQueueMark {
public:
	OutputQueueMark(
		OutputQueue& q,
		const BTString& rec,
		TReadId rdid,
		size_t threadId) :
		q_(q),
		rec_(rec),
		rdid_(rdid),
		threadId_(threadId)
	{
		q_.beginRead(rdid, threadId);
	}
	
	~OutputQueueMark() {
		q_.finishRead(rec_, rdid_, threadId_);
	}
	
protected:
	OutputQueue& q_;
	const BTString& rec_;
	TReadId rdid_;
	size_t threadId_;
};

#endif
