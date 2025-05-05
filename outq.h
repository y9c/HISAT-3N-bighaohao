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
#include "concurrentqueue.h"	//无锁队列
// 核心
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
		// 创建并启动一个线程来执行 get_output_from_queue_2
		// std::cout<<"start"<<std::endl;
    	// //std::thread output_thread(&OutputQueue::get_output_from_queue_2, this);
		// std::cout<<"end"<<std::endl;
		output_work=true;
	}

	void endoutput()
	{
		size_t temp_size=output_queue_2.size_approx();
		while(true)
		{
			temp_size=output_queue_2.size_approx();
			//std::cout<<"tempsize="<<temp_size<<" ";
			if(temp_size==0)
			{	
				std::this_thread::sleep_for(std::chrono::milliseconds(1)); //睡眠1ms
				if(temp_size==0)
				{output_work=false;break;}
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
		long int value = nfinished_2.load(std::memory_order_acquire);
		return value;
		//return nflushed_;
	}

	/**
	 * Return the number of records that have been started so far.
	 */
	TReadId numStarted() const {
		long int value = nstarted_2.load(std::memory_order_acquire);
		return value;
		//return nstarted_;
	}

	/**
	 * Return the number of records that have been finished so far.
	 */
	TReadId numFinished() const {
		long int value = nfinished_2.load(std::memory_order_acquire);
		return value;
		//return nfinished_;
	}

	/**
	 * Write already-committed lines starting from cur_.
	 */
	void flush(bool force = false, bool getLock = true);

	void get_output_from_queue_2()	//单开一个输出线程，负责从输出缓冲队列到输出缓冲区
	{	
		BTString temp;
		while(output_work)
		{
			while(output_queue_2.try_dequeue(temp))	//成功拿到
			{
				obuf_.writeString(temp);
			}
			//睡眠
			std::this_thread::sleep_for(std::chrono::milliseconds(1)); //睡眠1ms
		}
	}

protected:

	OutFileBuf&     obuf_;	//这是一个缓冲输出流的包装器，能够写入字符和其他数据类型。该类不是同步的；调用者需要负责同步。
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
	moodycamel::ConcurrentQueue<BTString> output_queue_2;	//输出缓冲队列 --》输出到输出缓冲区
	std::atomic<TReadId> nstarted_2;
	std::atomic<TReadId> nfinished_2;
	std::atomic<TReadId> nflushed_2;
	bool output_work;
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
