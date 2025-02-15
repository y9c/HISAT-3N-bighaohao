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

#include "outq.h"
#include <thread>  // 用于std::this_thread::sleep_for
#include <chrono>

/**
 * Caller is telling us that they're about to write output record(s) for
 * the read with the given id.
 */
void OutputQueue::beginRead(TReadId rdid, size_t threadId) {
	//ThreadSafe t(&mutex_m, threadSafe_);	//开头上锁
	nstarted_++;
	nstarted_2.fetch_add(1,std::memory_order_relaxed);
	if(reorder_) {
		ThreadSafe t(&mutex_m, threadSafe_);
		assert_geq(rdid, cur_);
		assert_eq(lines_.size(), finished_.size());
		assert_eq(lines_.size(), started_.size());
		if(rdid - cur_ >= lines_.size()) {
			// Make sure there's enough room in lines_, started_ and finished_
			size_t oldsz = lines_.size();
			lines_.resize(rdid - cur_ + 1);
			started_.resize(rdid - cur_ + 1);
			finished_.resize(rdid - cur_ + 1);
			for(size_t i = oldsz; i < lines_.size(); i++) {
				started_[i] = finished_[i] = false;
			}
		}
		started_[rdid - cur_] = true;
		finished_[rdid - cur_] = false;
	}									//t处理时解锁
	if(nstarted_%1000000==0)
	{
		long int value =numStarted();
		long int value2 =numFinished();
		long int value3 =numFlushed();
		printf("OutputQueue::beginRead nstarted=%d nstarted_2=%ld rdid=%d cur_=%d tid=%d nfinished_=%d   nfinished_2=%ld,nflush_2=%ld\n",nstarted_,value,rdid,cur_,threadId,nfinished_,value2,value3);
	}
	//printf("OutputQueue::beginRead nstarted=%d rdid=%d cur_=%d tid=%d\n",nstarted_,rdid,cur_,threadId);
	// 休眠 1 秒钟
    //std::this_thread::sleep_for(std::chrono::milliseconds(10));
}

/**
 * Writer is finished writing to 
 */
void OutputQueue::finishRead(const BTString& rec, TReadId rdid, size_t threadId) {
	//ThreadSafe t(&mutex_m, threadSafe_);
	if(reorder_) {
		ThreadSafe t(&mutex_m, threadSafe_);
		assert_geq(rdid, cur_);
		assert_eq(lines_.size(), finished_.size());
		assert_eq(lines_.size(), started_.size());
		assert_lt(rdid - cur_, lines_.size());
		assert(started_[rdid - cur_]);
		assert(!finished_[rdid - cur_]);
		lines_[rdid - cur_] = rec;
		nfinished_++;
		finished_[rdid - cur_] = true;
		flush(false, false); // don't force; already have lock
	} else {
		// obuf_ is the OutFileBuf for the output file
		//obuf_.writeString(rec);	//写ref --> 写缓存
		output_queue_2.enqueue(rec);	//先写入输出队列中
		nfinished_++;
		nflushed_++;
		nfinished_2.fetch_add(1,std::memory_order_relaxed);
		nflushed_2.fetch_add(1,std::memory_order_relaxed);
	}
}

/**
 * Write already-finished lines starting from cur_.
 */
void OutputQueue::flush(bool force, bool getLock) {
	// 获取当前时间（开始时间）
    auto start = std::chrono::high_resolution_clock::now();
	printf("using flush\n");
	std::cout<<reorder_<<std::endl;
	if(!reorder_) {
		return;
	}
	std::cout<<reorder_<<std::endl;
	ThreadSafe t(&mutex_m, getLock && threadSafe_);
	size_t nflush = 0;
	while(nflush < finished_.size() && finished_[nflush]) {
		assert(started_[nflush]);
		nflush++;
	}
	printf("nflush = %ld\n lines_=  %ld NFLUSH_THRESH=%ld",nflush,lines_.size(),NFLUSH_THRESH);
	// Waiting until we have several in a row to flush cuts down on copies
	// (but requires more buffering)
	//最后一次强制force
	if(force || nflush >= NFLUSH_THRESH) {
		std::cout<<force<<" forcing "<<nflush<<std::endl;
		for(size_t i = 0; i < nflush; i++) {
			assert(started_[i]);
			assert(finished_[i]);
			obuf_.writeString(lines_[i]);
		}
		lines_.erase(0, nflush);
		started_.erase(0, nflush);
		finished_.erase(0, nflush);
		cur_ += nflush;
		nflushed_ += nflush;
	}
	    // 获取当前时间（结束时间）
    auto end = std::chrono::high_resolution_clock::now();
	  // 计算操作的持续时间
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);  // 转换为毫秒
    std::cout << "操作耗时： " << duration.count() << " 毫秒" << std::endl;
}

#ifdef OUTQ_MAIN

#include <iostream>

using namespace std;

int main(void) {
	cerr << "Case 1 (one thread) ... ";
	{
		OutFileBuf ofb;
		OutputQueue oq(ofb, false);
		assert_eq(0, oq.numFlushed());
		assert_eq(0, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(1);
		assert_eq(0, oq.numFlushed());
		assert_eq(1, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(3);
		assert_eq(0, oq.numFlushed());
		assert_eq(2, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(2);
		assert_eq(0, oq.numFlushed());
		assert_eq(3, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(3, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(0);
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.finishRead(0);
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.flush(true);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.finishRead(2);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(2, oq.numFinished());
		oq.flush(true);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(2, oq.numFinished());
		oq.finishRead(1);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(3, oq.numFinished());
		oq.flush(true);
		assert_eq(3, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(3, oq.numFinished());
	}
	cerr << "PASSED" << endl;

	cerr << "Case 2 (one thread) ... ";
	{
		OutFileBuf ofb;
		OutputQueue oq(ofb, false);
		BTString& buf1 = oq.beginRead(0);
		BTString& buf2 = oq.beginRead(1);
		BTString& buf3 = oq.beginRead(2);
		BTString& buf4 = oq.beginRead(3);
		BTString& buf5 = oq.beginRead(4);
		assert_eq(5, oq.numStarted());
		assert_eq(0, oq.numFinished());
		buf1.install("A\n");
		buf2.install("B\n");
		buf3.install("C\n");
		buf4.install("D\n");
		buf5.install("E\n");
		oq.finishRead(4);
		oq.finishRead(1);
		oq.finishRead(0);
		oq.finishRead(2);
		oq.finishRead(3);
		oq.flush(true);
		assert_eq(5, oq.numFlushed());
		assert_eq(5, oq.numStarted());
		assert_eq(5, oq.numFinished());
		ofb.flush();
	}
	cerr << "PASSED" << endl;
	return 0;
}

#endif /*def ALN_SINK_MAIN*/
