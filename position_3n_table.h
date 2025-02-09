/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef POSITION_3N_TABLE_H
#define POSITION_3N_TABLE_H

#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <thread>
#include <cassert>
#include "alignment_3n_table.h"
#include <algorithm> // 需要包含头文件
#include <deque> //新增头文件
#include <queue>
#include <omp.h>
#include "readerwriterqueue.h"
#include "concurrentqueue.h"
#include <atomic>

using namespace std;

extern bool CG_only;
extern long long int loadingBlockSize;

/**
 * store unique information for one base information with readID, and the quality.
 */
class uniqueID
{
public:
    unsigned long long readNameID;
    bool isConverted;
    char quality;
    bool removed;

    uniqueID(unsigned long long InReadNameID,
             bool InIsConverted,
             char& InQual){
        readNameID = InReadNameID;
        isConverted = InIsConverted;
        quality = InQual;
        removed = false;
    }
};

/**
 * basic class to store reference position information
 */
class Position{
    mutex mutex_;
public:
    string chromosome; // reference chromosome name
    long long int location; // 1-based position
    char strand; // +(REF) or -(REF-RC)
    string convertedQualities; // each char is a mapping quality on this position for converted base.
    string unconvertedQualities; // each char is a mapping quality on this position for unconverted base.
    bool convertedQualities_flag;
    bool unconvertedQualities_flag;
    vector<uniqueID> uniqueIDs; // each value represent a readName which contributed the base information.
                              // readNameIDs is to make sure no read contribute 2 times in same position.
    //更改为deque以减少emplace操作用时
    //std::deque<uniqueID> uniqueIDs;
    /*Position 类的主要用途是存储参考基因组上某个位置的相关信息，包括：
    染色体名称。
    位置坐标。
    链的方向。
    转换和未转换碱基的映射质量。
    贡献该位置信息的 reads 的唯一标识符。*/


    void initialize() {
        chromosome.clear();
        location = -1;
        strand = '?';
        convertedQualities.clear();
        unconvertedQualities.clear();
        vector<uniqueID>().swap(uniqueIDs);
        //uniqueIDs.reserve(); // 内存预分配 1000 个元素的内存
        //std::deque<uniqueID>().swap(uniqueIDs); // 清空 deque
    }

    Position(){
        initialize();
    };

    /**
     * return true if there is mapping information in this reference position.
     */
    bool empty() {
        return convertedQualities.empty() && unconvertedQualities.empty();
    }

    /**
     * set the chromosome, location (position), and strand information.
     */

    void set (string& inputChr, long long int inputLoc) {
        chromosome = inputChr;
        location = inputLoc + 1;
    }

    //  取消复制,move的时候统一赋值
    void set_without_string(long long int inputLoc)
    {
        location = inputLoc + 1;
    }

    void set(char inputStrand) {
        strand = inputStrand;
    }

    /**
     * binary search of readNameID in readNameIDs.
     * always return a index.
     * if cannot find, return the index which has bigger value than input readNameID.
     */
    //二分法递归查找
    // int searchReadNameID (unsigned long long&readNameID, int start, int end) {
    //     if (uniqueIDs.empty()) {
    //         return 0;
    //     }
    //     if (start <= end) {
    //         int middle = (start + end) / 2;
    //         if (uniqueIDs[middle].readNameID == readNameID) {
    //             return middle;
    //         }
    //         if (uniqueIDs[middle].readNameID > readNameID) {
    //             return searchReadNameID(readNameID, start, middle-1);
    //         }
    //         return searchReadNameID(readNameID, middle+1, end);
    //     }
    //     return start; // return the bigger one
    // }
    //STL优化
    int searchReadNameID(unsigned long long& readNameID, int start, int end) {
        if (uniqueIDs.empty()) {
            return 0;
        }
        int middle=0;
        while (start <= end) {
            middle = (start + end) / 2;
            auto& midValue = uniqueIDs[middle].readNameID;

            if (midValue == readNameID) {
                return middle;
            }
            if (midValue > readNameID) {
                end = middle - 1;
            } else {
                start = middle + 1;
            }
        }
        return start;  // 返回大于等于 readNameID 的最小索引
    }

    /**
     * with a input readNameID, add it into readNameIDs.
     * if the input readNameID already exist in readNameIDs, return false.
     */
    //
    bool appendReadNameID(PosQuality& InBase, Alignment& InAlignment) {
        int idCount = uniqueIDs.size();
        if (idCount == 0 || InAlignment.readNameID > uniqueIDs.back().readNameID) {
            uniqueIDs.emplace_back(InAlignment.readNameID, InBase.converted, InBase.qual);
            return true;
        }
        int index = searchReadNameID(InAlignment.readNameID, 0, idCount);
        if (uniqueIDs[index].readNameID == InAlignment.readNameID) {
            // if the new base is consistent with exist base's conversion status, ignore
            // otherwise, delete the exist conversion status
            if (uniqueIDs[index].removed) {
                return false;
            }
            if (uniqueIDs[index].isConverted != InBase.converted) {
                uniqueIDs[index].removed = true;
                if (uniqueIDs[index].isConverted) {
                    for (int i = 0; i < convertedQualities.size(); i++) {
                        if (convertedQualities[i] == InBase.qual) {
                            convertedQualities.erase(convertedQualities.begin()+i);
                            return false;
                        }
                    }
                } else {
                    for (int i = 0; i < unconvertedQualities.size(); i++) {
                        if (unconvertedQualities[i] == InBase.qual) {
                            unconvertedQualities.erase(unconvertedQualities.begin()+i);
                            return false;
                        }
                    }
                }
            }
            return false;
        } else {
            uniqueIDs.emplace(uniqueIDs.begin()+index, InAlignment.readNameID, InBase.converted, InBase.qual);
            //std::cout<<uniqueIDs.size()<<std::endl;
            return true;
        }    //此处性能开销极大，vector插入操作
    }

    /**
     * append the SAM information into this position.
     */
    //用时巨大函数
    void appendBase (PosQuality& input, Alignment& a) {
        mutex_.lock();      //用时大，一个类只有一个锁？
        if (appendReadNameID(input,a)) {
            if (input.converted) {
                convertedQualities += input.qual;
            } else {
                unconvertedQualities += input.qual;
            }
        }
        mutex_.unlock();
    }
};

/**
 * store all reference position in this class.
 */
class Positions{
public:
    vector<Position*> refPositions; // the pool of all current reference position.
    string chromosome; // current reference chromosome name.
    long long int location; // current location (position) in reference chromosome.
    char lastBase = 'X'; // the last base of reference line. this is for CG_only mode.
    SafeQueue<string*> linePool; // pool to store unprocessed SAM line.
    SafeQueue<string*> freeLinePool; // pool to store free string pointer for SAM line.
    SafeQueue<Position*> freePositionPool; // pool to store free position pointer for reference position.
    SafeQueue<Position*> outputPositionPool; // pool to store the reference position which is loaded and ready to output.

    //更改输出池为普通队列，取消锁
    std::queue<Position*> freePositionPool_2;
    std::queue<Position*> outputPositionPool_2;
    // 更改为无锁队列
    moodycamel::ReaderWriterQueue<Position*> outputPositionPool_3;
    moodycamel::ConcurrentQueue<string*> linePool_3;
    moodycamel::ConcurrentQueue<string*> freeLinePool_3;

    bool working;
    bool output_thread_working; //输出线程是否工作
    mutex mutex_;
    long long int refCoveredPosition; // this is the last position in reference chromosome we loaded in refPositions.
    ifstream refFile;
    vector<mutex*> workerLock; // one lock for one worker thread.
    int nThreads = 1;
    ChromosomeFilePositions chromosomePos; // store the chromosome name and it's streamPos. To quickly find new chromosome in file.
    bool addedChrName = false;
    bool removedChrName = false;
    // vector<atomic<bool>> worker_finished;    //新声明标志量
    // atomic<bool> wk_f;
    //std::atomic<int> doneConsumers(0); //用于同步

    Positions(string inputRefFileName, int inputNThreads, bool inputAddedChrName, bool inputRemovedChrName) {
        working = true;
        output_thread_working=true;
        nThreads = inputNThreads;
        addedChrName = inputAddedChrName;
        removedChrName = inputRemovedChrName;
        for (int i = 0; i < nThreads; i++) {
            workerLock.push_back(new mutex);
            //worker_finished.push_back(std::atomic<bool>(false));
            //wk_f=true;
        }
        refFile.open(inputRefFileName, ios_base::in);
        LoadChromosomeNamesPos();
        freePositionPool_init();
    }

    ~Positions() {
        for (int i = 0; i < workerLock.size(); i++) {
            delete workerLock[i];
        }
        Position* pos;
        while(freePositionPool.popFront(pos)) {
            delete pos;
        }
        while(!freePositionPool_2.empty()) {
            pos=freePositionPool_2.front();
            freePositionPool_2.pop();
            delete pos;
        }
    }

    void freePositionPool_init()
    {
        std::mutex queueMutex;  // 保护 queue 的 mutex
        int positionsPerThread=100000000;
        const int temp_numThreads = 2;  // 启动 2 个线程
        #pragma omp parallel num_threads(temp_numThreads)
        {
            // 每个线程有一个局部的 Position 队列
            std::queue<Position*> localQueue;
            // 获取线程 ID
            int threadId = omp_get_thread_num();
            std::cout<<threadId<<std::endl;

            for (int i = 0; i < positionsPerThread; ++i) {
                Position* newPosition = new Position();  // 创建新对象
                localQueue.push(newPosition);  // 将指针推入该线程的局部队列
            }

            std::cout<<"begin push"<<std::endl;
            // // 合并局部队列到全局队列
            queueMutex.lock();
            while (!localQueue.empty()) {
                freePositionPool_2.push(localQueue.front());
                localQueue.pop();
            }
            queueMutex.unlock();

            std::cout<<"end push"<<std::endl;
        }
    }



    /**
     * given the target Position output the corresponding position index in refPositions.
     */
    int getIndex(long long int &targetPos) {
        int firstPos = refPositions[0]->location;
        return targetPos - firstPos;
    }

    /**
     * given reference line (start with '>'), extract the chromosome information.
     * this is important when there is space in chromosome name. the SAM information only contain the first word.
     */
    string getChrName(string& inputLine) {
        string name;
        for (int i = 1; i < inputLine.size(); i++)
        {
            char c = inputLine[i];
            if (isspace(c)){
                break;
            }
            name += c;
        }

        if(removedChrName) {
            if(name.find("chr") == 0) {
                name = name.substr(3);
            }
        } else if(addedChrName) {
            if(name.find("chr") != 0) {
                name = string("chr") + name;
            }
        }
        return name;
    }


    /**
     * Scan the reference file. Record each chromosome and its position in file.
     */
    void LoadChromosomeNamesPos() {
        std::cout<<"loadchr"<<std::endl;
        string line;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // this line is chromosome name
                chromosome = getChrName(line);
                streampos currentPos = refFile.tellg();
                chromosomePos.append(chromosome, currentPos);
            }
        }
        chromosomePos.sort();
        chromosome.clear();
    }

    /**
     * get a fasta line (not header), append the bases to positions.
     */
    void appendRefPosition(string& line) {
        Position* newPos;
        // check the base one by one
        char* b;
        // Position* newPos_[line.size()];
        // for(int i = 0; i < line.size(); i++)
        // {
        //     getFreePosition(newPos_[i]);
        // }

        // #pragma omp parallel for num_threads(2)
        // for(int i = 0; i < line.size(); i++)
        // {
        //     newPos_[i]->set(chromosome, location+i);
        // }

        // for(int i = 0; i < line.size(); i++)
        // {
        //     b = &line[i];
        //     if (CG_only) {
        //         if (lastBase == 'C' && *b == 'G') {
        //             refPositions.back()->set('+');
        //             newPos_[i]->set('-');
        //         }
        //     } else {
        //         if (*b == convertFrom) {
        //             newPos_[i]->set('+');
        //         } else if (*b == convertFromComplement) {
        //             newPos_[i]->set('-');
        //         }
        //     }
        //     refPositions.push_back(newPos_[i]);
        //     lastBase = *b;
        // }

        // location += line.size();



        for (int i = 0; i < line.size(); i++) {
            getFreePosition(newPos);
            //newPos->set(chromosome, location+i);
            newPos->set_without_string(location+i);
            b = &line[i];
            if (CG_only) {
                if (lastBase == 'C' && *b == 'G') {
                    refPositions.back()->set('+');
                    newPos->set('-');
                }
            } else {
                if (*b == convertFrom) {
                    newPos->set('+');
                } else if (*b == convertFromComplement) {
                    newPos->set('-');
                }
            }
            refPositions.push_back(newPos);
            lastBase = *b;
        }
        location += line.size();
    }

    /**
     * if we can go through all the workerLock, that means no worker is appending new position.
     */
    void appendingFinished() {
        for (int i = 0; i < nThreads; i++) {
            workerLock[i]->lock();
            workerLock[i]->unlock();
        }
    }

    /**
     * the output function for output thread.
     */
    void outputFunction(string outputFileName) {
        ostream* out_ = &cout;
        out_ = &cout;
        ofstream tableFile;
        if (!outputFileName.empty()) {
            tableFile.open(outputFileName, ios_base::out);
            out_ = &tableFile;
        }

        *out_ << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount\n";
        Position* pos;
        // 原代码
        // while (working) {
        //     if (outputPositionPool.popFront(pos)) {
        //         *out_ << pos->chromosome << '\t'
        //                   << to_string(pos->location) << '\t'
        //                   << pos->strand << '\t'
        //                   << pos->convertedQualities << '\t'
        //                   << to_string(pos->convertedQualities.size()) << '\t'
        //                   << pos->unconvertedQualities << '\t'
        //                   << to_string(pos->unconvertedQualities.size()) << '\n';
        //         returnPosition(pos);
        //     } else {
        //         this_thread::sleep_for (std::chrono::microseconds(1));
        //     }
        // }

        //修正版
        // while (working) {
        //     if (!outputPositionPool_2.empty() && output_thread_working) {
        //         pos=outputPositionPool_2.front();
        //         outputPositionPool_2.pop();
        //         *out_ << pos->chromosome << '\t'
        //                   << to_string(pos->location) << '\t'
        //                   << pos->strand << '\t'
        //                   << pos->convertedQualities << '\t'
        //                   << to_string(pos->convertedQualities.size()) << '\t'
        //                   << pos->unconvertedQualities << '\t'
        //                   << to_string(pos->unconvertedQualities.size()) << '\n';
        //         returnPosition(pos);
        //     } else {
        //         this_thread::sleep_for (std::chrono::microseconds(1));
        //     }
        // }

        // 无锁队列修正版 
        bool succeeded;
        while (working) {
            if (outputPositionPool_3.try_dequeue(pos)) { //尝试从队列中取出一个元素。如果队列为空，try_dequeue 会返回 false。
                //如果取出成功，返回的元素将保存在 number 中。
                *out_ << pos->chromosome << '\t'
                          << to_string(pos->location) << '\t'
                          << pos->strand << '\t'
                          << pos->convertedQualities << '\t'
                          << to_string(pos->convertedQualities.size()) << '\t'
                          << pos->unconvertedQualities << '\t'
                          << to_string(pos->unconvertedQualities.size()) << '\n';
                returnPosition(pos);
            } else {
                this_thread::sleep_for (std::chrono::microseconds(10000));
            }
        }
        tableFile.close();
    }
 
    /**
     * move the position which position smaller than refCoveredPosition - loadingBlockSize, output it.
     */
    void moveBlockToOutput() {
        if (refPositions.empty()) {
            return;
        }
        int index;
        for (index = 0; index < refPositions.size(); index++) {
            if (refPositions[index]->location < refCoveredPosition - loadingBlockSize) {
                if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                    returnPosition(refPositions[index]);
                } else {
                    //outputPositionPool.push(refPositions[index]);
                    //outputPositionPool_2.push(refPositions[index]);
                    outputPositionPool_3.enqueue(refPositions[index]);
                }
            } else {
                break;
            }
        }
        if (index != 0) {
            refPositions.erase(refPositions.begin(), refPositions.begin()+index);
        }
    }

    /**
     * move all the refPosition into output pool.
     */
    void moveAllToOutput() {
        if (refPositions.empty()) {
            return;
        }

        int index=0;
        int tempsize=refPositions.size();
        bool temp_empty_flag[tempsize];
        std::cout<<"refpositions size = "<<tempsize<<std::endl;
        std::mutex queueMutex;
        int count=0;
        #pragma omp parallel num_threads(8)
        {   
            // 局部池
            std::queue<Position*> localFreePool;
            #pragma omp for
            for (index = 0; index < tempsize; index++) {
                temp_empty_flag[index]=  (refPositions[index]->strand == '?' || refPositions[index]->empty());
                if (temp_empty_flag[index]) {
                    refPositions[index]->initialize();
                    localFreePool.push(refPositions[index]);
                } 
                else {
                    vector<uniqueID>().swap(refPositions[index]->uniqueIDs);
                    //outputPositionPool_2.push(refPositions[index]); 有序步骤放到最后
                }
            }
            int thread_id = omp_get_thread_num(); 
            //// 合并局部的 freePositionPool 到全局的 freePositionPool_2
            queueMutex.lock();
            //std::cout<<"now localfreepool size="<<localFreePool.size()<<" "<<thread_id<<std::endl;
            while(!localFreePool.empty())
            {
                freePositionPool_2.push(localFreePool.front());
                localFreePool.pop();
            }
            queueMutex.unlock();
        }

        string temp_chm=chromosome; //暂存
        //下面这个串行
        for (index = 0; index < tempsize; index++) {
            if (!temp_empty_flag[index]) {              //大部分是空，小部分是非空
                //outputPositionPool_2.push(refPositions[index]);
                refPositions[index]->chromosome=temp_chm;
                outputPositionPool_3.enqueue(refPositions[index]);
                count++;
            } 
        }
        std::cout<<"all !temp_empty_flag_count="<<count<<std::endl;


        // //原来的代码
        // for (int index = 0; index < refPositions.size(); index++) {
        //     if (refPositions[index]->strand == '?' || refPositions[index]->empty()) {   //短路求值，更改顺序，前者性能开销更小
        //         returnPosition(refPositions[index]);
        //     } else {
        //         vector<uniqueID>().swap(refPositions[index]->uniqueIDs);
        //         //outputPositionPool.push(refPositions[index]);
        //         outputPositionPool_2.push(refPositions[index]);
        //     }
        // }
        refPositions.clear();
    }

    /**
     * initially load reference sequence for 2 million bp
     */
    void loadNewChromosome(string targetChromosome) {
        refFile.clear();
        // find the start position in file based on chromosome name.
        streampos startPos = chromosomePos.getChromosomePosInRefFile(targetChromosome);
        chromosome = targetChromosome;
        refFile.seekg(startPos, ios::beg);
        refCoveredPosition = 2 * loadingBlockSize;
        string line;
        lastBase = 'X';
        location = 0;
        std::cout<<"using loadnewChr now free pool="<<freePositionPool_2.size()<<std::endl;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // this line is chromosome name
                return; // meet next chromosome, return it.
            } else {
                if (line.empty()) { continue; }
                // change all base to upper case
                for (int i = 0; i < line.size(); i++) {
                    line[i] = toupper(line[i]);
                }
                //std::cout<<line.size()<<std::endl;
                appendRefPosition(line);    //性能瓶颈          line='ACATTGTTGCCAAATATAAATAGTGAGAAAAGCATTTTATATTCCCTAAGGCTCCTTGAC'
                if (location >= refCoveredPosition) {
                    return;
                }
            }
        }
    }

    /**
     * load more Position (loadingBlockSize bp) to positions
     * if we meet next chromosome, return false. Else, return ture.
     */
    void loadMore() {
        refCoveredPosition += loadingBlockSize;
        string line;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // meet next chromosome, return.
                return ;
            } else {
                if (line.empty()) { continue; }

                // change all base to upper case
                for (int i = 0; i < line.size(); i++) {
                    line[i] = toupper(line[i]);
                }

                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return ;
                }
            }
        }
    }


    /**
     * add position information from Alignment into ref position.
     */
    //将 newAlignment（一个 Alignment 对象）的位置信息添加到 refPositions 中
    void appendPositions(Alignment& newAlignment) {
        //检查是否是有效的比对
        if (!newAlignment.mapped || newAlignment.bases.empty()) {
            return;
        }
        //计算比对的起始位置
        long long int startPos = newAlignment.location; // 1-based position
        // find the first reference position in pool.
        int index = getIndex(newAlignment.location);    //获取参考位置索引

        for (int i = 0; i < newAlignment.sequence.size(); i++) {    //遍历读取序列的每个碱基
            PosQuality* b = &newAlignment.bases[i];
            if (b->remove) {
                continue;
            }

            Position* pos = refPositions[index+b->refPos];  //将比对位置添加到参考位置
            //assert (pos->location == startPos + b->refPos); debug模式下保留，release可去除

            if (pos->strand == '?') {
                // this is for CG-only mode. read has a 'C' or 'G' but not 'CG'.
                continue;
            }
            pos->appendBase(newAlignment.bases[i], newAlignment); 
            //调用 Position 类的 appendBase 方法，将碱基信息添加到参考位置,主要耗时
        }
    }

    /**
     * get a string pointer from freeLinePool, if freeLinePool is empty, make a new string pointer.
     */
    void getFreeStringPointer(string*& newLine) {
       // if (freeLinePool.popFront(newLine)) {
        if (freeLinePool_3.try_dequeue(newLine)) {
            return;
        } else {
            newLine = new string();
        }
    }

    /**
     * get a Position pointer from freePositionPool, if freePositionPool is empty, make a new Position pointer.
     */
    void getFreePosition(Position*& newPosition) {          //loadnewch 加载新染色体主线程调用
        while (outputPositionPool_2.size() >= 10000000) {     //原来为outputpositiongpool
            this_thread::sleep_for (std::chrono::microseconds(1));
        }
        // if (freePositionPool.popFront(newPosition)) {
        //     return;
        // } else {
        //     newPosition = new Position();
        // }
        //修改为普通queue
        if (!freePositionPool_2.empty()) {
            newPosition=freePositionPool_2.front();
            freePositionPool_2.pop();
            return;
        } else {
            newPosition = new Position();
        }
    }

    /**
     * return the line to freeLinePool
     */
    void returnLine(string* line) {
        line->clear();
        //freeLinePool.push(line);
        freeLinePool_3.enqueue(line);
    }

    /**
     * return the position to freePositionPool.
     */
    // move all to output 主线程调用
    void returnPosition(Position* pos) {
        pos->initialize();
        //freePositionPool.push(pos);
        freePositionPool_2.push(pos);
    }

    // void set_worker_false()
    // {
    //     for(int i=0; i < nThreads; i++)
    //     {
    //         worker_finished[i]=false;
    //     }
    // }

    // bool get_worker_finished()
    // {
    //     bool flag=true;
    //     for(int i=0; i < nThreads; i++)
    //     {
    //         if(worker_finished[i]=false)
    //         {
    //             flag=false;
    //         }
    //     }
    //     return flag;
    // }
    
    /**
     * this is the working function.
     * it take the SAM line from linePool, parse it.
     */
    //工作线程：从linepool池中抓取数据并处理
    void append(int threadID) {
        string* line;
        Alignment newAlignment;

        string* lines[1000]; //size大小
        // int work_block_size=100;
        // int count=0;
        
        // //源代码
        // while (working) {
        //     workerLock[threadID]->lock();       //无用时
        //     //if(!linePool.popFront(line)) {
        //     if(!linePool_3.try_dequeue(line)) {     //成功返回true ！ -> false
        //         //worker_finished[threadID]=true;
        //         workerLock[threadID]->unlock();         //进入此步则表示已使用完
        //         this_thread::sleep_for (std::chrono::nanoseconds(1000));   //用时4%
        //         //wk_f=true;
        //         continue;
        //     }
        //     while (refPositions.empty()) {
        //         this_thread::sleep_for (std::chrono::microseconds(1000));
        //     }
        //     newAlignment.parse(line);   //用时较少  9%
        //     returnLine(line);   //解析完成返还资源到freepool    %2
        //     appendPositions(newAlignment);  //用时较多 30%
        //     workerLock[threadID]->unlock();
        // }

        // while (working) {
        //     workerLock[threadID]->lock();
        //     //避免竞争，块处理   r  
        //     while(count < work_block_size && !linePool_3.try_dequeue(line))
        //     {
        //         lines[count]=line;
        //         count++;
        //     }
        //     // 如果队列空了或者达到块大小:则下一步
        //     while (refPositions.empty()) {
        //         this_thread::sleep_for (std::chrono::microseconds(1000));
        //     }
        //     for(int i=0;i<count;i++)
        //     {
        //         newAlignment.parse(lines[i]);   //用时较少  9%
        //         returnLine(lines[i]);   //解析完成返还资源到freepool    %2
        //         appendPositions(newAlignment);  //用时较多 30%
        //     }
        //     count=0;
        //     workerLock[threadID]->unlock();
        // }

        std::size_t temp_count=0;
        int bulk_size=1000;
        while (working) {
            workerLock[threadID]->lock();       //无用时
            //从队列中抓取
            temp_count= linePool_3.try_dequeue_bulk(lines, bulk_size);
            //std::cout<<temp_count<<std::endl;
            if(temp_count==0)
            {
                workerLock[threadID]->unlock(); 
                this_thread::sleep_for (std::chrono::nanoseconds(10000));   //用时4%
                continue;
            }
            while (refPositions.empty()) {
                this_thread::sleep_for (std::chrono::microseconds(10000));
            }
            for(int i=0;i<temp_count;i++)
            {
                newAlignment.parse(lines[i]);   //用时较少  9%
                returnLine(lines[i]);   //解析完成返还资源到freepool    %2
            appendPositions(newAlignment);  //用时较多 30%
            }
            workerLock[threadID]->unlock();
        }
    }

};

#endif //POSITION_3N_TABLE_H
