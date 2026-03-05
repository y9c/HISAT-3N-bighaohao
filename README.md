# hisat-3n优化实现 （ASC25赛题）
本仓库包含初赛技术报告、初赛可行工作流、以及最新优化(决赛)的hisat-3n工具代码。决赛更多的是调度方面的优化，
[初赛报告链接.pdf](RNA-m5c-report.pdf)
[初赛可行工作流.tar.gz](RNA.tar.gz)
[初赛赛题](ASC25_Preliminary_Round_Announcement.pdf)
下面是初赛优化总结:
---
# RNA m5C Site Detection Workflow Acceleration (HISAT-3N / HISAT-3N-Table)
本仓库提供一套**可复现的 RNA m5C 修饰位点检测工作流**，并在保持检测流程正确性/可靠性的前提下，对工作流及关键工具进行了系统级性能优化。我们结合**软件插桩（instrumentation）**与 **Intel VTune Profiler** 定位瓶颈，围绕工作流中最耗时的两个组件 **`hisat-3n-table`** 与 **`hisat-3n`** 进行源码级与系统级加速，同时在调度层（workflow-level）挖掘可并行性，最终将端到端运行时间显著缩短。

---

## 1. 背景：为什么需要加速 m5C 检测工作流

RNA 会经历多种化学修饰，其中 **5-甲基胞嘧啶（m5C）** 在基因表达调控、RNA 稳定性、剪接与转运等方面具有重要作用，也与疾病（尤其癌症）机制相关。现有基于高通量测序（HTS）的检测方案（如 RNA-Bis-seq、UBS-seq）能做到单碱基分辨率，但由于依赖 C→T 转换信号，容易出现假阳性，同时其计算流程**步骤多、I/O 重、可扩展性差**，在真实数据上耗时往往成为瓶颈。

我们的目标是：**在不改变核心统计与筛选逻辑的前提下，最大化吞吐与多核可扩展性**，让该工作流更适用于批量数据处理与大规模研究。

---

## 2. 工作流概览：从 FASTQ 到最终 m5C 位点

工作流以测序得到的 FASTQ reads 为输入，整体包含以下阶段：

1. Reference index building（参考索引构建）
2. Data cleaning（数据清洗/截断过滤）
3. rRNA/tRNA filtering and genome alignment（先过滤非编码 RNA，再进行基因组比对）
4. Sorting and deduplication（排序与去重）
5. Site calling and filtering（位点识别与过滤）
6. Post-processing and Result Integration（跨重复整合与统计检验）

为了清晰描述与计时，我们将其划分为三段：

* **Part I：reference index building**（不计入计时区；可预先构建或下载索引）
* **Part II：数据处理与分析**（全流程主要耗时所在）
* **Part III：后处理与结果整合**（占比小，但仍可并行）

### 数据集

我们使用同一实验条件的三个 replicate：

* GSM7051146（replicate 1）
* GSM7051147（replicate 2）
* GSM7051148（replicate 3）

使用多个 replicate 可提升最终位点结果的可靠性，并允许更稳健的统计筛选（例如跨三个重复均满足阈值再保留）。

---

## 3. 性能剖析：瓶颈在哪里

我们先建立 baseline（不含 Part I）：

* **单节点串行 baseline：651.6 min**
* **replicate 并行 baseline：247.6 min**
* 三个 replicate 的 Part II 用时约：181 / 242 / 222 min；Part III 约 5.6 min

剖析结论很明确：

* **Part II 占据绝大多数运行时间**；
* Part II 中用时最多的两段是：

  * **rRNA/tRNA filtering and genome alignment**（核心工具：`hisat-3n` + `samtools`）
  * **site calling and filtering**（核心工具：`hisat-3n-table`，以及 `samtools view/bgzip` 等）

因此优化优先级是：
**先把工作流能并行的部分并行化 → 再深挖 `hisat-3n-table` 与 `hisat-3n` 的可扩展性与 I/O 瓶颈 → 最后处理 samtools 等 I/O 密集环节。**

---

## 4. 工作流层优化：把“独立性”转化为“并行度”

我们识别并利用了三类天然独立性：

### 4.1 replicate 级并行（跨节点/跨进程）

三个 replicate 在 Part II 完全独立，可分别在不同节点/进程上执行，最终汇总到 Part III，理论上可获得接近 **3×** 的加速（数据规模接近时更接近理想值）。

### 4.2 Part II 内部：Site calling & filtering 子任务并行

该阶段存在多条互不依赖的操作链（特别围绕 `hisat-3n-table` 的多个条件/输出），可由 Snakemake 自动并行或 bash 手工并行，实现该段约 **4×** 的局部加速。

### 4.3 Part III 并行

Part III 中 `join_pileup`、`stat_sample_background` 等可按样本并行，整体约 **1.6×** 加速。

---

## 5. 工具级优化 I：hisat-3n-table（~17×）

`hisat-3n-table` 是工作流中最耗时的组件之一，用于从比对结果生成统计表格/位点相关中间结果。其结构可概括为：
主线程加载 reference positions + 读入比对行 → 工作线程解析与更新 position → 主线程批量输出到输出池 → 输出线程写盘。

### 5.1 主要瓶颈

* 多个全局共享缓冲区（如 linepool/outputpool）带来大量 **lock/unlock** 与线程休眠（sleep），线程数增加反而可能更慢。
* 主线程在加载染色体与输出阶段承担过重工作（`loadNewChromosome` / `moveAllToOutput` 等热点）。
* 高频 malloc/free 造成内存管理开销放大（多线程下尤甚）。

### 5.2 核心优化策略

1. **调大缓冲区与 loadingBlockSize**：减少频繁同步与等待，降低主线程在 `appendingFinished / moveBlockToOutput / loadMore` 等路径上的触发频率。
2. **无锁队列替换安全队列**：对 `linepool/outputpool` 等共享队列使用 `moodycamel::ConcurrentQueue`，显著降低锁竞争；同时用原子计数弥补 `size_approx` 不精确的问题。对无并发需求的池（如 freepositionpool）改为普通队列。
3. **工作线程批量块处理**：一次取多行作为 block 处理，降低“多个线程同时打到相邻 position”引发的冲突概率。
4. **主线程热点加速**：

   * 延迟字符串赋值（chromosome 信息仅在最终输出需要时才赋值），减少主线程无效 string copy。
   * `moveAllToOutput` 中可并行的合法性检查与资源释放用 OpenMP 并行，缩短主线程停顿。
   * 将递归二分查找改为迭代版本，减少函数调用与栈开销（工作线程热点）。
5. **内存管理优化**：引入 `jemalloc` 改善多线程分配/释放性能；移除程序退出阶段冗余 delete（交给 OS 回收）；并在 RAM disk 上运行减少 I/O 影响。
6. **按染色体集合做进程级并行**：利用“染色体是最大独立单位”的特性切分输入，由多个进程各自处理不同染色体集合，最后合并输出，从而绕开主线程瓶颈并进一步提升可扩展性。

综合上述策略，我们在 `hisat-3n-table` 上获得约 **17×** 的加速收益。

---

## 6. 工具级优化 II：hisat-3n（~13.4×）

`hisat-3n` 用于过滤非编码 RNA 并进行比对，是 Part II 的另一大耗时项。其瓶颈集中在：多工作线程共享输入读取与输出写回接口，导致 I/O 相关锁竞争显著。

### 6.1 核心优化策略

1. **I/O 解耦：专用读线程 + 专用写线程**
   将读入与输出从工作线程中剥离，工作线程专注比对计算；读/写线程通过无锁队列与工作线程通信，显著减少锁竞争与同步开销。
2. **无锁输入/输出队列**
   读线程持续生产数据、写线程消费输出，工作线程只进行无锁出入队；并尽量在队列中传指针而非大对象，降低搬运成本；配合批量处理进一步提升吞吐。
3. **进程级并行：split fastq + 并行比对 + 合并输出**
   用 `seqkit split2` 切分输入 fastq，多进程并行跑 `hisat-3n`，最终用 `samtools cat` 合并结果。该策略在高核数机器上更容易接近线性扩展。

最终，`hisat-3n` 达到约 **13.4×** 的加速效果。

---

## 7. 其他组件：samtools 的 I/O 加速（~5×）

在本工作流中，`samtools view/sort/fastq` 等属于典型 **I/O 密集型** 操作。我们发现将中间大文件的读写放到 **RAM disk** 上执行，可显著减少磁盘访问延迟，从而在该类步骤上取得约 **5×** 的收益（示例以数据集 90 的 sort/fastq/view 计）。

---

## 8. 端到端结果：37.6 min

在整合工作流并行 + 关键工具加速后，端到端用时达到：**37.6 min**。

对比不同 baseline：

* 相对 replicate 并行 baseline（247.6 min）：约 **6.58×**
* 相对单节点串行 baseline（651.6 min）：约 **19×**

---

## 9. 复现与使用方式（README 友好版）

我们将完整工作流打包成可直接运行的结构（支持 bash / Snakemake 两种入口，Snakemake 版本可逐步完善）。典型目录包含数据、依赖工具、官方脚本以及一键运行脚本：

* `data_root/`：各样本 FASTQ 与中间文件
* `tool_root/`：hisat-3n、samtools、cutseq、umicollapse 等依赖
* `tempbin/`：整合处理脚本（含官方脚本）
* `bash/`：分步脚本与一键脚本
* `snakeversion/`：Snakemake 版本（配置与 Snakefile）

### Bash 运行（推荐快速上手）

* 若目录结构不变：直接 `bash ./test.sh` 一键跑完。
* 若路径不同：先 export 环境变量（如 `id_num / data_root / tool_root / tempbin / result_root`），再运行分步脚本；每个样本跑一次 `step1-14.sh`，最后 `step_final.sh` 合并结果。

---

## 10. 经验总结与可迁移点

本仓库的优化经验对“多线程生信工具/流水线”普遍适用：

* **先 profiling 再动手**：用 VTune/插桩将“感觉慢”转化为“热点函数与竞争点”。
* **并行优先级：workflow-level > tool-level**：能用任务并行吃掉的时间先吃掉。
* **锁竞争是可扩展性的第一杀手**：无锁队列 + 批量处理往往比盲目加线程更有效。
* **主线程热点决定上限**：当 worker 很快时，主线程加载/输出会成为硬瓶颈，需要结构性重构。
* **I/O 密集环节要从存储层下手**：RAM disk、管道、减少落盘次数往往立竿见影。
* **进程并行是多核机器的“第二条腿”**：当线程扩展受限时，按染色体/按 reads 切分进程并行非常有效。

---
