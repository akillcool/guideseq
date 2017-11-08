# guideseq: GUIDE-Seq分析工具包
[![travis badge](https://travis-ci.org/aryeelab/guideseq.svg?branch=master)](https://travis-ci.org/aryeelab/guideseq)

guideseq包实现了GUIDE-Seq数据的数据预处理和分析流程。 它将原始测序读取（FASTQ）和参数清单文件（.yaml）作为输入，并生成一个带注释的无偏离位点表作为输出。

## 目录
- [功能](#features)
- [依赖](#dependencies)
- [配置](#setup)
	- [安装依赖](#install_dependencies)
	- [下载参考基因组](#reference_genome)
	- [下载并设置guideseq](#guideseq_setup)
	- [配置一个MiSeq以输出索引读数](#miseq)
- [运行完整的分析流程](#full_pipeline)
	- [快速开始](#quickstart)
	- [写一个Manifest文件](#write_manifest)
	- [一个完整的Manifest文件示例](manifest_example)
	- [流程输出](#pipeline_output)
- [运行独立的分析步骤](#)
	- [解复用](#demultiplex)
	- [UMItag读数](#umitag)
	- [合并](#consolidate)
	- [对齐](#align)
	- [识别](#identify)
	- [筛除](#filter)
	- [可视化](#visualize)
- [测试guideseq软件包](#testing)
	- [单步回归测试](#regression_tests)
	- [全面测试](#full_test)
	- [手动测试](#manual_test)
- [常见问题](#FAQ)
	- [如何使用解复用数据来运行流程？](#demultiplexed_run)
	- [能否不通过UMIs来分析数据？](#no_umis)


## 功能<a name="features"></a>


该软件包实现了一个由读取预处理模块和一个无偏离识别模块组成的流程。 预处理模块使用混合多样本测序运行中的原始读取（FASTQ）作为输入。 读数（FASTQ）被分离成样本特异性FASTQ，并使用独特的分子指数（UMI）条形码信息去除PCR重复。

![guideseq_flowchart](guideseq_flowchart.png)

具体步骤如下:

1. **样品分离**: 根据样品特定的双重索引条码序列将集中的多样品测序运行分离成样品特定的读取文件
2. **PCR重复合并**:共享相同的UMI和相同的基因组序列的前六个碱基的读数被认为是来源于相同的pre-PCR分子，因此被整合为一个一致的读数，以改善GUIDE-Seq读取计数的定量解释。
3. **读数对齐**：使用具有默认参数的BWA-MEM算法（Li.H，2009）将解复用的，合并的配对末端读数与参照基因组对齐。
4. **候选位点识别**: 把用标签特异性引物扩增读数的起始映射位置在全基因组范围内列表。 起始映射位置使用10-bp滑动窗口进行合并。 具有映射到+和 - 链或同一条链但用正向和反向标签特异性引物扩增的读数窗口被标记为潜在DSB的位点。 在每个标记的窗口中最频繁出现的起始映射位置的任一侧检索25bp的参考序列。 使用Smith-Waterman局部对齐算法将检索的序列与预期的目标序列对齐。
5. **错误阳性过滤**：过滤掉与目标靶序列超过6个错配或存在于背景对照中的脱靶切割位点。
6. **报告**：通过GUIDE-Seq读取计数排序的已识别脱靶目标在最终输出表中注释。 预计GUIDE-Seq阅读计数将与解理率大致成线性关系 (Tsai et al., *Nat Biotechnol.* 2015).
7. **可视化**：检测到的脱靶位点的对齐通过颜色编码的序列栅格可视化，如下所示：

![guideseq_flowchart](EMX1_visualization.png)

## 依赖<a name="dependencies"></a>
* Python (2.7)
* 需要的基因fasta文件 ([示例](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta))
* [`bwa`](<http://bio-bwa.sourceforge.net/>) 对齐工具
* [`bedtools`](<http://bedtools.readthedocs.org/en/latest/>) 基因算法工具


## 配置<a name="setup"></a>

### 安装依赖<a name="install_dependencies"></a>

要运行guideseq, 你必须安装所有需要的依赖包:

- **Python 2.7**: 如果某个版本与您的操作系统不兼容，我们建议您使用 [Anaconda](https://www.continuum.io/downloads) 科学用Python包。
- **Burrows-Wheeler Aligner (bwa)**: 您可以使用软件包管理器来安装bwa（例如，OSX上的`brew`或Ubuntu / Debian上的`apt-get`）, 或者从[bwa项目主页](http://bio-bwa.sourceforge.net/) 下载并自行编译源码。
- **Bedtools**: 您可以使用软件包管理器来安装Bedtools（例如`brew`或`apt-get`）, 或者从[Bedtools项目主页](http://bedtools.readthedocs.org/en/latest/content/installation.html) 下载并自行编译源码。

对于bwa和bedtools，请确保知道各个可执行文件的路径，因为它们需要在管道清单文件中指定。


### 下载参考基因组<a name="reference_genome"></a>

guideseq软件包需要一个参考基因组来进行读数对照。 您可以使用您选择的任何基因组，但是对于我们所有的测试和原始GUIDE-seq分析（Tsai et al。*Nature Biotechnol* 2015），我们使用hg19基因组 ([下载](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta))。 如果您下载到的是一个压缩包，请确认先解压FASTA文件。

### 下载并设置guideseq<a name="guideseq_setup"></a>

所有依赖安装完成后，仅需几个简单的步骤就能够下载并设置guideseq软件包：

1. 获取一份guideseq软件包的源代码。 你可以从github的 [release page](https://github.com/aryeelab/guideseq/releases)下载后解压, 或者通过运行 `git clone --recursive https://github.com/aryeelab/guideseq.git`来获取。
2. 进入guideseq目录，运行 `pip install -r requirements.txt`来安装guideseq依赖。

安装完所有的guideseq依赖后，你就可以准备开始使用guideseq了。

### 配置一个MiSeq以输出索引读数<a name="miseq"></a>

guideseq软件包需要读数索引以从MiSeq测序运行中进行读取合并。 默认的MiSeq Reporter设置不会生成索引读数（I1，I2）。 该功能可以通过添加一行代码

```xml
<add key="CreateFastqForIndexReads" value="1"> 
```

到位于Miseq Reporter安装目录下的``Miseq Reporter.exe.config`` 文件来启用。 默认安装路径为 ``C:\Illumina\MiSeqReporter``。 编辑后的config文件应该是下面这样：


```xml
<appSettings>
    ... [LEAVE EXISTING LINES UNCHANGED] ...
    <add key="CreateFastqForIndexReads" value="1"> 
</appSettings>
```

MiSeq Reporter服务需要重启使变更生效。 后面的GenerateFASTQ运行流程 (也许包括其他的工作流程)除了``R1``和``R2``读数之外，还将会生成``I1``和``I2``读数。这四个读数文件都需要进行guideseq分析。

更多说明参见Miseq Reporter User Guide的第29页。


## 运行完整的分析流程<a name="full_pipeline"></a>

### 快速开始<a name="quickstart"></a>

要运行完整的guideseq分析流程，你必须先创建一个包含声明所有流程输入的YAML清单文件。YAML清单文件创建完成之后，你可以通过

```
python /path/to/guideseq.py all -m /path/to/manifest.yaml
```
来运行整个流程。以下是如何编写清单文件的具体说明。

If you wish to run an example on our abridged test data, you can simply run
如果你想通过我们的简略测试数据来运行一个例子，你可以通过在guideseq根目录中运行

```
python guideseq/guideseq.py all -m test/test_manifest.yaml
```
来完成。 `test_manifest`假设`bwa`和`bedtools`都从系统PATH路径下运行。你可以在`test/output`文件夹下看到流程输出。

### 写一个Manifest文件<a name="write_manifest"></a>
运行guideseq软件包的端到端分析功能时，需要多个输入参数。 为了简化这些输入参数的格式并鼓励再现性，这些参数通过格式化为YAML清单文件被输入到流程中。YAML文件可以轻松完成键值对的读取。这使我们可以轻松地指定我们的参数。 文件中需要以下字段：

- `reference_genome`: 参考基因组FASTA文件的绝对路径。
- `output_folder`: 所有流程输出文件的保存文件夹，绝对路径。
- `bwa`: `bwa` 可执行文件的绝对路径。
- `bedtools`: `bedtools` 可执行文件的绝对路径。
- `undemultiplexed`: 未解复用的配对末端测序文件的绝对路径。所需的参数是：
	- `forward`: 包含正向读取的FASTQ文件的绝对路径。
	- `reverse`: 包含反向读取的FASTQ文件的绝对路径。
	- `index1`: 包含正向索引读取的FASTQ文件的绝对路径。
	- `index2`: 包含反向索引读取的FASTQ文件的绝对路径。

`undemultiplexed` 字段的示例：

```
undemultiplexed:
    forward: ../test/data/undemux.r1.fastq.gz
    reverse: ../test/data/undemux.r2.fastq.gz
    index1: ../test/data/undemux.i1.fastq.gz
    index2: ../test/data/undemux.i2.fastq.gz
```

- `samples`: 包含每个样本细节的嵌套字段。 必须至少指定两个样品：“对照（control）”样品（用于过滤背景脱靶位点）和至少一个处理样品。 所需的参数是：
	- `target`: 样品靶位点
	- `barcode1`: 正向条码
	- `barcode2`: 反向条码
	- `description`: 样品描述

`samples` 字段的示例：

```
samples:
    control:
        target:
        barcode1: CTCTCTAC
        barcode2: CTCTCTAT
        description: Control

    [SAMPLENAME]:
        target: GAGTCCGAGCAGAAGAAGAANGG
        barcode1: TAGGCATG
        barcode2: TAGATCGC
        description: EMX1
```

### 一个完整的Manifest文件示例<a name="manifest_example"></a>

下面是一个完整的YAML清单文件的例子。 随意复制它，并用您自己的实验数据替换参数。 请记住，您可以输入不止一个处理样本（例如下面的“EMX1”数据）。

```
reference_genome: test/test_genome.fa
output_folder: test/output

bwa: bwa
bedtools: bedtools

demultiplex_min_reads: 1000

undemultiplexed:
    forward: test/data/undemultiplexed/undemux.r1.fastq.gz
    reverse: test/data/undemultiplexed/undemux.r2.fastq.gz
    index1: test/data/undemultiplexed/undemux.i1.fastq.gz
    index2: test/data/undemultiplexed/undemux.i2.fastq.gz

samples:
    control:
        target:  
        barcode1: CTCTCTAC
        barcode2: CTCTCTAT
        description: Control

    EMX1:
        target: GAGTCCGAGCAGAAGAAGAANGG
        barcode1: TAGGCATG
        barcode2: TAGATCGC
        description: EMX_site1

```

### 流程输出<a name="pipeline_output"></a>

运行整个流程时，每一步的结果会在输出目录（即`output_folder`字段中设定的文件夹）中生成一个独立的文件夹。输出的文件夹及其各自的内容如下：


#### 输出文件夹
- `output_folder/demultiplexed`: 包含每个样本的四个未解复用的读取文件（正向，反向，index1，index2）。
- `output_folder/umitagged`: 包含每个样本的两个umitagged读取文件（正向，反向）。
- `output_folder/consolidated`: 包含每个样本的两个合并读取文件（正向，反向）。
- `output_folder/aligned`: 包含每个样本的一个对齐`.sam`文件。
- `output_folder/identified`: 包含每个样本的制表符分隔的`.txt`文件，每行有一个off-target标识。
- `output_folder/filtered`: 包含每个样品的一个制表符分隔的`.txt`文件，其中有已标识的DSB的每个样本，这些DSB是背景位点（而非脱靶位点）
- `output_folder/visualization`: 包含一个`.svg`矢量图像，表示每个样本的所有检测到的脱靶目标与目标点的对齐。


最终检测到的脱靶位点被放置在`output_folder/identified`文件夹中，每个YAML清单中指定的每个样本都有一个`.txt`文件。 这些脱靶样本文件的每一行中填充的字段如下所示：

#### 输出脱靶 `.txt` 字段：

- `BED Chromosome`: 窗口染色体
- `BED Min.Position`: 基于窗口0的开始位点
- `BED Max.Position`: 基于窗口0的结束位点
- `BED Name`: 窗口的名称
- `Filename`: 分析中使用的当前`.SAM`文件的名称。
- `WindowIndex`: 窗口的索引号
- `Chromosome`: 与在窗口中具有最大读数的位置相对应的染色体(与 `BED Chromosome`相匹配)
- `Position`: 窗口中读取次数最多的位点
- `Sequence`: 窗口序列，位于`Chromosome:Position`上下游25bp之间
- `+.mi`: 具有不同分子指数的正向读数
- `-.mi`: 具有不同分子指数的反向读数
- `bi.sum.mi`: `+.mi`和`-.mi`字段的总和 (GUIDE-seq读数)
- `bi.geometric_mean.mi`: `+.mi`和`-.mi`字段的几何平均值
- `+.total`: 正向映射读取的总数
- `-.total`: 反向映射读取的总数
- `total.sum`: `+.total`和`-.total`字段的总和。
- `total.geometric_mean`: `+.total`和`-.total`字段的几何平均值
- `primer1.mi`: 具有不同分子指数的正向引物扩增的读数数目
- `primer2.mi`: 具有不同分子指数的反向引物扩增的读数数目
- `primer.geometric_mean`: `primer1.mi`和`primer2.mi`字段的几何平均值
- `position.stdev`: 基因组窗口内位置的标准偏差
- `Off-Target Sequence`: 从参考基因组衍生的脱靶序列
- `Mismatches`: 预期目标序列与脱靶序列之间的错配数目
- `Length`: 目标序列的长度
- `BED off-target Chromosome`: 脱靶染色体
- `BED off-target start`: 脱靶0基准起始位置
- `BED off-target end`: 脱靶0基准结束位置
- `BED off-target name`: 脱靶名称
- `BED Score`: 符合标准BED格式的字段
- `Strand`: 表示检测到的脱离目标站点的链。 `+` 表示正向链， `-` 表示反向链
- `Cells`: 细胞类型
- `Target site`: 目标位点名称
- `Target Sequence`: 预期的目标位点序列（包括PAM）

解释这个输出和识别脱靶点的关键字段是： `BED off-target Chromosome`, `BED off-target start`, `BED off-target end`, `BED off-target name`, `BED off-target strand`, `Off-Target Sequence`, `bi.sum.mi`

#### 输出可视化

输出的可视化格式为`.svg`矢量格式，它是一个开放的图像标准，可以在任何现代的网页浏览器（例如Google Chrome，Apple Safari，Mozilla Firefox）中查看，并且可以在任何矢量图像编辑应用程序中查看和编辑（例如Adobe Illustrator）。 因为输出的可视化是矢量图像，所以它们可以无限放大或缩小，而不会损失质量，还可以轻松地编辑为形状。 这使得guideseq包所制作的图像非常适合海报，演示文稿和论文。

## 运行独立的分析步骤<a name="individual_steps"></a>

除了端到端的流程分析功能之外，guideseq软件包还允许分析的每一步都能单独运行。 在这里我们详细介绍了每一步所需的输入和预期的输出。 对于每一个步骤，我们都包含一个可以从guideseq根目录执行的命令实例，以对包含的样本数据运行该步骤。 这些命令的输出在`test/output`文件夹中。

### `解复用` 多样本合并测序(需要YAML清单文件)<a name="demultiplex"></a>

- **功能**: 对给定的在YAML清单中未复用的序列文件和条码样本，以FASTQ格式输出解复用的样本特定读数。 YAML清单中每个样本的正向，反向和两个索引文件被输出到`output_folder/consolidated`文件夹。
- **所需参数**:
	- `-m 或 --manifest`: 指定YAML清单文件的路径
- **运行示例**:
	- `python guideseq/guideseq.py demultiplex -m test/test_manifest.yaml`

### `umitag` 读数<a name="umitag"></a>

- **功能**: 对给定文件夹`output_folder/undemultiplexed`中的解复用文件（其输出路径（`output_folder`）在清单中指定），通过将UMI条形码序列添加到FASTQ读数名称头中来“标签化”读数，为随后的PCR重复读取整合步骤进行准备。 清单中每个样本的正向和反向文件被输出到 `output_folder/umitagged`文件夹。
- **所需参数**:
	- `--read1`: 正向解复用读数文件（FASTQ）的路径
	- `--read2`: 反向解复用读数文件（FASTQ）的路径
	- `--index1`: index1解复用读数文件（FASTQ）的路径
	- `--index2`: index2解复用读数文件（FASTQ）的路径
	- `--outfolder`: 输出文件文件夹的路径
- **运行示例**:

	```
	python guideseq/guideseq.py umitag --read1 test/data/demultiplexed/EMX1.r1.fastq \
	--read2 test/data/demultiplexed/EMX1.r2.fastq \
	--index1 test/data/demultiplexed/EMX1.i1.fastq \
	--index2 test/data/demultiplexed/EMX1.i2.fastq \
	--outfolder test/output/
	```

### `合并` PCR复制<a name="consolidate"></a>

- **功能**: 给出在清单中指定的未复用的序列文件和示例条形码，将合并的正向和反向读取输出到 `outfolder`。
- **所需参数**:
	- `--read1`: 正向umitagged读数文件（FASTQ）的路径
	- `--read2`: 反向umitagged读数文件（FASTQ）的路径
	- `--outfolder`: 输出文件文件夹的路径
- **可选参数**:
	- `--min_quality`: 合并中要考虑的最小读数质量
	- `--min_frequency`: 合并位点的最小读数频率
- **运行示例**:

	```
	python guideseq/guideseq.py consolidate --read1 test/data/umitagged/EMX1.r1.umitagged.fastq \
	 --read2 test/data/umitagged/EMX1.r2.umitagged.fastq \
	 --outfolder test/output/
	```

### `对齐` 位点到基因<a name="align"></a>

- **功能**: 给出合并后的正向和反向读数，使用`bwa`包执行序列与参考基因组的双末端对齐。将 `.sam` 对齐文件输出到 `outfolder`.
- **所需参数**:
	- `--bwa`: `bwa`的执行路径
	- `--genome`: 参考基因组FASTA文件的路径
	- `--read1`: 合并后的正向读数文件（FASTQ）的路径
	- `--read2`: 合并后的反向读数文件（FASTQ）的路径
	- `--outfolder`: 输出文件文件夹的路径
- **运行示例**:

	```
	python guideseq/guideseq.py align --bwa bwa --genome test/test_genome.fa\
	 --read1 test/data/consolidated/EMX1.r1.consolidated.fastq\
	 --read2 test/data/consolidated/EMX1.r2.consolidated.fastq\
	 --outfolder test/output/
	```

### `识别` 脱靶候选位点<a name="identify"></a>

- **功能**: 对指定的位点、参考基因组和目标序列给出其对齐sam文件，输出一个制表符分隔的、包含已识别的脱靶位点的`.txt`文件
- **所需参数**:
	- `--aligned`: 特定位点的对齐`.sam`文件的路径。
	- `--genome`: 参考基因组FASTA文件的路径
	- `--outfolder`: 输出文件文件夹的路径
	- `--target_sequence`: 样本中的指定序列 (对照样本（control）为空白)
- **可选参数**:
	- `--description`: 指定关于样本的附加信息。
- **运行示例**:

	```
	python guideseq/guideseq.py identify --aligned test/data/aligned/EMX1.sam\
	 --genome test/test_genome.fa --outfolder test/output/\
	 --target_sequence GAGTCCGAGCAGAAGAAGAANGG --description EMX1
	```

### `筛除` 背景DSB位点<a name="filter"></a>

- **功能**: 给定处理和控制样本的识别位点`.txt`文件，输出一个`.txt`文件，其格式与`识别`步骤输出的格式相同，其中包含筛选出的误报位点。
- **所需参数**:
	- `--bedtools`:  `bedtools`可执行文件的路径
	- `--identified`: 处理样本在`识别`步骤输出的`.txt`文件的路径。
	- `--background`: 对照样本（control）在`识别`步骤输出的`.txt`文件的路径。
	- `--outfolder`: 输出文件文件夹的路径
- **运行示例**:

	```
	python guideseq/guideseq.py filter --bedtools bedtools\
	 --identified test/data/identified/EMX1_identifiedOfftargets.txt\
	 --background test/data/identified/control_identifiedOfftargets.txt\
	 --outfolder test/output/
	```

### `可视化` 检测到的脱靶位点<a name="visualize"></a>

- **功能**: 给出一个已识别的脱靶位点`.txt`文件，输出对齐的脱靶位点可视化图。
- **所需参数**:
	- `--infile`:  输入的脱靶位点 `.txt` 文件路径
	- `--outfolder`: 包含输出的`.svg`图形的输出文件夹的路径
- **可选参数**:
	- `--title`: 指定可视化图形的标题，标题显示在图形顶部。用于海报和演示文稿。
- **运行示例**:

	```
	python guideseq/guideseq.py visualize --infile test/data/identified/EMX1_identifiedOfftargets.txt\
	 --outfolder test/output/ --title EMX1
	```

## 测试guideseq软件包<a name="testing"></a>

本着测试驱动开发的精神，我们为流程的每一步都编写了端到端的测试。 这些测试可以用来确保软件的功能按照期望运行。

注意：由于不同版本的`bwa`软件包之间存在差异，您必须使用`bwa v0.7.9a`来进行这些测试。 我们还建议您在运行这些测试时使用`bedtools v2.25.0`来保证一致性。

### 单步回归测试<a name="regression_tests"></a>

对于正在进行的测试和开发，我们为流程的每个步骤创建了一组简略的输入数据和预期的输出数据。 这样就可以快速测试流程的更改以进行特征回归。

要运行这些测试，你必须首先安装`nose`测试Python包。

```
pip install nose
```

然后，只需在guideseq根目录运行
```
nosetests
```
每个流程步骤的回归测试将会开始运行。

### 全面测试<a name="full_test"></a>

如果您有更多的时间，我们准备了一个bash脚本，用于从源代码下载和编译所有依赖项，下载新的参考基因组以及完整的GUIDE-seq测序运行数据，并对整个流程进行全面测试。 这个测试需要很长时间，但是我们需要在我们添加新版本之前运行它。

要运行全面测试，进入`guideseq/test`目录并运行

```
./large_test.sh
```

然后，坐看全面测试自动在终端中运行。

### 手动测试<a name="manual_test"></a>

如果你想通过分析流程运行一个完整的GUIDE-Seq数据集，你可以在下面找到，附带一个测试YAML清单（取决于你的依赖位置）：

```
http://aryee.mgh.harvard.edu/guideseq/data/guideseq_test_fastq.zip
```

上述文件应该与以下参考基因组一起使用：

```
http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta
```

## 常见问题<a name="FAQ"></a>

### 如何使用解复用数据来运行流程？<a name="demultiplexed_run"></a>

如果您已经有解复用数据，则可以在单独对数据解复用后再运行流程的每个步骤，如上面`运行独立的分析步骤`一节中所述。 请确认按以下顺序运行各个步骤：

- `umitag`
- `consolidate`
- `align`
- `identify`
- `filter`
- `visualize`

### 能否不通过UMIs来分析数据？<a name="no_umis"></a>

是。 如果您的读数不含UMI，那么可以按照上面的“单独运行分析步骤”一节中所述，从`对齐`步骤开始，对先前解复用的数据运行流程。  **请注意，我们还没有分析过这样的数据！** 我们怀疑PCR复制偏倚可能会影响GUIDE-Seq读数的定量解释，但尚未探索这一点。
