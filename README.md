# guideseq: GUIDE-Seq分析工具包
[![travis badge](https://travis-ci.org/aryeelab/guideseq.svg?branch=master)](https://travis-ci.org/aryeelab/guideseq)

guideseq包实现了GUIDE-Seq数据的数据预处理和分析管道。 它将原始测序读取（FASTQ）和参数清单文件（.yaml）作为输入，并生成一个带注释的无偏离位点表作为输出。

## 目录
- [功能](#features)
- [依赖](#dependencies)
- [配置](#setup)
	- [安装依赖](#install_dependencies)
	- [下载参考基因组](#reference_genome)
	- [下载并设置guideseq](#guideseq_setup)
	- [配置一个MiSeq以输出索引读数](#miseq)
- [Running the Full Analysis Pipeline](#full_pipeline)
	- [Quickstart](#quickstart)
	- [Writing A Manifest File](#write_manifest)
	- [A Full Manifest File Example](manifest_example)
	- [Pipeline Outputs](#pipeline_output)
- [Running Analysis Steps Individually](#)
	- [Demultiplex](#demultiplex)
	- [UMItag](#umitag)
	- [Consolidate](#consolidate)
	- [Align](#align)
	- [Identify](#identify)
	- [Filter](#filter)
	- [Visualize](#visualize)
- [Testing the guideseq Package](#testing)
	- [Single-Step Regression Tests](#regression_tests)
	- [Full Large Test](#full_test)
	- [Manual Testing](#manual_test)
- [Frequently Asked Questions](#FAQ)
	- [How do I Run the Pipeline with Demultiplexed Data?](#demultiplexed_run)
	- [Can I analyze data without UMIs?](#no_umis)


## 功能<a name="features"></a>


该软件包实现了一个由读取预处理模块和一个无偏离识别模块组成的流水线。 预处理模块使用混合多样本测序运行中的原始读取（FASTQ）作为输入。 读数（FASTQ）被分离成样本特异性FASTQ，并使用独特的分子指数（UMI）条形码信息去除PCR重复。

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


## Running the Full Analysis Pipeline<a name="full_pipeline"></a>

### Quickstart<a name="quickstart"></a>

To run the full guideseq analysis pipeline, you must first create a manifest YAML file that describes all pipeline inputs. Once you have done so, you can simply run

```
python /path/to/guideseq.py all -m /path/to/manifest.yaml
```
to run the entire pipeline. Below are specific instructions detailing how to write the manifest file.

If you wish to run an example on our abridged test data, you can simply run

```
python guideseq/guideseq.py all -m test/test_manifest.yaml
```
from the guideseq root directory. The `test_manifest` assumes that both the `bwa` and `bedtools`executables are in your system PATH. You will see the pipeline results outputted to the `test/output` folder.

### Writing A Manifest File<a name="write_manifest"></a>
When running the end-to-end analysis functionality of the guideseq package, a number of inputs are required. To simplify the formatting of these inputs and to encourage reproducibility, these parameters are inputted into the pipeline via a manifest formatted as a YAML file. YAML files allow easy-to-read specification of key-value pairs. This allows us to easily specify our parameters. The following fields are required in the manifest:

- `reference_genome`: The absolute path to the reference genome FASTA file.
- `output_folder`: The absolute path to the folder in which all pipeline outputs will be saved.
- `bwa`: The absolute path to the `bwa` executable
- `bedtools`: The absolute path to the `bedtools` executable
- `undemultiplexed`: The absolute paths to the undemultiplexed paired end sequencing files. The required parameters are:
	- `forward`: The absolute path to the FASTQ file containing the forward reads.
	- `reverse`: The absolute path to the FASTQ file containing the reverse reads.
	- `index1`: The absolute path to the FASTQ file containing the forward index reads.
	- `index2`: The absolute path to the FASTQ file containing the reverse index reads.

An example `undemultiplexed` field:

```
undemultiplexed:
    forward: ../test/data/undemux.r1.fastq.gz
    reverse: ../test/data/undemux.r2.fastq.gz
    index1: ../test/data/undemux.i1.fastq.gz
    index2: ../test/data/undemux.i2.fastq.gz
```

- `samples`: A nested field containing the details of each sample. At least two samples must be specified: a "control" sample (to be used to filter out background off-target sites) and at least one treatment sample. The required parameters are:
	- `target`: The sample targetsites
	- `barcode1`: The forward barcode
	- `barcode2`: The reverse barcode
	- `description`: A description of the sample

An example `samples` field:

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

### A Full Manifest File Example<a name="manifest_example"></a>

Below is an example of a full manifest file. Feel free to copy it and replace the parameters with your own experiment data. Remember that you can input more than just one treatment sample (e.g. the "EMX1" data below).

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

### Pipeline Output<a name="pipeline_output"></a>

When running the full pipeline, the results of each step are outputted to the `output_folder` in a separate folder for each step. The output folders and their respective contents are as follows:


#### Output Folders
- `output_folder/demultiplexed`: Contains the four undemultiplexed reads files (forward, reverse, index1, index2) for each sample.
- `output_folder/umitagged`: Contains the two umitgged reads files (forward, reverse) for each sample.
- `output_folder/consolidated`: Contains the two consolidated reads files (forward, reverse) for each sample.
- `output_folder/aligned`: Contains an alignment `.sam` file for each sample.
- `output_folder/identified`: Contains a tab-delimited `.txt` file for each sample with an identified off-target in each row.
- `output_folder/filtered`: Contains a tab-delimited `.txt` file for each sample containing the identified DSBs that are background sites (not off-targets)
- `output_folder/visualization`: Contains a `.svg` vector image representing an alignment of all detected off-targets to the targetsite for each sample.


The final detected off-target sites are placed in the `output_folder/identified` folder, with one `.txt` file for each sample specified in the manifest. The fields that are populated in each row of these off-target files are specified below:

####Output Off-Targets `.txt` Fields:

- `BED Chromosome`: Window chromosome
- `BED Min.Position`: Window 0-based start position
- `BED Max.Position`: Window 0-based end position
- `BED Name`: Name of window 
- `Filename`: The name of the current `.SAM` file used in analysis.
- `WindowIndex`: Index number of window
- `Chromosome`: Chromosome corresponding to position with maximum reads in window (matches `BED Chromosome`)
- `Position`: Position with maximum number of reads in window
- `Sequence`: The window sequence, starting 25 bp upstream and ending 25 bp downstream of `Chromosome:Position`
- `+.mi`: Number of forward reads with distinct molecular indices
- `-.mi`: Number of reverse reads with distinct molecular indices
- `bi.sum.mi`: Sum of the `+.mi` and `-.mi` fields (GUIDE-seq Read Count)
- `bi.geometric_mean.mi`: Geometric mean of the `+.mi` and `-.mi` fields
- `+.total`: Total number of forward mapping reads 
- `-.total`: Total number of reverse mapping reads 
- `total.sum`: Sum of `+.total` and `-.total` fields
- `total.geometric_mean`: Geometric mean of the `+.total` and `-.total` fields
- `primer1.mi`: Number of reads amplified by forward primer with distinct molecular indices
- `primer2.mi`: Number of reads amplified by reverse primer with distinct molecular indices
- `primer.geometric_mean`: Geometric mean of the `primer1.mi` and `primer2.mi` fields
- `position.stdev`: Standard deviation of positions within genomic window
- `Off-Target Sequence`: Off-target sequence derived from genome reference
- `Mismatches`: Number of mismatches between the intended target sequence and the off-target sequence
- `Length`: Length of the target sequence
- `BED off-target Chromosome`: Off-target chromosome
- `BED off-target start`: Off-target 0-based start position
- `BED off-target end`: Off-target 0-based end position
- `BED off-target name`: Off-target name
- `BED Score`: Field to conform to standard BED format
- `Strand`: Indicates the strand of detected off-target site. `+` for forward strand and `-` for reverse strand
- `Cells`: Cell type
- `Target site`: Targetsite name
- `Target Sequence`: Intended target site sequence (including PAM)

The key fields for interpreting this output and identifying off-target sites are: `BED off-target Chromosome`, `BED off-target start`, `BED off-target end`, `BED off-target name`, `BED off-target strand`, `Off-Target Sequence`, `bi.sum.mi`

#### Output Visualizations

The outputted visualizations are in the `.svg` vector format, which is an open image standard that can be viewed in any modern web browser (e.g. Google Chrome, Apple Safari, Mozilla Firefox), and can be viewed and edited in any vector editing application (e.g. Adobe Illustrator). Because the output visualizations are vector images, they can be scaled up or down infinitely without a loss in quality, and can also be edited as shapes with ease. This makes the images produced by the guideseq package ideal for posters, presentations, and papers.

## Running Analysis Steps Individually<a name="individual_steps"></a>

In addition to end-to-end pipeline analysis functionality, the guideseq package also allows for every step fo the analysis to be run individually. Here we have detailed the required inputs and expected outputs of each step. For each step, we have included a "runnable example" command that can be executed from the guideseq root directory to run that step on the included sample data. These "runnable example" snippets put their output in the `test/output` folder.

### `demultiplex` Pooled Multi-Sample Sequencing (Manifest Required)<a name="demultiplex"></a>

- **Functionality**: Given undemultiplexed sequence files and sample barcodes specified in the manifest, output the demultiplexed sample-specific reads in FASTQ format. The forward, reverse, and two index files for each sample in the manifest	 are outputted to the `output_folder/consolidated` folder.
- **Required Parameters**:
	- `-m or --manifest`: Specify the path to the manifest YAML file
- **Runnable Example**:
	- `python guideseq/guideseq.py demultiplex -m test/test_manifest.yaml`

### `umitag` Reads<a name="umitag"></a>

- **Functionality**: Given the demultiplexed files in the folder `output_folder/undemultiplexed` (where `output_folder` is specified in the manifest), 'tag' the reads by adding the UMI barcode sequence to the FASTQ read name header in preparation for the subsequent PCR duplicate read consolidation step. The forward and reverse files for each sample in the manifest are outputted to the `output_folder/umitagged` folder.
- **Required Parameters**:
	- `--read1`: Path to the forward demultiplexed reads file (FASTQ)
	- `--read2`: Path to the reverse demultiplexed reads file (FASTQ)
	- `--index1`: Path to the index1 demultiplexed reads file (FASTQ)
	- `--index2`: Path to the index2 demultiplexed reads file (FASTQ)
	- `--outfolder`: Path to the folder in which the output files will be saved
- **Runnable Example**:

	```
	python guideseq/guideseq.py umitag --read1 test/data/demultiplexed/EMX1.r1.fastq \
	--read2 test/data/demultiplexed/EMX1.r2.fastq \
	--index1 test/data/demultiplexed/EMX1.i1.fastq \
	--index2 test/data/demultiplexed/EMX1.i2.fastq \
	--outfolder test/output/
	```

### `consolidate` PCR Duplicates<a name="consolidate"></a>

- **Functionality**: Given undemultiplexed sequence files and sample barcodes specified in the manifest, output the consolidated forward and reversed reads to the `outfolder`.
- **Required Parameters**:
	- `--read1`: Path to the forward umitagged reads file (FASTQ)
	- `--read2`: Path to the reverse umitagged reads file (FASTQ)
	- `--outfolder`: Path to the folder in which the output files will be saved
- **Optional Parameters**:
	- `--min_quality`: The minimum quality of a read for it to be considered in the consolidation
	- `--min_frequency`: The minimum frequency of a read for the position to be consolidated
- **Runnable Example**:

	```
	python guideseq/guideseq.py consolidate --read1 test/data/umitagged/EMX1.r1.umitagged.fastq \
	 --read2 test/data/umitagged/EMX1.r2.umitagged.fastq \
	 --outfolder test/output/
	```

### `align` Sites to Genome<a name="align"></a>

- **Functionality**: Given the consolidated forward and reverse reads, execute a paired-end mapping of the sequences to the reference genome using the `bwa` package. Outputs an alignment `.sam` file to the `outfolder`.
- **Required Parameters**:
	- `--bwa`: Path to the `bwa` executable
	- `--genome`: Path to the reference genome FASTA file
	- `--read1`: Path to the consolidated forward read FASTQ file
	- `--read2`: Path to the consolidated reverse read FASTQ file
	- `--outfolder`: Path to the folder in which the output files will be saved
- **Runnable Example**:

	```
	python guideseq/guideseq.py align --bwa bwa --genome test/test_genome.fa\
	 --read1 test/data/consolidated/EMX1.r1.consolidated.fastq\
	 --read2 test/data/consolidated/EMX1.r2.consolidated.fastq\
	 --outfolder test/output/
	```

### `identify` Off-target Site Candidates<a name="identify"></a>

- **Functionality**: Given the alignment samfile for a given site, a reference genome, and a target sequence, output a tab-delimited `.txt` file containing the identified off-target sites.
- **Required Parameters**:
	- `--aligned`: Path to the site-specific alignment `.sam` file.
	- `--genome`: Path to the reference genome FASTA file.
	- `--outfolder`: Path to the folder in which the output files will be saved.
	- `--target_sequence`: The sequence targeted in the sample (blank for control sample)
- **Optional Parameters**:
	- `--description`: Specify additional information about the sample.
- **Runnable Example**:

	```
	python guideseq/guideseq.py identify --aligned test/data/aligned/EMX1.sam\
	 --genome test/test_genome.fa --outfolder test/output/\
	 --target_sequence GAGTCCGAGCAGAAGAAGAANGG --description EMX1
	```

### `filter` Background DSB Sites<a name="filter"></a>

- **Functionality**: Given the identified site `.txt` files for a treatment and control samples, output a `.txt` file in the same format as outputted by the `identify` step containing the sites filtered out as false-positives.
- **Required Parameters**:
	- `--bedtools`: Path to the `bedtools` executable
	- `--identified`: Path to the `.txt` file outputted by the `identify` step for a treatment sample.
	- `--background`: Path to the `.txt` file outputted by the `identify` step for a control sample.
	- `--outfolder`: Path to the folder in which the output files will be saved.
- **Runnable Example**:

	```
	python guideseq/guideseq.py filter --bedtools bedtools\
	 --identified test/data/identified/EMX1_identifiedOfftargets.txt\
	 --background test/data/identified/control_identifiedOfftargets.txt\
	 --outfolder test/output/
	```

### `visualize` Detected Off-Target Sites<a name="visualize"></a>

- **Functionality**: Given an identified off-target sites `.txt` file, output an alignment visualization of the off-target sites.
- **Required Parameters**:
	- `--infile`:  Path to the input `.txt.` off-targets file
	- `--outfolder`: Path to the outputted folder containing the outputted `.svg` graphic
- **Optional Parameters**:
	- `--title`: Specify the title of the visualization, to be printed at the top of the graphic. Useful for posters and presentations.
- **Runnable Example**:

	```
	python guideseq/guideseq.py visualize --infile test/data/identified/EMX1_identifiedOfftargets.txt\
	 --outfolder test/output/ --title EMX1
	```

## Testing the guideseq Package<a name="testing"></a>

In the spirit of Test-Driven Development, we have written end-to-end tests for each step of the pipeline. These can be used to ensure that the software is running with expected functionality.

NOTE: Due to differences in sorting between different versions of the `bwa` package, you must be using `bwa v0.7.9a` for these tests to work. We also recommend that you use `bedtools v2.25.0` when running these tests for consistency's sake.

### Single-Step Regression Tests<a name="regression_tests"></a>

For ongoing testing and development, we have created an abridged set of input data and expected output data for each step of the pipeline. This way, changes to the pipeline can be quickly tested for feature regression.

To run these tests, you must first install the `nose` testing Python package.

```
pip install nose
```

Then, from the guideseq root directory, simply run

```
nosetests
```

and the regression tests for each pipeline step will be run.

### Full Large Test<a name="full_test"></a>

If you have more time, we have prepared a bash script that downloads and compiles all dependencies from source, downloads a fresh reference genome and a full GUIDE-seq sequencing data run, and performs a full test of the entire pipeline. This test takes a long time, but we require that it be run before we tag a new release.

To run the full large test, enter the `guideseq/test` folder and run

```
./large_test.sh
```

Then, sit back and watch the full large testing process unfold automatically in the terminal.

### Manual Testing<a name="manual_test"></a>

If you wish to run a full GUIDE-Seq dataset through the analysis pipeline, you may find it and a test manifest (to be altered depending on your dependency locations) here:

```
http://aryee.mgh.harvard.edu/guideseq/data/guideseq_test_fastq.zip
```

which should be used with the following reference genome:

```
http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta
```

## Frequently Asked Questions<a name="FAQ"></a>

### How do I Run the Pipeline with Demultiplexed Data?<a name="demultiplexed_run"></a>

If you already have demultiplexed data, you can run the pipeline on the data by running each step after demultiplexing individually, as described in the "Running Analysis Steps Individually" section above. Be sure to run the individual steps in the following orders:

- `umitag`
- `consolidate`
- `align`
- `identify`
- `filter`
- `visualize`

### Can I analyze data without UMIs?<a name="no_umis"></a>

Yes. If your reads do not have UMIs, you can run the pipeline on previously demultiplexed data as described in the "Running Analysis Steps Individually" section above, starting with the `align` step. **Note that we have not analyzed such data ourselves!** We suspect that PCR duplication bias may affect the quantitative interpretion of GUIDE-Seq read counts, but have not explored this.
