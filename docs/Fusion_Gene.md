### Fusion gene analyisis in RNA-Seq

### 1) Introduction

融合基因（Fusion gene）是指两个单独的基因全部或部分序列融合而形成的嵌合基因。 融合伴侣可能贡献5’UTR、编码区或3’多腺苷酸信号。形成融合基因的原因有：染色体间易位、缺失、倒位等。

最为经典的案例是**慢性粒细胞白血病**中**BCR-ABL1**的基因融合，治疗慢性粒细胞白血病的药物，**伊马替尼/格列卫**，其作用靶点就是该融合基因。

近期发现，**ZBTB16-ABL1**融合基因是导致形成急性T淋巴细胞白血病的关键基因。

> **Bing Chen**, et al. [Identification of fusion genes and characterization of transcriptome features in T-cell acute lymphoblastic leukemia](https://www.pnas.org/content/115/2/373). **_PNAS_**,  Jan 2018, 115 (2) 373-378.

STAR-Fusion是一款利用RNA-Seq数据检测人类融合基因的软件。

STAR-Fusion提供了Docker镜像，以方便用户使用。

```bash
docker pull trinityctat/starfusion
```
在本示例中，我们使用STAR-Fusion的Docker镜像进行分析。如果您不使用Docker镜像而是自行安装，请查看[STAR-Fusion的安装指南](https://github.com/STAR-Fusion/STAR-Fusion/wiki/installing-star-fusion)。

- [STAR-Fusion的GitHub主页](https://github.com/STAR-Fusion/STAR-Fusion/wiki) 有详细的软件使用方法说明

> **Brian J. Haas**, et al. [STAR-Fusion: Fast and Accurate Fusion Transcript Detection from RNA-Seq.](https://www.biorxiv.org/content/10.1101/120295v1) **bioRxiv**, 2017.


其他可以用于分析融合基因的软件有：[Prada](http://bioinformatics.mdanderson.org/main/PRADA:Overview), [FusionCatcher](http://biorxiv.org/content/early/2014/11/19/011650), [SoapFuse](http://soap.genomics.org.cn/soapfuse.html), [TophatFusion](http://ccb.jhu.edu/software/tophat/fusion_index.html), [DISCASM/GMAP-Fusion](https://github.com/DISCASM/DISCASM/wiki)。


### 2) Running steps

#### (1) 下载数据库

我们需要从Broad Institute数据库网站下载STAR-Fusion所需要的参考基因组与注释文件，选择“plug-n-play”压缩文件进行下载。下载地址如下：

https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/

下载后将其命名为CTAT_resource_lib.tar.gz ，解压。

#### (2) 运行STAR-Fusion

STAR-Fusion可以直接以Fastq为输入文件进行融合基因分析；也可以使用STAR的Chimeric.out.junction文件作为STAR-Fusion的输入文件。
下面分别介绍使用这2种不同输入文件进行分析的方法。

##### 2.1 输入文件为Fastq

由于STAR运行时会占用较大内存（RAM），约20～30G；如果STAR-Fusion加了`--FusionInspector validate `参数可能会使内存总占用达到～40G，因此当我们从fastq开始使用STAR-fusion分析时需要合理控制并行运行的STAR-Fusion任务数量。

```bash
#假设CTAT_resource_lib文件夹与reads_1.fq.gz、reads_2.fq.gz都在当前目录
docker run -v `pwd`:/data \ #将当前目录挂载为Docker的/data目录
--rm trinityctat/starfusion \ #当分析任务结束后，立即删除容器
/usr/local/src/STAR-Fusion/STAR-Fusion \
    --left_fq /data/reads_1.fq.gz \
    --right_fq /data/reads_2.fq.gz \
    --genome_lib_dir /data/ctat_genome_lib/ \
    -O /data/StarFusionOut \

```

##### 2.2 输入文件为Chimeric.out.junction

使用STAR将Fastq比对到参考基因组上，输出Chimeric.out.junction文件

```bash
echo STAR start `date`
source /BioII/lulab_b/containers/singularity/wrappers/bashrc
/BioII/lulab_b/chenyinghui/software/STAR/STAR-2.7.3a/bin/Linux_x86_64_static/STAR \
 --runThreadN 4 \
 --genomeDir /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/ctat_genome_lib_build_X_docker/ref_genome.fa.star.idx \
 --readFilesIn /Share2/home/lulab/zhuyumin/share/zhuyumin/test/docker/StarFusionOut/cutadapt/SRR5712523_1.fastq.gz  /Share2/home/lulab/zhuyumin/share/zhuyumin/test/docker/StarFusionOut/cutadapt/SRR5712523_2.fastq.gz \
 --outFileNamePrefix /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/SRR5712523_X/SRR5712523. \
 --outReadsUnmapped None \
 --readFilesCommand "gunzip -c" \
 --outSAMstrandField intronMotif \
 --outSAMunmapped Within \
 --chimSegmentMin 12 \
 --chimJunctionOverhangMin 12 \
 --chimOutJunctionFormat 1 \
 --alignSJDBoverhangMin 10 \
 --alignMatesGapMax 100000 \
 --alignIntronMax 100000 \
 --alignSJstitchMismatchNmax 5 -1 5 5 \
 --outSAMattrRGline ID:SRR5712523 \
 --chimMultimapScoreRange 3 \
 --chimScoreJunctionNonGTAG -4 \
 --chimMultimapNmax 20 \
 --chimNonchimScoreDropMin 10 \
 --peOverlapNbasesMin 12 \
 #--twopassMode Basic 在cnode上运行，该参数加上以后内存会溢出，建议注释掉
 --peOverlapMMp 0.1 

echo STAR end `date`

```

以Chimeric.out.junction为输入文件，用STAR-Fusion进行融合基因分析

```bash
echo starfusion start `date`
docker run -v /BioII:/BioII --rm trinityctat/starfusion /usr/local/src/STAR-Fusion/STAR-Fusion \
--CPU 2 \
--genome_lib_dir /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/ctat_genome_lib_build_X_docker \
-J /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/SRR5712523_X/SRR5712523.Chimeric.out.junction \
--output_dir /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/SRR5712523_fusion_X_docker


echo starfusion end `date`
```
如果运行STAR-Fusion时候有如下报错：

ERROR, ctat genome lib: ctat_genome_lib_build_X does not validate and may require indices to be rebuilt.

那么需要运行如下脚本，来重建index。因为STAR-Fusion的index在不同系统之间是不通用的，详见[该网址](https://github.com/STAR-Fusion/STAR-Fusion/wiki/rebuild_ctat_genome_lib_indices)

```bash
docker run -v /BioII:/BioII --rm trinityctat/starfusion /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/util/rebuild_indices.pl /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/ctat_genome_lib_build_X_docker
```
### 2) Results

```bash
#FusionName	JunctionReadCount	SpanningFragCount	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots
CLIC2--AC234781.1	3	0	ONLY_REF_SPLICE	CLIC2^ENSG00000155962.13	chrX:155334371:-	AC234781.1^ENSG00000224216.1	chrX:155335249:-	YES_LDAS	0.1327	GT	1.8062	AG	1.8323	["INTRACHROMOSOMAL[chrX:0.00Mb]","LOCAL_REARRANGEMENT:-:[201]"]
```

__JunctionReads__: 跨越融合位点(junction site)的reads数量   
   
__SpanningFragCount__: 对于某一成对的reads（mate pair），如果其中1条read完整地比对在融合位点上游的融合基因内，另一条read完整地比对在融合位点下游的融合基因内，则该成对的reads记为spanning fragments。SpanningFragCount表示spanning fragments数量。   
根据融合断点的位置和reads的长度，有时可能会出现SpanningFragCount 等于0，所有反映融合事件的reads都以JunctionReads形式出现的情况。   
   
__FFPM__: 支持Chimeric RNA的reads数目取决于融合转录的表达量以及测序的reads数量。随着测序总体数据量的升高，测序导致的重复的Chimeric RNA的reads数量也会升高。因此，我们需要对支持融合事件的reads数量进行标准化，即 FFPM (fusion fragments per million total reads)。STAR-Fusion中默认的FFPM阈值为0.1，意味着在总reads数为10 Million前提下，至少需要有1条read支持该Chimeric RNA。该默认参数可以有效过滤假阳性结果。如果需要更改FFPM过滤标准，可以在STAR-Fusion `--min_FFPM`参数中进行设置。
   
__LargeAnchorSupport__: 如果融合位点两侧25bp被junction reads覆盖，则记为"YES_LDAS"(LDAS = long double anchor support)。如果融合事件缺少LargeAnchorSupport以及spanning fragment支持，那么该Chimeric RNA结果很可能是假阳性。   
   
__SpliceType__: 断点是否出现在参考转录本结构注释文件提供的外显子连接点(exon junction site)。   
   
__LeftBreakEntropy, RightBreakEntropy__: 融合位点两侧15bp的外显子碱基的香农熵([Shannon entropy](http://bearcave.com/misl/misl_tech/wavelets/compression/shannon.html))。最大熵为2，表示复杂度最高。最低为0(表示15个碱基为单一核苷酸)。低熵位点通常应被视为不太可靠的融合位点。   
   
__annot__: 最后一列"annots"利用 FusionAnnotator(与STAR-Fusion绑定的程序)为Chimeric RNA提供了一个简化的注释。对于人类基因组，fusion注释信息基于[CTAT_HumanFusionLib](https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki)数据库，其中包括许多已知与癌症相关的融合基因的信息资源，以及一些在正常组织中也频繁出现的、与癌症无关的融合基因（这些注释信息有助于我们过滤掉假阳性的融合基因）。
