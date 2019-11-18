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

> [STAR-Fusion的GitHub主页](https://github.com/STAR-Fusion/STAR-Fusion/wiki) 有详细的软件使用方法说明
>
> **Brian J. Haas**, et al. [STAR-Fusion: Fast and Accurate Fusion Transcript Detection from RNA-Seq.](https://www.biorxiv.org/content/10.1101/120295v1) **bioRxiv**, 2017.


其他可以用于分析融合基因的软件有：[Prada](http://bioinformatics.mdanderson.org/main/PRADA:Overview), [FusionCatcher](http://biorxiv.org/content/early/2014/11/19/011650), [SoapFuse](http://soap.genomics.org.cn/soapfuse.html), [TophatFusion](http://ccb.jhu.edu/software/tophat/fusion_index.html), [DISCASM/GMAP-Fusion](https://github.com/DISCASM/DISCASM/wiki)。


### 2)running steps

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
/BioII/lulab_b/chenyinghui/software/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR \
 --runThreadN 2 \
 --genomeDir /BioII/lulab_b/chenyinghui/database/Homo_sapiens/GRCh38_STAR_index \
 --readFilesIn /Share2/home/lulab/zhuyumin/share/zhuyumin/test/docker/StarFusionOut/cutadapt/SRR5712523_1.fastq.gz  /Share2/home/lulab/zhuyumin/share/zhuyumin/test/docker/StarFusionOut/cutadapt/SRR5712523_2.fastq.gz \
 --outFileNamePrefix /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/SRR5712523/SRR5712523. \
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

docker run -v /BioII:/BioII -v /Share2:/Share2 --rm trinityctat/starfusion /usr/local/src/STAR-Fusion/STAR-Fusion \
--CPU 2 \
--genome_lib_dir /Share2/home/lulab/zhuyumin/share/zhuyumin/test/docker/StarFusionOut/ctat_genome_lib_build_dir \
-J /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/SRR5712523/SRR5712523.Chimeric.out.junction \
--output_dir /BioII/lulab_b/chenyinghui/project/Docker/STAR-Fusion/SRR5712523_fusion

echo starfusion end `date`
```
