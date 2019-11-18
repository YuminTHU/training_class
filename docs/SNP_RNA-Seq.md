### SNP analysis in RNA-Seq
### 1) Introduction

本示例使用STAR将RNA测序数据比对到参考基因组后，使用GATK(4.0以上版本)进行SNP/INDEL的检测，最后使用ANNOVAR对SNP/INDEL进行注释。

#### (1) STAR

目前，可以用于将RNA测序数据（reads）比对到参考基因组的软件有：Bowtie、TopHat、HISAT、STAR等。

STAR在运行时候占用机器的内存较大，一般可达到20~30G，因此需要根据可用内存控制同时运行的STAR任务数量。


- [STAR的Github主页](https://github.com/alexdobin/STAR)
- 参考文献： **Alexander Dobin**, et al. [STAR: ultrafast universal RNA-seq aligner](https://academic.oup.com/bioinformatics/article/29/1/15/272537) _Bioinformatics_. 2012. 29(1): 15-21.

#### (2) GATK

GATK是Broad Institute开发的一款用于检测变异（SNP/INDEL）的软件，拥有较高的引用率（已有上万次引用）。

- [GATK的主页](https://software.broadinstitute.org/gatk/)
- [GATK Forum](https://gatkforums.broadinstitute.org/gatk/) （GATK开发人员与在该论坛回答用户疑问）
- 参考文献：**Aaron McKenna**, et al. [The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data.](https://genome.cshlp.org/content/20/9/1297.long) _Genome Research_. 2010. 20: 1297-1303.

#### (3) ANNOVAR

本示例中使用ANNOVAR进行变异位点信息注释。ANNOVAR是一款优秀的变异注释软件，注释速度快，且可以免费使用。用户可以选择下载公共数据库进行注释，也可以用自己制作的数据库文件（ANNOVAR接受BED/VCF格式）进行注释。

可以使用ANNOVAR提供的Perl脚本下载数据库，如下：

```bash
perl /BioII/lulab_b/chenyinghui/software/annovar/annovar/annotate_variation.pl \
-buildver hg38 \
-downdb \
-webfrom annovar \
refGene \ #下载的数据库名称
/BioII/lulab_b/chenyinghui/database/Homo_sapiens/annovar  #下载数据库存放路径
```

- 参考文献： **Wang K**, et al. [ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data](http://nar.oxfordjournals.org/content/38/16/e164) _Nucleic Acids Research_. 2010. 38:e164.

### 2) Running steps

#### (1) Alignment

GATK要求输入的SAM/BAM文件中有Read Group信息，因此我们需要在`--outSAMattrRGline`中填写相应的Read Group信息，其中`ID`为必填项。Read Group格式详见[SAM格式](http://samtools.github.io/hts-specs/SAMv1.pdf)

```bash
echo alignment start `date`
source /BioII/lulab_b/containers/singularity/wrappers/bashrc
STAR \
--genomeDir /BioII/lulab_b/chenyinghui/project/Docker/SNP/reference/Homo_sapiens_GRCh38_ch1_STAR_Index \
--runThreadN 1 \
--readFilesIn /BioII/lulab_b/chenyinghui/project/Docker/SNP/chr1.fq \ 
--readFilesCommand "gunzip -c" \
--outSAMtype BAM SortedByCoordinate \
--outSAMattrRGline ID:1 LB:exoRBase PL:ILLUMINA PU:unit1 SM:SRR5714908 \
--outFileNamePrefix /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/1.STAR_alignment/SRR5714908.
echo alignment end `date`
```

#### (2) MarkDuplicates

在建库的PCR过程中会形成一些重复的DNA/cDNA片段序列，这些重复序列被称为PCR duplicates。另外，在测序仪进行光学测量时候，也会形成一些光学重复，optical duplicates。如果变异位点位于这些重复的序列中，可能导致变异频率偏高，因此需要对重复序列进行标记，使得后续变异检测软件可以识别这些重复序列。

我们使用GATK包含的MarkDuplicates功能模块进行重复序列标记。（该功能原属于Picard软件，后被GATK收录）

```bash
echo 2.MarkDuplicates start `date`
/BioII/lulab_b/chenyinghui/software/GATK/gatk-4.1.3.0/gatk  MarkDuplicates \
--java-options '-Xmx2G'
--INPUT  /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/1.STAR_alignment/SRR5714908.Aligned.sortedByCoord.out.bam \
--OUTPUT /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/2.MarkDuplicates/SRR5714908.sorted.MarkDup.bam \
--CREATE_INDEX true \
--METRICS_FILE /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/2.MarkDuplicates/SRR5714908.sorted.MarkDup.bam.metrix \
--VALIDATION_STRINGENCY SILENT
echo 2.MarkDuplicates end `date`
```

#### (3) SplitNCigarReads

SAM/BAM文件的第6列为CIGAR表达式，用来表示该序列各个位置的碱基的比对情况。

由于RNA测序得到的reads中的一部分可能会跨越不同的exon，对于这些reads的CIGAR表达式中会出现`N`值。使用GATK中的SplitNCigarReads功能将这种reads切分为k+1个reads（K为CIGAR中的`N`的数量）

```bash
echo 3.SplitNCigarReads start `date`

/BioII/lulab_b/chenyinghui/software/GATK/gatk-4.1.3.0/gatk  SplitNCigarReads \
--java-options '-Xmx2G'
-R /BioII/lulab_b/chenyinghui/project/Docker/SNP/reference/Homo_sapiens.GRCh38.ch1.faa \
-I /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/2.MarkDuplicates/SRR5714908.sorted.MarkDup.bam \
-O /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/3.SplitNCigarReads/SRR5714908.sorted.MarkDup.SplitNCigar.bam

echo 3.SplitNCigarReads end `date`

```

#### (4) HaplotypeCaller

该步骤是正式使用GATK进行变异检测的步骤。

```bash
echo 4.HaplotypeCaller start `date`

/BioII/lulab_b/chenyinghui/software/GATK/gatk-4.1.3.0/gatk HaplotypeCaller \
--java-options '-Xmx2G' \
--input  /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/3.SplitNCigarReads/SRR5714908.sorted.MarkDup.SplitNCigar.bam \
--output /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/4.HaplotypeCaller/SRR5714908.raw.vcf.gz \
--reference /BioII/lulab_b/chenyinghui/project/Docker/SNP/reference/Homo_sapiens.GRCh38.ch1.fa \
-dont-use-soft-clipped-bases \
--standard-min-confidence-threshold-for-calling 20

echo 4.HaplotypeCaller end `date`

```



#### (5) VariantFiltration

我们可以根据变异的聚集程度、变异的链偏好性、变异的平均质量水平、位点测序深度等指标进行过滤。

```bash
echo 5.VariantFiltration start `date`

/BioII/lulab_b/chenyinghui/software/GATK/gatk-4.1.3.0/gatk VariantFiltration \
--java-options '-Xmx1G' \
--variant /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/4.HaplotypeCaller/SRR5714908.raw.vcf.gz \
--output /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/5.VariantFiltration/SRR5714908.filtered.vcf.gz \
--reference /BioII/lulab_b/chenyinghui/project/Docker/SNP/reference/Homo_sapiens.GRCh38.ch1.fa \
--window 35 \
--cluster 3 \
--filter-name 'FS' \
--filter 'FS > 30.0' \
--filter-name 'QD' \
--filter 'QD < 2.0' \
--filter-name 'DP' \
--filter 'DP < 5'
zcat /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/5.VariantFiltration/SRR5714908.filtered.vcf.gz  | awk -F '\t' '{if ($0 ~ "#" || $7 == "PASS") print $0 }' - >/BioII/lulab_b/chenyinghui/project/Docker/SNP/output/5.VariantFiltration/SRR5714908.filtered.clean.vcf

echo 5.VariantFiltration end `date`

```

`--window 35 \ --cluster3` :  （变异的聚集程度）if there are more than 3 variants cluster in 35 bp window, these variants will be filtered.

`--filter-name 'FS' \ --filter 'FS > 30.0' `:  （[变异的链偏好性（strand bias）](https://gatkforums.broadinstitute.org/gatk/discussion/8056/fisher-s-exact-test)）Phred-scaled p-value （-log10(p-value) ） using Fisher's exact test to detect strand bias. Phred-score closer to 0 means there is a lower chance of there being bias. Higher FS values therefore indicate more bias.

`--filter-name 'QD'`: （变异的平均质量水平）Variant Confidence/Quality by Depth.   

`--filter-name 'DP'`：（位点的测序深度）read depth that cover this site.



`--filter-name`中的过滤参数的具体含义可以参考GATK生成的VCF文件中开头的注释部分。用户可以根据VCF文件中出现的其他指标对变异进行更多样的过滤筛选。

值得指出的是，满足用户所设置的过滤表达式（如平均质量QD低于2: `--filter 'QD < 2.0'`）的变异才是我们需要过滤的变异。这些需要被过滤的“不合格”变异仍然会被保留在VCF文件中，但是在VCF第6列 `QUAL`中会被标注过滤的原因（平均质量QD太低，则标记为`QD`），通过筛选的、合格的变异位点会被标记`PASS`。
我们可以用`awk`等命令去除VCF中不合格变异，保留合格变异。
 

#### (6) Annotation

我们可以使用软件以及数据库对得到的变异进行注释，可以获得注释信息包括：变异的位置（位于哪个基因？ 位于exon/intron/UTR ?）、在人群（例如东亚人群）中的频率，临床意义（Pathogenic/Benign）等等。这些注释信息可以帮助研究人员对变异的重要性作出判断。

ANNOVAR将注释分为gene-based annotation、filter-based annotation、Region-based annotation等，在下方的`--operation`参数中分别对应`g`与`f`。

可以用于ANNOVAR注释的公共数据库以及对应的注释类型可以参考[ANNOVAR的官网](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)中的解释。

```bash
echo 6.Annotation start `date`
perl /BioII/lulab_b/chenyinghui/software/annovar/annovar/table_annovar.pl \
/BioII/lulab_b/chenyinghui/project/Docker/SNP/output/5.VariantFiltration/SRR5714908.filtered.clean.vcf  \ #需要注释的VCF文件
/BioII/lulab_b/chenyinghui/database/Homo_sapiens/annovar/ \ #注释使用的数据库所处的文件夹路径
--buildver hg38 \ #比对所用的人类参考基因组的版本，决定了变异坐标是基于何种版本参考基因组
--outfile /BioII/lulab_b/chenyinghui/project/Docker/SNP/output/6.Annotation/SRR5714908.annotated.variants \  #输出文件路径与文件名前缀
--remove \  #自动删除中间文件
--nastring . \  #注释空缺值用.代替
--vcfinput \  #输入文件为VCF格式
--thread 2 \
--protocol ensGene,cosmic70,gnomad211_genome,exac03,clinvar_20190305,avsnp150 \ #使用的数据库
--operation g,f,f,f,f,f \  #数据库对应的注释类型

echo 6.Annotation end `date`
```
