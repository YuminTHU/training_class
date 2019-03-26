### pro-processing
#### 0. create annotation
```
scripts/create_annotation.sh -G annotation_fly/dmel-all-r6.18.gtf -f annotation_fly/dmel-all-chromosome-r6.18.fasta  -o annotation_fly  -s scripts;
```
##### input files
1. <annotation.gtf> : the annotation gtf should contain start_codon and stop_codon information,eg: dmel-all-r6.18.gtf
2. <genome.fasta> : genome fasta ,eg: dmel-all-chromosome-r6.18.fasta
3. <annotation_dir> : the directory for all the annotation output
4. <scripts_dir> : the directory of all the scripts in the package

##### output files
annotation directory, including :
1. start_codon.bed : the bed file annotating start codon
2. final.ORFs : all identified ORFs, eg: FBtr0300105_0_31_546 where FBtr0300105 refers to the transcript, 0 refers to the reading frame relative to the start of transcript, 31 refers to the start site, 546 refers to the stop codon.

#### 1. P-site determination
核糖体上具有一系列与蛋白质合成有关的结合位点与催化位点，分别为A位点(aminoacyl-site，A-site)，P位点(peptidyl-site，P-site)和E位点(exit-site，E- site)先后于tRNA发生结合。P位点是肽段翻译延长的主要场所，在该位点上tRNA将携带的氨基酸移交给旁边的肽段从而使得肽段序列发生延长。为了能够更加明显的观察到3-nt的周期性，在处理Ribo-seq数据时我们参考之前已发表的方法，对每一条Ribo-seq比对上的测序片段转换为其对应的P-site位点。

This step determines the P-site position for each Ribo-seq reads length by overlapping with the annotated start codons from previous step
```
scripts/P-site_determination.sh -i GSE52799/SRR1039770.sort.bam -S annotation_fly/start_codon.bed -o GSE52799 -n SRR1039770 -s scripts;
```
##### input files
1. <Ribo_bam> : secondary alignment removed to ensure one genomic position per aligned read and sorted
2. annotation :
  <start_codon.bed> : annotated start site start_codon.bed. It is generated in the create_annotation.sh step.
3. <out_dir> : the directory of the output result, eg: GSE52799
4. <study_name> : the name of all the output file, default: test. eg: SRR1039770
5. <scripts_dir>	: the directory of all the scripts in the package

##### output files
P-site directory, including :
1. name.psite1nt.txt : the Ribo-seq reads length and its corresponding P-sites position(= offset + 1). It may look this this :
```
30  13
```
2. name.psite.pdf : the PDF displaying the histogram of aggregated reads
我们收集了所有已被注释的起始密码子并将 这些起始密码子和Ribo-seq 比对上的序列进行重合，分别计算Ribo-seq序列的5’端偏离起始密码子的第一个碱基A的距离(offset)。根据 Ribo-seq测序片段长度的不同，我们进一步将Ribo-seq片段分成多个组分。在每个长度对应的组分里，作出Ribo-seq片段5’端偏离起始密码子A的距离(offset)的直方图。

贴个图

每一行代表不同长度的 Ribo-se0q 测序片段的直方图。例如在 29nt 长度的 Ribo-seq 片段中，我们可以明显的看到10950在0 距离为 13nt 的位点含有一个峰值(peak)。鉴于大部分核糖体会在翻译起始位点21000停滞较多的时间，因此对于 29nt 长的 Ribo-seq片段，其 P-site 位点的定义应该4代06500表直方图中绝大多数的核糖体，因此我们将 P-site 位点应该定义为峰值最高的第 13 个碱基(13nt)的位置。

#### 2.Generating P-site track
