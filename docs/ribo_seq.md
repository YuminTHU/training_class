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

#### 2.Generating P-site track
