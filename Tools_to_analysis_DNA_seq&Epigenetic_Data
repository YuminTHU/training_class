### 2) DNA-seq  
#### (2.1) Mapping and QC  
- Remove adaptor: [cutadapt](https://cutadapt.readthedocs.io/en/stable/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- Mapping: [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [STAR](https://github.com/alexdobin/STAR)  
- QC: [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#### (2.2) Mutation
- Mutation discovery: [GATK](https://gatk.broadinstitute.org/hc/en-us), [Varscan](http://dkoboldt.github.io/varscan/)
- Mutation annotation: [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)
#### (2.3) Assembly
- denovo assembly software: [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
#### (2.4) CNV
- Whole Genome Seq: [Control-FREEC](http://boevalab.inf.ethz.ch/FREEC/)  
- Whole exome Seq: [CONTRA](http://contra-cnv.sourceforge.net/), [ExomeCNV](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3179661/)
#### (2.5) SV (structural variation)
- [Breakdancer](http://breakdancer.sourceforge.net/)

### 5) Epigenetic Data
#### (5.1)DNA Methylation 
 - Bi-sulfate data:
   - Review: [Katarzyna Wreczycka, et al. Strategies for analyzing bisulfite sequencing data. Journal of Biotechnology. 2017.](https://www.sciencedirect.com/science/article/pii/S0168165617315936)
   - Mapping: [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/), [BSMAP](https://github.com/zyndagj/BSMAPz)
   - Differential Methylation Regions (DMRs) detection: [methylkit](https://bioconductor.org/packages/release/bioc/html/methylKit.html), [ComMet](https://github.com/yutaka-saito/ComMet)
   - Segmentation of the methylome, Classification of Fully Methylated Regions (FMRs), Unmethylated Regions (UMRs) and Low-Methylated Regions (LMRs): [MethylSeekR](http://www.bioconductor.org/packages/release/bioc/html/MethylSeekR.html)
   - Annotation of DMRs: [genomation](https://bioconductor.org/packages/release/bioc/html/genomation.html), [ChIPpeakAnno](https://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html)  
   - Web-based service: [WBSA](http://wbsa.big.ac.cn/)
 - IP data:
   - Overview to CHIP-Seq: https://github.com/crazyhottommy/ChIP-seq-analysis
   - peak calling: [MACS2](https://github.com/taoliu/MACS/wiki/Advanced:-Call-peaks-using-MACS2-subcommands)
   - Peak annotation: [HOMER annotatePeak](http://homer.ucsd.edu/homer/ngs/annotation.html), [ChIPseeker](http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
   - Gene set enrichment analysis for ChIP-seq peaks: [GREAT](http://bejerano.stanford.edu/great/public/html/)
   
#### (5.2)DNAase-seq: 
- review : [Yongjing Liu, et al. Brief in Bioinformatics, 2019.](https://academic.oup.com/bib/article-abstract/20/5/1865/5053117?redirectedFrom=fulltext)
- peek calling:  [F-Seq](http://fureylab.web.unc.edu/software/fseq/)
- peek annotation: [ChIPpeakAnno](https://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html)
- Motif analysis: [MEME-ChIP](http://meme-suite.org/doc/meme-chip.html?man_type=web)
#### (5.3)ATAC-seq:
 - pipeline recommended by [Harward informatics](https://github.com/harvardinformatics/ATAC-seq)
 - peek calling: [MACS2](https://github.com/taoliu/MACS/wiki/Advanced:-Call-peaks-using-MACS2-subcommands)  
 - peak annotation: [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
 - Motif discovery: [HOMER](http://homer.ucsd.edu/homer/introduction/basics.html)
 
