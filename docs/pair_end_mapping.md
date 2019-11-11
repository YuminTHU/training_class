### 3a\) mapping

```bash
STAR --genomeDir ./STAR_TAIR10_index --readFilesIn Ath_1.fq Ath_2.fq --outFileNamePrefix ./output/ --outSAMtype BAM SortedByCoordinate

```

* `--genomeDir`  specifies path to the genome directory where genome indices where generated.
* `--readFilesIn`  name(s) (with path) of the files containing the sequences to be mapped (e.g.
RNA-seq FASTQ files). If using Illumina paired-end reads, the read1 and read2 files have to
be supplied. 
* `--outFileNamePrefix`  all output files are written in the current directory.
* `--outSAMtype BAM SortedByCoordinate` output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools
sort command.
