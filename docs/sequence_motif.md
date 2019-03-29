### sequence motif analysis
#### 1. get UTR or promoter sequence
##### 1.1 install R package GenomicFeatures and biozhuoer tools (cnode)
GenomicFeatures package used to extract needed sequence

biozhuoer tools used to concat sequences of the same UTR or promoter
```
source("http://www.bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")
library(GenomicFeatures)

if (!('devtools' %in% .packages(T))) install.packages('devtools');
devtools::install_github('dongzhuoer/biozhuoer');
```
##### 1.2 generate txdb object
There are many functions for us to get genme annotation file:
```
gtf_file="/BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/gencode.v27.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
```
##### 1.3 get 3'UTR & 5'UTR site range
```
utr5p = fiveUTRsByTranscript(txdb, use.names=T)
utr3p = threeUTRsByTranscript(txdb, use.names=T)

utr3p.df=as.data.frame(utr3p)
utr5p.df=as.data.frame(utr5p)

write.table(utr3p.df, "utr3p.info", row.names=FALSE, sep='\t',quote=FALSE )
write.table(utr5p.df, "utr5p.info", row.names=FALSE, sep='\t' ,quote=FALSE)
```
##### 1.4 get promoter site range
```
promoter=promoters(txdb)
promoter.df=as.data.frame(promoter)
write.table(promoter.df, "promoter.info", row.names=FALSE, sep='\t' ,quote=FALSE)
```

#### 2. intersect with interested genes
##### 2.1 interested 3'UTR
```
sort -t $'\t' -k 2 utr3p.info|join -o 1.3 2.1 1.2 1.9 1.4 1.5 1.6 1.7 1.8 1.10 -t $'\t' -1 2 -2 2 - \
  <(cut -f 1 SC2_SF2.ct.dn.1_0.01.protein_coding |sort |join -t $'\t' -1 1 -2 1 - <(sort -t $'\t' -k 1 \
  <(grep -o -P -e "gene_id.*; transcript_id.*?;" /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/gencode.v27.annotation.gtf |sort |uniq|sed -e 's/gene_id "//' -e 's/"; transcript_id "/\t/' -e 's/";//'))|sort -t $'\t' -k 2  ) |\
  sort -t $'\t' -k 1 >interested_three_prime_UTR.info
```
explanation for interested_three_prime_UTR.info:
```
column1:chr
column2:gene
column3:transprict
column4:exon_name
column5:start
column6:end
column7:width(seq length)
column8:strand
column9:exon_id
column10:exon_rank
```
##### 2.2 interested promoter
```
sort -t $'\t' -k 7 promoter.info|join -o 1.1 2.1 1.7 1.2  1.3 1.4 1.5 1.6 -t $'\t' -1 7 -2 2 - \
  <(cut -f 1 SC2_SF2.ct.dn.1_0.01.protein_coding |sort |join -t $'\t' -1 1 -2 1 - \
  <(sort -t $'\t' -k 1 <(grep -o -P -e "gene_id.*; transcript_id.*?;" /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/gencode.v27.annotation.gtf |sort \
  |uniq|sed -e 's/gene_id "//' -e 's/"; transcript_id "/\t/' -e 's/";//' ))|sort -t $'\t' -k 2  ) |\
  sort -t $'\t' -k 1 >interested_promoter.info
```

#### 3. convert to bed format
#### 4. get genome sequence
#### 5. generate random sequence as background sequence
#### 6. motif enrichment


