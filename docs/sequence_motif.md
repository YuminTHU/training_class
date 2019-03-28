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
#### 2. intersect with interested genes
#### 3. convert to bed format
#### 4. get genome sequence
#### 5. generate random sequence as background sequence
#### 6. motif enrichment


