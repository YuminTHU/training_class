## sequence motif analysis
### workflow
![](../assets/structure_motif.pipeline.png)
### 0. get interested sequence and control sequence as sequence motif analysis
### 1. GraphProt
GraphProt supports input sequences in fasta format. The viewpoint mechanism sets viewpoints to all nucleotides in uppercase letters, nucleotides in **lowercase letters are only used for RNA structure predictions.**
#### 1.1 get lowercase letters
```
awk '{if ($1~/^>/) print;else print tolower($1)}' interested_promoter.fa \
  >interested_promoter_lowercase.fa
```
#### 1.2 get control lowercase letters
```
awk '{if ($1~/^>/) {print;} else {print tolower($1)}}' interested_promoter.control \
  >interested_promoter_lowercase.control
```
#### 1.3 run graphprot
##### 1.3.1 get optimized parameters
```
GraphProt.pl \
--action ls \
--fasta    interested_promoter_lowercase.fa \
--negfasta interested_promoter_lowercase.control \
--prefix   interested_promoter_optim
```
##### 1.3.2 train a model
```
GraphProt.pl \
--action train \
--fasta    interested_promoter_lowercase.fa \
--negfasta interested_promoter_lowercase.control \
--prefix   interested_promoter
--params   interested_promoter_optim.params
```
##### 1.3.3 get structure motif
```
GraphProt.pl \
--action motif \
--model  interested_promoter.model \
--params interested_promoter_optim.params \
--fasta  interested_promoter_lowercase.fa \
--prefix interested_promoter
```

### 2. BEAM
#### 2.1 transfer to RNA sequence
```
awk '{if ($1~/^>/) print;else {gsub(/T/,"U",$1);print;}}' interested_promoter.fa >interested_promoter_rna.fa
awk '{if ($1~/^>/) print;else {gsub(/T/,"U",$1);print;}}' interested_promoter.control >interested_promoter_rna.control
```
#### 2.2 prepare fastB format
##### 2.1 prepare dot-bracket format file by RNAfold
Compute the best (MFE) structure for this sequence (primary sequence with dot-bracket)
```
RNAfold  < interested_promoter_rna.fa >dot.fa
```
screenshot of output:dot.fa
```
```
##### 2.2 get file with BEAR notation ---> fastB (fastBEAR)
```
awk '/^>/ {print; getline; print; getline; print $1}' dot.fa >dot_to_encode.fa
java -jar BearEncoder.new.jar dot_to_encode.fa  BEAMready.fa
```
screenshot of output: BEAMready.fa


#### 2.3 get structure motif
```
java -jar BEAM_release1.6.1.jar -f BEAMready.fa -g bg.fa -w 10 -W 40 -M 3 
```
screenshot of output

#### 2.4 visualize motif
#### 2.4.1 install weblogo
```
pip install weblogo
```
##### 2.4.2 plot with suggested command
we should not only specify -o option with format suffix, we should also specify -F with format
```
# -o:Output file (default: stdout)
# -F:Format of output: eps (default), png, png_print, pdf, jpeg, svg, logodata
weblogo -a 'ZAQXSWCDEVFRBGTNHY' -f BEAMready_m1_run1_wl.fa -D fasta \
-o out.jpeg -F jpeg --composition="none" \
-C red ZAQ 'Stem' -C blue XSW 'Loop' -C forestgreen CDE 'InternalLoop' \
-C orange VFR 'StemBranch' -C DarkOrange B 'Bulge' \
-C lime G 'BulgeBranch' -C purple T 'Branching' \
-C limegreen NHY 'InternalLoopBranch'
```
##### 2.4.3 example output

###other recommended tools
RNApromo
RNAcontext
