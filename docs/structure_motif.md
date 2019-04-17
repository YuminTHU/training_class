## sequence motif analysis
### workflow
![](../assets/structure_motif.pipeline.png)
### 0. get interested sequence and control sequence as sequence motif analysis
### 1. GraphProt:modelling binding preferences of RNA-binding proteins
#### 1.1 get optimized parameters

#### 1.2 train a model
#### 1.3 get structure motif with graphprot



### 2. RNApromo
#### 2.1 download
https://genie.weizmann.ac.il/pubs/rnamotifs08/64bit_exe_rnamotifs08_motif_finder.tar.gz
#### 2.2 predict structure motifs
```
rnamotifs08_motif_finder.pl -positive_seq input_pos_seq.fa -output_dir Output
```
##### input
Positive sequences - a fasta format file containing the sequences to predict motifs on.
```
>Pos_1
ATAAGAGACCACAAGCGACCCGCAGGGCCAGACGTTCTTCGCCGAGAGTCGTCGGGGTTTCCTGCTTCAACAGTGCTTGGACGGAACCCGGCGCTCGTTCCCCACCCCGGCCGGCCGCCCATAGCCAGCCCTCCGTCACCTCTTCACCGCACCCTCGGACTGCCCCAAGGCCCCCGCCGCCGCTCCA
```
#### 2.3 find known motifs
After learning a motif, you can search a database of sequences to find positions that match the motif you learned. To do that you need to first match a secondary structure to each of the input sequences in your database, either using existing structure prediction algorithms, or using some other information.
#### 2.3.1 Produce a likelihood score for each sequence in the database.
```
rnamotifs08_motif_match.pl database.tab -cm model.cm
```
##### input
The database is then specified in the following format: <id> <sequence> <structure>
database.tab
```
seq_1	AUAAGAGACCACAAGCGACCCGCAGGGCCAGACGUUCUUCGCCGAGAGUCGUCGGGGUUUCCUGCUUCAACAGUGCUUGGACGGAACCCGGCGCUCGUUCCCCACCCCGGCCGGCCGCCCAUAGCCAGCCCUCCGUCACCUCUUCACCGCACCCUCGGACUGCCCCAAGGCCCCCGCCGCCGCUCCA	.............((((..(.((.(((((.(((((((((....))))).))))(((((.((((((...(((((((((.((......)).)))))).))).........(((.(((........))).)))..................))).....)))..)))))..)))))..)).).))))...
```
##### output
```
seq_2:0    19.7698
seq_7:0    19.3706
seq_3:0    19.1064
seq_1:0    18.073
seq_5:0    16.5508
seq_9:0    14.5906
seq_4:0    10.3685
seq_10:0   9.15077
seq_6:0    6.81294
seq_8:0    0.233537
```
#### 2.3.2 Produce a likelihood score for the best motif position in each sequence in the database, and the position itself.
##### output
```
seq_2:0    19.7698    33      48      UUCAACAGUGUUUGGA        (((((......)))))        <<<<<,,,,,,>>>>>
seq_7:0    19.3706    104     119     GGGAGCAGUGUCUUCC        (((((......)))))        <<<<<,,,,,,>>>>>
seq_3:0    19.1064    16      31      GUCCUCAGUGCAGGGC        (((((......)))))        <<<<<,,,,,,>>>>>
seq_1:0    18.073     30      52      GACGUUCUUCGCCGAGAGUCGUC (((((((((....))))).)))) <<<<<<<<<,--->>>>>,>>>>
seq_5:0    16.5508    104     119     AGCUACAGUGUUAGCU        (((((......)))))        <<<<<,,,,,,>>>>>
seq_9:0    14.5906    32      47      GAGCCAGUGUGUUUCU        ((((......))))..        <<<-,,,,,,->>>,,
seq_4:0    10.3685    7       21      UUGUCAGUGCACAAA         ((((......)))).         <<<<,,,,,,>>>>,
seq_10:0   9.15077    133     152     CAACCUCCACCUUCUGGGUU    .(((((.........)))))    ,<----,,,,,,,,,---->
seq_6:0    6.81294    1       16      UAUGGAGAUUUCCAUA        (((((......)))))        <<<<<,,,,,,>>>>>
seq_8:0    0.233537   95      115     ACACCCCAGCCCUGCAGUGUA   ((((..((....))..)))).   <<<<,,--,,,,---->>>>,
```
