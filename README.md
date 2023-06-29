# ConPosER - CONserved POSition identifiER  <img align="right" src="https://github.com/kevinsblake/ConPosER/blob/main/photos/hex/hex.png" width=300>
Identify conserved positions in amino acid multiple sequence alignments in R.

The sequence-structure-function paradigm follows that a protein’s primary sequence dictates its three-dimensional structure, which in turn influences function. Sequence determinants are generally restricted to a limited number of residues, evidenced by structurally similar proteins having low sequence identity and the ability to experimentally mutate many positions along the protein sequence with no impact on activity. Those positions that do affect activity when mutated tend to be located at the protein's core or at functional sites. It follows that residues which are conserved across a protein family are likely determinants of function, as any mutations at these sites abolish function which would exclude these mutants from selection and sequencing.

`ConPosER` generates a multiple sequence alignment (MSA) of input amino acid (AA) sequences, then identifies positions 100% conserved across all sequences. It can then plot the 2D location of these positions on a gene map. MSAs are generated using the [`msa`](https://bioconductor.org/packages/release/bioc/html/msa.html) package. 

`ConPosER` includes two functions:
1. `conposer_id()`: Takes input AA sequences, and outputs list of 100% conserved positions and the residue observed there.
2. `conposer_plot()`: Plots the 100% conserved positions on 2D gene map, plus a barplot of the number of AAs observed at each position. Calls a sub-function `conposer_geneplot()` which includes the formatting for the geneplot.

# Content

## Install

Install directly from GitHub:
```r
source("https://raw.github.com/kevinsblake/ConPosER/main/ConPosE.R")
```

Alternatively, can download and then install using the filepath:
```r
source("dir1/dir2/pathotype.R")
```

## Functions

### `conposer_id()`

#### Description

Function to list each conserved position and the amino acid at that positions.

#### Usage

```r
library(Biostrings)
library(msa)
library(dplyr)

conposer_id(seqs.file="sequences.fasta", msa=c("ClustalW", "ClustalOmega", "Muscle"), gap.lim=0.30)
```

#### Arguments

`seqs.file`	Filepath to fasta file containing input sequences.

`msa`		Specifies the multiple sequence alignment algorithm used to align the input sequences.

`gap.lim`	Specifies the fraction of sequences that can have a gap at a given position. Columns with gaps greater than this fraction are ignored. (Default = 0.30) 	

#### Examples

```
> df <- conposer_id("sequences.fasta", msa="ClustalOmega")
using Gonnet
> head(df)
  pos AA
1  39  G
2  44  G
3  51  L
4  62  E
5  70  R
6  73  G
```

### `conposer_plot()`

#### Description

Function for generating a gene plot of the MSA indicating the locations of the conserved positions, plus a barplot of the number of AAs observed at each position.

#### Usage

```r
library(Biostrings)
library(msa)
library(dplyr)

conposer_plot(seqs.file="sequences.fasta", filename="geneplot.pdf", linecol="black", gap.lim=0.30)
```

#### Arguments

`seqs.file`	Path and filename of fasta file containing input sequences.

`filename`	Path and filename of output geneplot. Default is "geneplot.pdf"

`linecol`	Color of lines indicating conserved positions. Default is "black."

`gap.lim`	Specifies the fraction of sequences that can have a gap at a given position. Columns with gaps greater than this fraction are ignored. (Default = 0.30) 

#### Examples

```r
> conposer_plot(seqs.file="sequences.fasta", filename="geneplot.pdf", linecol="red")
using Gonnet
pdf 
  2
```

![Ex1](https://github.com/kevinsblake/ConPosER/blob/main/photos/examples/example_01.png)

## References

- Bodenhofer _et al._ msa: an R package for multiple sequence alignment. _Bioinformatics._ 2015 Dec 15; 31(24):3997-9. doi: 10.1093/bioinformatics/btv494.
