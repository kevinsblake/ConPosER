# ConPosER - CONserved POSition identifiER
Identify conserved positions in amino acid multiple sequence alignments in R.

The sequence-structure-function paradigm follows that a protein’s primary sequence dictates its three-dimensional structure, which in turn influences function. Sequence determinants are generally restricted to a limited number of residues, evidenced by structurally similar proteins having low sequence identity and the ability to experimentally mutate many positions along the protein sequence with no impact on activity. Those positions that do affect activity when mutated tend to be located at the protein's core or at functional sites. It follows that residues which are conserved across a protein family are likely determinants of function, as any mutations at these sites abolish function which would exclude these mutants from selection and sequencing.

`ConPosER` generates a multiple sequence alignment (MSA) of input amino acid (AA) sequences, then identifies positions 100% conserved across all sequences. It can then plot the 2D location of these positions on a gene map. MSAs are generated using the [`msa`](https://bioconductor.org/packages/release/bioc/html/msa.html) package. 

`ConPosER` includes two functions:
1. `conposer_id()`: Takes input AA sequences, and outputs list of 100% conserved positions and the residue observed there.
2. `conposer_plot()`: Plots the 100% conserved positions on 2D gene map. Calls a sub-function `conposer_geneplot()` which includes the formatting for the geneplots.

# Content


## Install

Install directly from GitHub:
```r
source("https://raw.github.com/kevinsblake/ConPosER/main/ConPosE.R")
```

## Functions

### `conposer_id()`

#### Description

Function to list each conserved position and the amino acid at that positions.

#### Usage

```r
library(Biostrings)
library(msa)

seqs.file <- readAAStringSet("data/sequences/sequences.fasta")

conposer_id(seqs.file, msa=c("ClustalW", "ClustalOmega", "Muscle"))
```

#### Arguments

`seqs.file`	Path and filename fasta file containing input sequences

`msa`	Specifies the multiple sequence alignment algorithm used to align the input sequences.

#### Examples

### `conposer_plot()`

#### Description

Function for generating a gene plot of the MSA, indicating the positions of the conserved positions in 2D space.

#### Usage

```r
library(Biostrings)
library(msa)

seqs.file <- readAAStringSet("data/sequences/sequences.fasta")

conposer_plot(seqs.file, filename="geneplot.pdf", linecol="black", isFilter=FALSE)

```

#### Arguments

`seqs.file`	Path and filename of fasta file containing input sequences
`filename`	Path and filename of output geneplot. Default is "geneplot.pdf"
`linecol`	Color of lines indicating conserved positions. Default is "black"
`isFilter`	Specifies if positions with gaps in >5% of the sequences in the alignment are masked. Masking these positions can avoid spuriously calling a conserved position that is really a gap in most sequences. Default = FALSE, which shows every position.

#### Examples

## References

- Bodenhofer _et al._ msa: an R package for multiple sequence alignment. _Bioinformatics._ 2015 Dec 15; 31(24):3997-9. doi: 10.1093/bioinformatics/btv494.
