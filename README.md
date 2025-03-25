# LC-UMI:LOESS regression-corrected UMI count

## Overview
LCUMI is a method that adjusts UMI counts to correct cross-cell amplification variability, enabling absolute gene expression quantification at the single-cell level.

This is an R package for enhancing the precision and reliability of single-cell transcriptomics. 

## Installation
LCUMI is available in Bioconductor. In addition, one can install the development version from the Github repository:
``` r

## To install the development version from the Github repo:
install.packages("devtools")
devtools::install_github("RulinHua/LC-UMI")
library("LCUMI")
```
