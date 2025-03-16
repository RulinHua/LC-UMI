# LC-UMI:LOESS regression-corrected UMI count

## Overview
LcUMI is a method that adjusts UMI counts to correct cross-cell amplification variability, enabling absolute gene expression quantification at the single-cell level.

This is an R package for enhancing the precision and reliability of single-cell transcriptomics. 

## Installation
LcUMI is available in Bioconductor. In addition, one can install the development version from the Github repository:
``` r
## To install the package from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LcUMI")

## To install the development version from the Github repo:
devtools::install_github("RulinHua/LcUMI")
```
