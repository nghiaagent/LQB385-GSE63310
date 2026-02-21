# A webR and Jupyter-based transcriptome analysis learning module


## Introduction

This repository contains the source code to reproduce a webR +
Jupyter-based bioinformatics learning module.

The module is a guided walkthrough for analysing the
[GSE63310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63310)
RNA-seq dataset, containing sorted mice mammary cells from across three
cell types in the mammary luminal lineage. The module guides students
through exploratory analysis, differential expression analysis, and gene
set analysis.

The module is distributed as two files:

- `analysed_GSE63310.Rdata`: Prepared dataset to minimise setup time.
- `LQB385_BW3.ipynb`: Jupyter notebook containing code for all analyses.

## Prerequisites

To reproduce the dataset, make sure the following are installed on your
system:

- R 4.5.0
- (recommended) Positron.

## Reproduciblity

To reproduce the dataset:

- Install all prerequisites software.
- Clone the repository; Positron is recommended for this purpose.
- Run the following in an R terminal to download all packages:

``` r
renv::restore()
```

- Run the following in an R terminal to perform all analyses:

``` r
source(here::here("R/99_build_dataset.R"))
```

## Student instructions

1.  Go to <https://jupyter.r-wasm.org>
2.  Upload the provided `LQB385_BW3.ipynb` and `analysed_GSE63310.Rdata`
    files.
3.  Open the `LQB385_BW3.ipynb` notebook and run all cells sequentially
    to reproduce the analyses.

## Limitations

- As the webR-enabled Jupyter interface has no internet access, images
  cannot be hyperlinked from the web. Instead, all images were
  base64-encoded and embedded directly into the notebook.
- This solution cannot render HTML tables from R code cells (unlike
  e.g. Quarto Live). Instead, they were rendered as text using `kable`.
