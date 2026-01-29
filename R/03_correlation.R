here::i_am("R/03_correlation.R")

# Import packages
library("limma")
library("magrittr")
library("tidyverse")
library("WGCNA")

# Load data
dge_voom <- readRDS(here::here("output/dge_voom.rds"))
dge <- readRDS(here::here("output/dge.rds"))

# Extract expression data
dge_correlations <- dge_voom$E %>%
  # Change row names to gene symbols
  set_rownames(dge_voom$genes$SYMBOL) %>%
  # Transpose so genes are on columns
  t() %>%
  cor()

sel <- c("Xkr4", "Sulf1", "Gpc1")

x <- dge_correlations[
  rownames(dge_correlations) %in% sel,
  colnames(dge_correlations) %in% sel
]

# Save data
saveRDS(dge_correlations, file = "output/dge_correlations.rds")
