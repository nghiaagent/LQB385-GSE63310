# Import packages
library("edgeR")
library("limma")
library("tidyverse")

# Load data
dge <- readRDS(here::here("output/dge.rds"))
design <- readRDS(here::here("output/design.rds"))

# Process data prior to analysis
dge <- dge %>%
  # Remove lowly expressed genes
  .[filterByExpr(.), , keep.lib.sizes = FALSE] %>%
  # Calculate TMM normalization factors
  calcNormFactors(method = "TMM")
