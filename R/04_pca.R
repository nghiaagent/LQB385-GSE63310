here::i_am("R/04_pca.R")

# Import packages
library(PCAtools)
library(limma)
library(tidyverse)

# Load data
dge_voom <- readRDS(here::here("output/dge_voom.rds"))
dge <- readRDS(here::here("output/dge.rds"))

# Perform PCA
plot_pca <- dge_voom$E %>%
  pca(
    metadata = dge$samples,
    removeVar = 0.9
  ) %>%
  biplot(
    colby = "group",
    shape = "lane",
    lab = NULL,
    pointSize = 3,
    legendPosition = "right"
  ) +
  guides(
    color = guide_legend("Cell type"),
    shape = guide_legend("Sequencing batch")
  )

# Save data
saveRDS(plot_pca, here::here("output/pca/plot_pca.rds"))