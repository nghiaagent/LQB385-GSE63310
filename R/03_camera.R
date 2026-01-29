here::i_am("R/03_camera.R")

# Import packages
library("edgeR")
library("limma")
library("tidyverse")

# Load data
dge_voom <- readRDS(here::here("output/dge_voom.rds"))
msigdb_mh <- readRDS(here::here("input/MSigDB/msigdb_mh.rds"))
design <- readRDS(here::here("output/design.rds"))
contrasts <- readRDS(here::here("output/contrasts.rds"))

contrasts_coefs <- c(1:ncol(contrasts)) %>%
  set_names(colnames(contrasts))

# Prepare gene sets
msigdb_mh_ready <- msigdb_mh %>%
  geneIds() %>%
  ids2indices(identifiers = rownames(dge_voom))

# Perform camera
camera_all <- contrasts_coefs %>%
  map(\(coef) {
    camera_res <- camera(
      y = dge_voom$E,
      index = msigdb_mh_ready,
      design = design,
      contrast = contrasts[, coef]
    )

    return(camera_res)
  })

# Save data
saveRDS(camera_all, here::here("output/camera_all.rds"))
