here::i_am("R/99_build_dataset.R")

# Run all scripts
source(here::here("R/01_load_gse63310.R"))
source(here::here("R/01_load_msigdb.R"))
source(here::here("R/02_dge_limma.R"))
source(here::here("R/03_camera.R"))
source(here::here("R/03_correlation.R"))
source(here::here("R/04_pca.R"))

# Rename objects
GSE_logCPM <- dge_voom
GSE_fitted <- dge_fit
GSE_camera <- camera_all
GSE_pca <- dge_pca

# Save data
save(
  GSE_logCPM,
  GSE_fitted,
  GSE_camera,
  GSE_pca,
  file = "analysed_GSE63310.Rdata"
)
