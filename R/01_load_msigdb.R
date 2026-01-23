here::i_am("R/01_load_msigdb.R")

# Import packages
library("GSEABase")
library("tidyverse")

# Perform process if the appropriate RDS file is missing
if (!file.exists(here::here("input/MSigDB/msigdb_mh.rds"))) {
  # Convert to RDS for camera
  msigdb_mh <- getGmt(
    here::here("input/MSigDB/mh.all.v2025.1.Mm.entrez.gmt"),
    geneIdType = EntrezIdentifier()
  )

  # Save data
  saveRDS(msigdb_mh, here::here("input/MSigDB/msigdb_mh.rds"))
}
