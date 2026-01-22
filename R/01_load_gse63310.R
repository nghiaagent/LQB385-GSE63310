# Import packages
library("AnnotationDbi")
library("edgeR")
library("here")
library("limma")
library("Mus.musculus")
library("RNAseq123")
library("R.utils")
library("tidyverse")

# Check if files are all DL'd
files <- c(
  "GSE63310_RAW.tar",
  "GSM1545535_10_6_5_11.txt",
  "GSM1545536_9_6_5_11.txt",
  "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt",
  "GSM1545540_JMS8-3.txt",
  "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt",
  "GSM1545544_JMS9-P7c.txt",
  "GSM1545545_JMS9-P8c.txt"
)

files_sel <- c(
  "GSM1545535_10_6_5_11.txt",
  "GSM1545536_9_6_5_11.txt",
  "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt",
  "GSM1545540_JMS8-3.txt",
  "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt",
  "GSM1545544_JMS9-P7c.txt",
  "GSM1545545_JMS9-P8c.txt"
)

if (all(file.exists(here::here("input/GSE63310/", files)))) {
  message("All files already downloaded.")
  is_downloaded <- TRUE
} else {
  message("Some files are missing. Downloading...")
  is_downloaded <- FALSE
}

# Download GSE63310 data only if it is not already downloaded
if (!is_downloaded) {
  # Download GSE data and unzip
  download.file(
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file",
    destfile = here::here("input/GSE63310/GSE63310_RAW.tar"),
    mode = "wb"
  )

  untar(
    here::here("input/GSE63310/GSE63310_RAW.tar"),
    exdir = here::here("input/GSE63310")
  )

  # Unzip samples
  files_sel %>%
    str_c("input/GSE63310/", ., ".gz", sep = "") %>%
    walk(gunzip)
}

# Create DGEList
dge <- readDGE(
  files = here::here("input/GSE63310", files_sel),
  columns = c(1, 3)
)

# Add samplesheet
## Extract column names, replace original column names with simpler names
colnames(dge) <- colnames(dge) %>%
  str_replace(".*/GSM\\d+_", "")

## Create samplesheet
dge$samples$group <- factor(c(
  "LP",
  "ML",
  "Basal",
  "Basal",
  "ML",
  "LP",
  "Basal",
  "ML",
  "LP"
))

dge$samples$lane <- factor(rep(c("L004", "L006", "L008"), c(3, 4, 2)))

# Add gene info
geneinfo <- AnnotationDbi::select(
  Mus.musculus,
  keys = rownames(dge),
  columns = c("SYMBOL", "TXCHROM"),
  keytype = "ENTREZID"
) %>%
  distinct(ENTREZID, .keep_all = TRUE)

dge$genes <- data.frame(ENTREZID = rownames(dge))
dge$genes <- left_join(dge$genes, geneinfo, by = "ENTREZID")

# Create design matrix
design <- model.matrix(~ 0 + group + lane, data = dge$samples)
colnames(design) <- colnames(design) %>%
  str_replace("group", "")

# Create contrast matrix
contrasts <- makeContrasts(
  lp_vs_basal = LP - Basal,
  ml_vs_basal = ML - Basal,
  ml_vs_lp = ML - LP,
  levels = colnames(design)
)

# Save data
saveRDS(dge, here::here("output/dge.rds"))
saveRDS(design, here::here("output/design.rds"))
saveRDS(contrasts, here::here("output/contrasts.rds"))
