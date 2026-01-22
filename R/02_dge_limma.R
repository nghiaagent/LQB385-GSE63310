# Import packages
library("edgeR")
library("limma")
library("tidyverse")

# Load data
dge <- readRDS(here::here("output/dge.rds"))
design <- readRDS(here::here("output/design.rds"))
contrasts <- readRDS(here::here("output/contrasts.rds"))
contrasts_coefs <- c(1:ncol(contrasts)) %>%
  set_names(colnames(contrasts))
# Process data prior to analysis
dge <- dge %>%
  # Remove lowly expressed genes
  .[filterByExpr(.), , keep.lib.sizes = FALSE] %>%
  # Calculate TMM normalization factors
  calcNormFactors(method = "TMM")

# Perform voom CPM transform
dge_voom <- dge %>%
  voom(design = design)

# Perform linear fit
dge_fit <- dge_voom %>%
  lmFit(design = design) %>%
  contrasts.fit(contrasts = contrasts) %>%
  eBayes()

# Output top tables
dge_toptables <- contrasts_coefs %>%
  map(\(coef) {
    # Get top 10 upregulated genes
    upreg <- topTable(
      dge_fit,
      coef = coef,
      number = Inf
    ) %>%
      filter(adj.P.Val < 0.05, logFC > 0) %>%
      arrange(desc(logFC)) %>%
      top_n(10)

    # Get bottom 10 upregulated genes
    downreg <- topTable(
      dge_fit,
      coef = coef,
      number = Inf
    ) %>%
      filter(adj.P.Val < 0.05, logFC < 0) %>%
      arrange(desc(logFC)) %>%
      top_n(-10)

    # Merge and return
    all <- rbind(upreg, downreg)
    return(all)
  })

# Save data
saveRDS(dge_voom, here::here("output/dge_voom.rds"))
saveRDS(dge_fit, here::here("output/dge_fit.rds"))
saveRDS(dge_toptables, here::here("output/dge_toptables.rds"))

# Write DE gene list
contrasts_coefs %>%
  iwalk(\(coef, name) {
    genes <- topTable(
      dge_fit,
      coef = coef,
      number = Inf
    ) %>%
      filter(adj.P.Val < 0.05) %>%
      .$SYMBOL

    write_lines(genes, file = str_c("output/dge_list/", name, ".txt"))
  })
