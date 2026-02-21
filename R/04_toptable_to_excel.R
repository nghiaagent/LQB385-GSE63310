here::i_am("r/04_toptable_to_excel.R")

# Import packages
library("limma")
library("tidyverse")
library("openxlsx")

# Load data
dge_fit <- readRDS(here::here("output/dge_fit.rds"))
contrasts <- readRDS(here::here("output/contrasts.rds"))
contrasts_coefs <- c(1:ncol(contrasts)) %>%
  set_names(colnames(contrasts))

# Build list of DEGs
toptables <- contrasts_coefs %>%
  map(\(coef) {
    # Get top upregulated genes
    upreg <- topTable(
      dge_fit,
      coef = coef,
      number = Inf
    ) %>%
      filter(adj.P.Val < 0.05, logFC > 0) %>%
      arrange(desc(logFC))

    # Get top 10 downregulated genes
    downreg <- topTable(
      dge_fit,
      coef = coef,
      number = Inf
    ) %>%
      filter(adj.P.Val < 0.05, logFC < 0) %>%
      arrange(desc(logFC))

    # Return merged table
    combined <- list(
      "upregulated" = upreg,
      "downregulated" = downreg
    )
    return(combined)
  }) %>%
  unlist(recursive = FALSE)

# Export to excel
write.xlsx(
  toptables,
  file = "output/dge_list/GSE63310_degs.xlsx",
  asTable = TRUE
)
