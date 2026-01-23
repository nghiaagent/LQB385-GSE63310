here::i_am("R/04_volcano.R")

# Import packages
library("EnhancedVolcano")
library("limma")
library("tidyverse")

# Load data
dge_fit <- readRDS(here::here("output/dge_fit.rds"))
contrasts <- readRDS(here::here("output/contrasts.rds"))
contrasts_coefs <- c(1:ncol(contrasts)) %>%
  set_names(colnames(contrasts))

# Draw volcano plots for each contrast
plots_volcano <- contrasts_coefs %>%
  imap(\(coef, name) {
    toptable <- topTable(
      fit = dge_fit,
      coef = coef,
      number = Inf
    )

    volcano_plot <- EnhancedVolcano(
      toptable = toptable,
      lab = toptable$SYMBOL,
      x = "logFC",
      y = "adj.P.Val",
      ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
      pCutoff = 0.05,
      FCcutoff = 1,
      pointSize = 1.0,
      labSize = 4.0,
      title = str_c("Volcano Plot - ", name),
      subtitle = "",
      caption = NULL,
      legendPosition = "right",
      legendLabSize = 12,
      legendIconSize = 4.0,
      colAlpha = 1.0,
      gridlines.major = FALSE,
      gridlines.minor = FALSE
    )

    return(volcano_plot)
  })

# Save data
saveRDS(
  object = plots_volcano,
  file = here::here("output/volcano/plots_volcano.rds")
)
