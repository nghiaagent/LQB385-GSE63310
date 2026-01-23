here::i_am("R/04_heatmap.R")

# Import packages
library("circlize")
library("cluster")
library("ComplexHeatmap")
library("limma")
library("tidyverse")
library("viridis")

# Load data
dge_voom <- readRDS(here::here("output/dge_voom.rds"))
dge_fit <- readRDS(here::here("output/dge_fit.rds"))
dge <- readRDS(here::here("output/dge.rds"))

# Set color scheme, breaks, number of gene clusters
## For gene expression (rlog z-log counts per million)
col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)
n_kmeans <- 3

# Get top 200 DEGs by limma F-test
## Also k-means cluster to split this into 3 clusters
top_genes <- topTable(
  fit = dge_fit,
  coef = NULL,
  number = 200,
  sort.by = "F"
) %>%
  .$ENTREZID

# Get expression data, extract top 1000 DEGs, Z normalise for heatmap
heat_expression <- dge_voom$E[top_genes, ] %>%
  t() %>%
  scale() %>%
  t()

## Get kmeans
heat_kmeans <- heat_expression %>%
  pam(k = n_kmeans)

heat_kmeans$clustering <- heat_kmeans$clustering %>%
  factor(
    levels = c(1:n_kmeans),
    labels = str_c("Cluster ", c(1:n_kmeans))
  )

# Get column annotation data
heat_column_annotation <- dge$samples %>%
  dplyr::select(group, lane) %>%
  mutate(
    group = group %>%
      factor(
        levels = c(
          "Basal",
          "LP",
          "ML"
        ),
        labels = c(
          "Basal",
          "Luminal progenitor",
          "Mature luminal"
        )
      )
  )

# Define column annotation object
heat_column_annotation <- HeatmapAnnotation(
  df = heat_column_annotation,
  which = "col",
  annotation_height = 0.6,
  annotation_legend_param = list(
    `group` = list(
      nrow = 3,
      title = "Cell type",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    `lane` = list(
      nrow = 3,
      title = "Sequencing batch",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    )
  )
)

# Create heatmap object
heat <- Heatmap(
  heat_expression,
  name = "Gene\nZ-\nscore",
  col = colorRamp2(breaks, col),
  border = FALSE,

  # parameters for the colour-bar that represents gradient of expression
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    legend_width = unit(8, "cm"),
    legend_height = unit(5.0, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontface = "bold")
  ),

  # row (gene) parameters
  row_split = heat_kmeans$clustering,
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 90,
  show_row_names = FALSE,

  # column (sample) parameters
  column_title = NULL,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_column_names = TRUE,

  # specify top and bottom annotations
  top_annotation = heat_column_annotation
)

# Export heatmap
png(
  file = here::here("output/heatmap/heatmap_log_counts_per_million.png"),
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

plot(heat)

dev.off()
