library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

path <- "dev/project-area/data/crohns_scrnaseq/scvi_tools_output/"
rds_path <- file.path(path, "Integrated_05_label/query_concat_curated.Rds")
seurat_obj <- readRDS(rds_path)

seurat_obj[["RNA"]] <- seurat_obj[["originalexp"]]
DefaultAssay(seurat_obj) <- "RNA"

# Compute per-cluster averages
gene_order <- top5 %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

cluster_means <- AverageExpression(
  seurat_obj,
  features = gene_order,
  assays = "RNA",
  slot = "data"
)$RNA

cluster_means@Dimnames[[2]] <- gsub("^g", "", cluster_means@Dimnames[[2]])

# Z-score per gene
cluster_means_z <- t(scale(t(cluster_means)))
summary(as.vector(cluster_means_z))


# Define color scale
low <- -2
high <- 4
neg_breaks <- seq(low, 0, length.out = 50)
pos_breaks <- seq(0, high, length.out = 51)
breaks <- c(neg_breaks, pos_breaks[-1])

cols <- c(
  colorRampPalette(c("darkblue", "white"))(length(neg_breaks)),
  colorRampPalette(c("white", "darkred"))(length(pos_breaks) - 1)
)

### Plot heatmap
marker_heatmap <- function(data, title) {
  pdf(file.path(path, paste0(title, "_markers_heatmap.pdf")), width = 7, height = nrow(data)*0.5)

  ht <- pheatmap(
    data,
    color = cols,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = expression("Top Marker Genes per Cluster"),
    angle_col = "90",
    silent = TRUE,
    fontsize = 9
  )
  
  # Draw heatmap 
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = unit(0.9, "npc"), height = unit(0.95, "npc"), x = 0.48))
  grid::grid.draw(ht$gtable)
  
  grid::grid.text(
    "Cluster",
    y = unit(-1, "lines"), x = 0.4,
    gp = grid::gpar(cex = 1)
  )
  
  grid::grid.text(
    "Marker gene",
    x = 1, y = 0.5,
    rot = 90,
    gp = grid::gpar(cex = 1)
  )
  
  grid::grid.text(
    expression("Gene expression (TP10K+1) Z-score"),
    x = 1.025, y = 0.925, rot = 90,
    gp = grid::gpar(cex = 1)
  )
  
  dev.off()
}

marker_heatmap(cluster_means_z, "cluster")


###Top 5 marker genes per curated cell type###
Idents(seurat_obj) <- "curated"

markers_curated <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5 <- markers_curated %>%
  filter(!grepl("^(MT-|RPS|RPL)", gene, ignore.case = TRUE)) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

gene_order <- top5 %>%
  left_join(category_map, by = c("cluster" = "curated")) %>%
  arrange(category, cluster, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

cluster_means <- AverageExpression(
  seurat_obj,
  features = gene_order,
  assays = "RNA",
  group.by = "curated"
)$RNA

ordered_celltypes <- category_map %>%
  arrange(category, curated) %>%
  pull(curated)

cluster_means <- cluster_means[, ordered_celltypes, drop = FALSE]
cluster_means_z <- t(scale(t(cluster_means)))

category_map_clean <- category_map %>%
  filter(curated %in% colnames(cluster_means_z)) %>%
  distinct(curated, category)

ordered_celltypes <- category_map_clean %>%
  arrange(category, curated) %>%
  pull(curated)

cluster_means_z <- cluster_means_z[, ordered_celltypes, drop = FALSE]

summary(as.vector(cluster_means_z))

facet_marker_heatmap <- function(data, meta, title) {
  pdf(file.path(path, paste0(title, "_facet_markers_heatmap.pdf")), width = 10, height = 20)
  
  # Align metadata to data order — keep existing order
  meta_ordered <- meta %>%
    filter(curated %in% colnames(data)) %>%
    distinct(curated, category) %>%
    mutate(curated = factor(curated, levels = colnames(data)))
  
  meta_ordered <- meta_ordered[order(meta_ordered$curated), ]
  
  # Subset the matrix exactly as given (no reordering)
  submat <- data[, colnames(data), drop = FALSE]
  
  # Compute category boundaries according to *existing* column order
  category_boundaries <- meta_ordered %>%
    group_by(category) %>%
    summarise(boundary = n(), .groups = "drop") %>%
    mutate(boundary = cumsum(boundary))
  
  ht <- pheatmap(
    submat,
    color = cols,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = expression("Top Marker Genes per Cell Type"),
    angle_col = "90",
    fontsize = 9,
    silent = TRUE
  )
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = unit(0.9, "npc"), height = unit(0.95, "npc"), x = 0.48))
  grid::grid.draw(ht$gtable)
  
  total_cols <- ncol(submat)
  
  grid::grid.text("Cell type", y = unit(-1, "lines"), x = 0.5, gp = grid::gpar(cex = 1))
  grid::grid.text("Marker gene", x = unit(-3, "lines"), y = 0.5, rot = 90, gp = grid::gpar(cex = 1))
  grid::grid.text(
    expression("Gene expression (TP10K + 1) Z-score)"),
    x = 1.03, y = 0.925, rot = 90,
    gp = grid::gpar(cex = 0.9)
  )
  
  dev.off()
}

facet_marker_heatmap(cluster_means_z, category_map, "curated")