library(Seurat)
library(DoubletFinder)
library(dplyr)
library(patchwork)


path <- "dev/project-area/data/crohns_scrnaseq/scvi_tools_output/"
rds_path <- file.path(path, "Integrated_05_label/query_concat_curated.Rds")
seurat_obj <- readRDS(rds_path)


# Ensure correct assay is active
seurat_obj[["RNA"]] <- seurat_obj[["originalexp"]]
DefaultAssay(seurat_obj) <- "RNA"


seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:24, k.param = 30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.9)

category_map <- seurat_obj@meta.data %>%
  group_by(curated, category) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n)) %>%
  slice_max(n, n = 1) %>%  # dominant category per curated type
  select(curated, category)


df <- seurat_obj@meta.data %>%
  select(curated, category, originalexp_snn_res.0.9, Diagnosis) %>%
  mutate(cluster = as.factor(originalexp_snn_res.0.9),
         curated = as.factor(curated))

prop_df <- df %>%
  group_by(Diagnosis, cluster, curated) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Diagnosis, curated) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

categorised_prop_df <- prop_df %>%
  left_join(category_map, by = "curated")

# calculate doublets 
doublet_by_cluster <- seurat_doublet@meta.data %>%
  group_by(!!sym("Diagnosis"), !!sym("originalexp_snn_res.0.9")) %>%
  summarise(
    total_cells = n(),
    doublets = sum(DoubletFinder == "Doublet"),
    doublet_frac_cluster = doublets / total_cells,
    .groups = "drop"
  )
doublet_by_celltype <- seurat_doublet@meta.data %>%
  group_by(!!sym("Diagnosis"), !!sym("curated")) %>%
  summarise(
    total_cells = n(),
    doublets = sum(DoubletFinder == "Doublet"),
    doublet_frac_celltype = doublets / total_cells,
    .groups = "drop"
  )

prop_df <- left_join(
  prop_df,
  doublet_by_cluster,
  by = c("Diagnosis", "cluster" = "originalexp_snn_res.0.9")
)


prop_df <- left_join(
  prop_df,
  doublet_by_celltype,
  by = c("Diagnosis", "curated")
)
categorised_prop_df <- prop_df %>%
  left_join(category_map, by = "curated")

# generate heatmap
heatmap <- function(data, title) {
  data$cluster <- factor(data$cluster, levels = sort(unique(as.numeric(as.character(data$cluster)))))
  main <-ggplot(data, aes(x = curated, y = cluster, fill = prop)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "#FFFDD0",   # cream for low values
      high = "darkblue", # navy for high values
      name = "Proportion",
      limits = c(0, 1)
    ) +
    facet_wrap(
      ~ category,
      ncol = 9,
      scales = "free_x",
      space = "free_x",
      strip.position = "top"  
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "#FFFDD0", colour = NA),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(
        size = 10,
        angle = 90,
        hjust = 0
      ),
      strip.clip = "off"
    ) +
    coord_cartesian(clip = "off") +  # ✅ prevents rotated labels from being cropped
    labs(
      x = "Curated Cell Type",
      y = "Cluster (pca_vst_top5000_k30_res0.9)",
      title = title
    )
}



# --- Apply to subsets ---
crohns_categorised_prop_df <- subset(categorised_prop_df, Diagnosis == "Crohn's Disease")
normal_categorised_prop_df <- subset(categorised_prop_df, Diagnosis == "Normal")

pdf(file.path(path, "cluster_celltype_heatmap.pdf"), width = 14, height = 8)
heatmap(normal_categorised_prop_df, "Normal") + heatmap(crohns_categorised_prop_df, "Crohn's Disease")
dev.off()

categorised_prop_df$cluster <- factor(
  categorised_prop_df$cluster,
  levels = sort(unique(as.numeric(as.character(categorised_prop_df$cluster))))
)

# Plot single-column heatmap
pdf(file.path(path, "cluster_doublets_heatmap.pdf"), width = 10, height = 3)
ggplot(categorised_prop_df, aes(y = Diagnosis, x = curated, fill = doublet_frac_cluster)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "#FFFDD0",  # cream
    high = "darkblue",    # blue
    name = "Doublet Fraction",
    limits = c(0, 1)
  ) +
  facet_wrap(
    ~ category,
    ncol = 9,
    scales = "free_x",
    space = "free_x",
    strip.position = "top"  
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#FFFDD0", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(
      size = 10,
      angle = 90,
      hjust = 0
    ),
    strip.clip = "off"
  )
dev.off()
