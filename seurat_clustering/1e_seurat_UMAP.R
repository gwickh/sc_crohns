source(file.path("sc_crohns/seurat_clustering", ".Rprofile"))
CLUSTERING_OUTPUT_PATH <- file.path(SEURAT_OBJECT_PATH, "clustering_stats")

UMAP_PATH <- file.path(SEURAT_OBJECT_PATH, "UMAP_plots")
dir.create(UMAP_PATH, showWarnings = FALSE, recursive = TRUE)

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

# ----------------- Function: Run and Plot UMAP -----------------
RunUMAPandPlot <- function(
    seurat_object, 
    disp, 
    k, 
    res, 
    cases
) {
  # Define paths and names
  dr <- paste0("pca_", cases)
  umap_dr <- paste0("umap_", dr)
  pc_file <- file.path(CLUSTERING_OUTPUT_PATH, paste(dr, "num_PCs.txt", sep = "_"))
  outpath <- file.path(UMAP_PATH, dr)
  dir.create(outpath, showWarnings = FALSE)
  
  # Get number of PCs
  nPCs <- drop(as.matrix(read.table(pc_file)))
  
  # Run UMAP if not already present
  if (!(umap_dr %in% names(seurat_object@reductions))) {
    seurat_object <- RunUMAP(
      seurat_object, 
      reduction = dr, 
      dims = 1:nPCs, 
      reduction.name = umap_dr, 
      reduction.key = paste0("UMAP_", dr, "_")
    )
  }
  
  # Internal function for saving plots
  ggsave_plot <- function(plot_obj, filename, width = 7.4, height = 6) {
    plot_obj <- plot_obj + theme_minimal() + theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(filename, plot = plot_obj, width = width, height = height)
  }
  
  # Plot DimPlots by metadata groups
  dimplot_configs <- list(
    list(group.by = "orig.ident", file_suffix = "samples"),
    list(group.by = "Phase", file_suffix = "ccScore")
  )
  
  for (cfg in dimplot_configs) {
    if (cfg$group.by %in% colnames(seurat_object@meta.data)) {
      p <- DimPlot(seurat_object, group.by = cfg$group.by, reduction = umap_dr)
      ggsave_plot(
        p, 
        file.path(outpath, paste0("UMAP_plot_", cfg$file_suffix, ".png"))
      )
    }
  }
  
  # Plot FeaturePlot and density of nCount_RNA
  if ("nCount_RNA" %in% colnames(seurat_object@meta.data)) {
    fp <- FeaturePlot(seurat_object, features = "nCount_RNA", reduction = umap_dr, cols = c("blue", "red"))
    ggsave_plot(
      fp, 
      file.path(outpath, paste("UMAP_plot_UMIcount.png", sep = "_"))
    )
    
    dp <- ggplot(seurat_object@meta.data, aes(x = nCount_RNA)) +
      theme_classic() + geom_density() + ylab("density")
    ggsave_plot(
      dp, 
      file.path(outpath, "density_plot_UMIcount.png"), 
      width = 3.5, 
      height = 3.5
    )
  }
  
  # Plot clustering results
  for (kk in k) {
    for (r in res) {
      cl_ident <- paste("clusters", dr, "k", kk, "res", r, sep = "_")
      print(paste("Looking for cluster column:", cl_ident))
      print(cl_ident %in% colnames(seurat_object@meta.data))
      
      if (cl_ident %in% colnames(seurat_object@meta.data)) {
        p <- DimPlot(seurat_object, group.by = cl_ident, reduction = umap_dr)
        ggsave_plot(p, file.path(outpath, paste0("UMAP_plot_", cl_ident, ".png")))
      } else {
        warning(paste("Column not found in metadata:", cl_ident))
      }
    }
  }
  
  return(seurat_object)
}



# ----------------- Execute -----------------
disp_cases <- paste0("mean.var.plot_disp_", disps)
vst_cases  <- paste0("vst_top_", n_features)
cases <- c(disp_cases, vst_cases)

theme_set(theme_cowplot())

for (case in cases) {
  seurat_object <- RunUMAPandPlot(
    seurat_object = seurat_object,
    disp = disp, 
    k = neighbors, 
    res = res, 
    cases = case
  )
}

saveRDS(seurat_object, file = file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))

