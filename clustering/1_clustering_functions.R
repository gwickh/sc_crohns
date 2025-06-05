DrPlot <- function(
    seurat_object, 
    dr = "pca", 
    dr_plot = "pca", 
    pcs = 50, 
    pval = 1e-05, 
    non.rand.sd.frac = 0.5, 
    k, 
    res, 
    OUTDIR, 
    dr_type = "PCA"
) {
  ggsave_plot <- function(plot_obj, filename, width = 7.4, height = 6) {
    ggsave(filename, plot = plot_obj, width = width, height = height)
  }
  dimplot_configs <- list(
    list(group.by = "orig.ident", file_suffix = "samples"),
    list(group.by = "Phase", file_suffix = "ccScore")
  )
  for (cfg in dimplot_configs) {
    if (cfg$group.by %in% colnames(seurat_object@meta.data)) {
      p <- DimPlot(seurat_object, group.by = cfg$group.by, reduction = dr_plot)
      ggsave_plot(p, file.path(OUTDIR, paste0(dr_type, "_plot_", cfg$file_suffix, ".pdf")))
    }
    if ("nCount_RNA" %in% colnames(seurat_object@meta.data)) {
      fp <- FeaturePlot(seurat_object, features = "nCount_RNA", reduction = dr_plot, cols = c("blue", "red"))
      ggsave_plot(fp, file.path(OUTDIR, paste0(dr_type, "_plot_UMIcount.pdf")))
      
      dp <- ggplot(seurat_object@meta.data, aes(x = nCount_RNA)) +
        theme_classic() + geom_density() + ylab("density")
      ggsave_plot(dp, file.path(OUTDIR, "density_plot_UMIcount.pdf"), width = 3.5, height = 3.5)
    }
  }
  for (kk in k) {
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep = "")
      if (cl_ident %in% colnames(seurat_object@meta.data)) {
        p <- DimPlot(seurat_object, group.by = cl_ident, reduction = dr_plot)
        ggsave_plot(p, file.path(OUTDIR, paste0(dr_type, "_plot_", cl_ident, ".pdf")))
      }
    }
  }
}
run_UMAP <- function() {
  case <- paste("mean.var.plot_disp", ymin, sep = "")
  dr <- paste("pca", case, sep = "_")
  tsne_dr <- paste0("tsne_", dr)
  umap_dr <- paste0("umap_", dr)
  in_dir <- file.path(SEURAT_OBJECT_LOC, case)
  WRITEOUT <- file.path(OUTDIR, dr)
  dir.create(WRITEOUT, showWarnings = FALSE)
  nPCs <- GetNPCs(in_dir)
  if (!(tsne_dr %in% names(seurat_object@reductions)))
    seurat_object <- RunTSNE(seurat_object, reduction = dr, dims = 1:nPCs, reduction.name = tsne_dr, reduction.key = paste0("tSNE_", dr, "_"))
  if (!(umap_dr %in% names(seurat_object@reductions)))
    seurat_object <- RunUMAP(seurat_object, reduction = dr, dims = 1:nPCs, reduction.name = umap_dr, reduction.key = paste0("UMAP_", dr, "_"))
  DrPlot(seurat_object, dr = dr, dr_plot = dr, k = k, res = res, OUTDIR = WRITEOUT, dr_type = "PCA")	
  DrPlot(seurat_object, dr = dr, dr_plot = tsne_dr, k = k, res = res, OUTDIR = WRITEOUT, dr_type = "TSNE")
  DrPlot(seurat_object, dr = dr, dr_plot = umap_dr, k = k, res = res, OUTDIR = WRITEOUT, dr_type = "UMAP")
}