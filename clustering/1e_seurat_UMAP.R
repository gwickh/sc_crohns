# Set static variables
SCRIPT_DIR <- "sc_crohns/clustering"
SEURAT_OBJECT_LOC <- "project-area/data/crohns_scrnaseq/clustering_output"
OUTDIR <- SEURAT_OBJECT_LOC

# Source .Rprofile from sc_crohns/clustering/
source(file.path(SCRIPT_DIR, ".Rprofile"))

# Function to get significant PCs
GetNPCs <- function(in_dir) {
  pc_file <- file.path(in_dir, "num_PCs.txt")
  nPCs <- drop(as.matrix(read.table(pc_file)))
  
  nPCs
}

DrPlot <- function(seurat_object, dr = "pca", dr_plot = "pca", pcs = 50, pval = 1e-05, non.rand.sd.frac = 0.5, k, res, OUTDIR, dr_type = "PCA") {
  
  pdf(file.path(OUTDIR, paste0(dr_type, "_plot_samples.pdf")), width = 3.7*2, height = 3*2)
  print(DimPlot(seurat_object, group.by = "sample.name", reduction = dr_plot))
  dev.off()
  
  pdf(file.path(OUTDIR, paste0(dr_type, "_plot_ccScore.pdf")), width = 3.7*2, height = 3*2)
  print(DimPlot(seurat_object, group.by = "Phase", reduction = dr_plot))
  dev.off()
  
  pdf(file.path(OUTDIR, paste0(dr_type, "_plot_UMIcount.pdf")), width = 3.7*2, height = 3*2)
  print(FeaturePlot(seurat_object, features = "nCount_RNA", reduction = dr_plot, cols = c("blue","red")))
  dev.off()
  
  g <- ggplot(data = seurat_object@meta.data, aes(x = nCount_RNA)) + theme_classic() + geom_density() + ylab("density")
  pdf(file.path(OUTDIR, "density_plot_UMIcount.pdf"), width = 3.5, height = 3.5)
  print(g)
  dev.off()
  
  for (kk in k) {
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
      pdf(file.path(OUTDIR, paste0(dr_type, "_plot_", cl_ident, ".pdf")), width = 3.7*2, height = 3*2)
      print(DimPlot(seurat_object, group.by = cl_ident, reduction = dr_plot))
      dev.off()
    }
  }
}

clustering_params <- list(
  disps = 0.5,
  nfeats = 5000,
  k = 30,
  res = c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2)
)

disps <- clustering_params$disps
nfeats <- clustering_params$nfeats
k <- clustering_params$k
res <- clustering_params$res

options(future.globals.maxSize = 60 * 1024^3)
max_workers <- parallelly::availableCores()
cores <- min(max_workers, detectCores() - 1)
plan("multisession", workers = cores, future.seed = TRUE)
seurat_object <- readRDS(file.path(SEURAT_OBJECT_LOC, "clustered_object.Rds"))

theme_set(theme_cowplot())

for (ymin in disps) {
  case <- paste("mean.var.plot_disp", ymin, sep = "")
  dr <- paste("pca", case, sep = "_")
  tsne_dr <- paste0("tsne_", dr)
  umap_dr <- paste0("umap_", dr)
  in_dir <- file.path(SEURAT_OBJECT_LOC, case)
  OUTDIR <- file.path(OUTDIR, dr)
  dir.create(OUTDIR, showWarnings = FALSE)
  nPCs <- GetNPCs(in_dir)
  if (!(tsne_dr %in% names(seurat_object@reductions)))
    seurat_object <- RunTSNE(seurat_object, reduction = dr, dims = 1:nPCs, reduction.name = tsne_dr, reduction.key = paste0("tSNE_", dr, "_"))
  if (!(umap_dr %in% names(seurat_object@reductions)))
    seurat_object <- RunUMAP(seurat_object, reduction = dr, dims = 1:nPCs, reduction.name = umap_dr, reduction.key = paste0("UMAP_", dr, "_"))
  DrPlot(seurat_object, dr = dr, dr_plot = dr, k = k, res = res, OUTDIR = OUTDIR, dr_type = "PCA")	
  DrPlot(seurat_object, dr = dr, dr_plot = tsne_dr, k = k, res = res, OUTDIR = OUTDIR, dr_type = "TSNE")
  DrPlot(seurat_object, dr = dr, dr_plot = umap_dr, k = k, res = res, OUTDIR = OUTDIR, dr_type = "UMAP")	
}

for (nfeat in nfeats) {
  case <- paste("vst_top", nfeat, sep = "")
  dr <- paste("pca", case, sep = "_")
  tsne_dr <- paste0("tsne_", dr)
  umap_dr <- paste0("umap_", dr)
  in_dir <- file.path(SEURAT_OBJECT_LOC, case)
  OUTDIR <- file.path(OUTDIR, dr)
  dir.create(OUTDIR, showWarnings = FALSE)
  nPCs <- GetNPCs(in_dir)
  if (!(tsne_dr %in% names(seurat_object@reductions)))
    seurat_object <- RunTSNE(seurat_object, reduction = dr, dims = 1:nPCs, reduction.name = tsne_dr, reduction.key = tsne_dr)
  if (!(umap_dr %in% names(seurat_object@reductions)))
    seurat_object <- RunUMAP(seurat_object, reduction = dr, dims = 1:nPCs, reduction.name = umap_dr, reduction.key = umap_dr)
  DrPlot(seurat_object, dr = dr, dr_plot = dr, k = k, res = res, OUTDIR = OUTDIR, dr_type = "PCA")	
  DrPlot(seurat_object, dr = dr, dr_plot = tsne_dr, k = k, res = res, OUTDIR = OUTDIR, dr_type = "TSNE")
  DrPlot(seurat_object, dr = dr, dr_plot = umap_dr, k = k, res = res, OUTDIR = OUTDIR, dr_type = "UMAP")
}

saveRDS(seurat_object, file = file.path(OUTDIR, "UMAP_projected_object.Rds"))

