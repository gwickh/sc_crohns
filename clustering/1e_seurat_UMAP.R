source(file.path("sc_crohns/clustering", ".Rprofile"))

# Function to get significant PCs
GetNPCs <- function(in_dir) {
  pc_file <- file.path(in_dir, "num_PCs.txt")
  nPCs <- drop(as.matrix(read.table(pc_file)))
  
  nPCs
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

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_LOC, "seurat_object.Rds"))
}

theme_set(theme_cowplot())

  for (ymin in disps) {
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

for (nfeat in nfeats)

saveRDS(seurat_object, file = file.path(OUTDIR, "seurat_object.Rds"))

