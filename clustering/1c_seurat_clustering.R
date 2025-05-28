# Run clustering
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

for (ymin in disps) {
  case <- paste("mean.var.plot_disp", ymin, sep="")
  dr <- paste("pca", case, sep="_")
  out_dir <- file.path(OUTDIR, case)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  seurat_object <- ClustDr(
    seurat_object, 
    dr = dr, 
    k = k, 
    res = res, 
    file = file.path(out_dir, "num_PCs.txt")
  )
}

for (nfeat in nfeats) {
  case <- paste("vst_top", nfeat, sep="")
  dr <- paste("pca", case, sep="_")
  out_dir <- file.path(OUTDIR, case)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  seurat_object <- ClustDr(
    seurat_object, 
    dr = dr, 
    k = k, 
    res = res, 
    file = file.path(out_dir, "num_PCs.txt")
  )
}

saveRDS(seurat_object, file = file.path(OUTDIR, "seurat_object.Rds"))
