source(file.path("sc_crohns/clustering", ".Rprofile"))

# Function to identify number of informative PCs --------------------------
# based on standard deviation cutoff and cluster with Louvain at varying on
# resolution and kNN
Clust_dr <- function(
    seurat_object, 
    dr, 
    pcs = 50, 
    pval = 1e-05, 
    non.rand.sd.frac = 0.5, 
    k, 
    res, 
    file
) {
  # Compute minimum std a PC must have to be considered “non-random” assuming
  # last 10 PCs approximate random noise to exclude SVD artifacts
  min_stdev <- (1+non.rand.sd.frac)*mean(seurat_object@reductions[[dr]]@stdev[(pcs-10):pcs])
  
  # Loops through PCs to determine nPCs, the number of informative PCs that 
  # explains > 0.05 variance with respect to the next
  for (nPCs in 0:(length(seurat_object@reductions[[dr]]@stdev)-1)) {
    if (seurat_object@reductions[[dr]]@stdev[nPCs+1] < min_stdev) 
      break
  }
  write(nPCs, file = file)
  
  # For each k build a NN graph using the top nPCs components
  for (kk in k) {
    seurat_object <- FindNeighbors(
      seurat_object, 
      k.param = kk, 
      reduction = dr, 
      dims = 1:nPCs
    )
    # For each resolution, run Louvain clustering
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
      seurat_object <- FindClusters(seurat_object, resolution = r)
      names(seurat_object@meta.data)[names(seurat_object@meta.data) == "seurat_clusters"] <- cl_ident
      seurat_object@active.ident <- as.factor(seurat_object$orig.ident)
    }
  }
  
  return(seurat_object)
}


# Run clustering
disps <- seq(0.5, 1.5, 0.5)
n_features <- c(1000, 2000, 5000)
k <- seq(10, 50, 5)
res <- seq(0.5, 1.2, 0.1)

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

for (disp in disps) {
  case <- paste("mean.var.plot_disp", disp, sep="")
  dr <- paste("pca", case, sep="_")
  out_dir <- file.path(SEURAT_OBJECT_PATH, case)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  seurat_object <- ClustDr(
    seurat_object, 
    dr = dr, 
    k = k, 
    res = res, 
    file = file.path(out_dir, "num_PCs.txt")
  )
}

for (n in n_features) {
  case <- paste("vst_top", n, sep="")
  dr <- paste("pca", case, sep="_")
  out_dir <- file.path(SEURAT_OBJECT_PATH, case)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  seurat_object <- ClustDr(
    seurat_object, 
    dr = dr, 
    k = k, 
    res = res, 
    file = file.path(out_dir, "num_PCs.txt")
  )
}

saveRDS(seurat_object, file = file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
