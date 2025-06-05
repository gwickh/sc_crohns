source(file.path("sc_crohns/clustering", ".Rprofile"))

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

# Function to identify number of informative PCs --------------------------
# based on standard deviation cutoff and cluster with Louvain at varying on
# resolution and kNN
Clust_dr <- function(
    seurat_object, 
    dr, 
    pcs = 50, 
    non.rand.sd.frac = 0.5, 
    K, 
    res, 
    file
) {
  # Compute minimum std a PC must have to be considered “non-random” assuming
  # last 10 PCs approximate random noise to exclude SVD artifacts
  stdevs <- Stdev(seurat_object[[dr]])
  last10_idx <- max(1, pcs - 9):pcs
  mean_stdev_lastPCs <- mean(stdevs[last10_idx])
  min_stdev <- (1 + non.rand.sd.frac) * mean_stdev_lastPCs
  
  # Loops through PCs to determine nPCs, the number of informative PCs that 
  # explains > 0.05 variance with respect to the next
  nPCs <- 0
  for (i in seq_along(stdevs)) {
    if (stdevs[i] < min_stdev) break
    nPCs <- i
  }
  write(nPCs, file = file)
  
  # For each k build a NN graph using the top nPCs components
  for (k in K) {
    seurat_object <- FindNeighbors(
      seurat_object, 
      k.param = k, 
      reduction = dr, 
      dims = 1:nPCs
    )
    # For each resolution, run Louvain clustering
    for (r in res) {
      cl_ident <- paste("clusters", dr, "k", k, "res", r, sep="_")
      seurat_object <- FindClusters(seurat_object, resolution = r)
      seurat_object[[cl_ident]] <- Idents(seurat_object)
    }
  }
  
  return(seurat_object)
}

# Run clustering ----------------------------------------------------------
# Create directory to store PCA stats
CLUSTERING_OUTPUT_PATH <- file.path(SEURAT_OBJECT_PATH, "clustering_stats")
dir.create(CLUSTERING_OUTPUT_PATH, showWarnings = FALSE)

# Define params
disps <- seq(0.5, 1.5, 0.5)
n_features <- c(1000, 2000, 5000)
neighbors <- seq(10, 50, 5)
res <- seq(0.5, 1.2, 0.1)

for (disp in disps) {
  filename <- paste("mean.var.plot_disp", disp, sep="_")
  dr_name <- paste("pca", filename, sep="_")
  out_dir <- file.path(SEURAT_OBJECT_PATH, dr_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  seurat_object <- Clust_dr(
    seurat_object, 
    dr = dr_name, 
    K = neighbors, 
    res = res, 
    file = file.path(
      CLUSTERING_OUTPUT_PATH, 
      paste(dr_name, "num_PCs.txt", sep = "_")
    )
  )
}

for (n in n_features) {
  filename <- paste("vst_top", n, sep="")
  dr_name <- paste("pca", filename, sep="_")
  out_dir <- file.path(SEURAT_OBJECT_PATH, dr_name)
  seurat_object <- Clust_dr(
    seurat_object, 
    dr = dr_name, 
    K = neighbors, 
    res = res, 
    file = file.path(
      CLUSTERING_OUTPUT_PATH, 
      paste(dr_name, "num_PCs.txt", sep = "_")
    )
  )
}

saveRDS(seurat_object, file = file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
