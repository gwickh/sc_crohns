# Set static variables
SCRIPT_DIR <- "sc_crohns/clustering"
SEURAT_OBJECT_LOC <- "project-area/data/crohns_scrnaseq/clustering_output"
OUTDIR <- SEURAT_OBJECT_LOC

# Source .Rprofile from sc_crohns/clustering/
source(file.path(SCRIPT_DIR, ".Rprofile"))

# Function to identify number of informative PCs based on standard deviation cutoff
# and cluster with Louvain at varying on resolution and kNN
ClustDr <- function(
    seurat_object, 
    dr = "pca", 
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

# Perform clustering based on defined gene list
ClustGenes <- function(seurat_object, genes, genes_list_label = "hvg", k, res) {
  # For each k build a NN graph using the top nPCs components
  for (kk in k) {
    seurat_object <- FindNeighbors(
      seurat_object, 
      k.param = kk, 
      features = genes, 
      dims = NULL
    )
    for (r in res) {
      cl_ident <- paste("clusters_", genes_list_label, "_k", kk, "_res", r, sep="")
      seurat_object <- FindClusters(seurat_object, resolution = r)
      names(seurat_object@meta.data)[names(seurat_object@meta.data) == "seurat_clusters"] <- cl_ident
      seurat_object@active.ident <- as.factor(seurat_object$orig.ident)
    }
  }
  
  return(seurat_object)
}

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

options(future.globals.maxSize = 60 * 1024^3)
max_workers <- parallelly::availableCores()
cores <- min(max_workers, detectCores() - 1)
plan("multisession", workers = cores)
seurat_object <- readRDS(file.path(SEURAT_OBJECT_LOC, "pca_reduced_object.Rds"))

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
  #	genes_file <- file.path(OUT_DIR, case, "features.txt")
  #	seurat_object <- ClustGenes(
  #    seurat_object, 
  #    genes = c(as.matrix(read.table(genes_file))), 
  #    genes_list_label = case, 
  #    k = k, 
  #    res = res
  #  )
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
  # genes_file <- file.path(OUT_DIR, case, "features.txt")
  #	seurat_object <- ClustGenes(
  #    seurat_object, 
  #    genes = c(as.matrix(read.table(genes_file))), 
  #    genes_list_label = case, 
  #    k = k, 
  #    res = res
  #  )
  }

saveRDS(seurat_object, file = file.path(OUTDIR, "clustered_object.Rds"))
