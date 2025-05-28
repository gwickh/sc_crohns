# compute cluster composition 
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

# Create list of cases and dimensionality reductions
cases <- c()
for (ymin in disps)
  cases <- c(cases, paste("mean.var.plot_disp", ymin, sep=""))
for (nfeat in nfeats)
  cases <- c(cases, paste("vst_top", nfeat, sep=""))
drs <- paste("pca", cases, sep="_")

# Create directories for each dim red
for (ll in drs) 
  dir.create(file.path(ll), showWarnings = FALSE, recursive = TRUE)

# Perform cluster composition computation
for (dr in drs) {
  print(paste("dr =", dr))
  for (kk in k) {
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
      out_prefix <- file.path(dr, cl_ident)
      ClustersComposition(
        seurat_object = seurat_object, 
        out.prefix = out_prefix, 
        cl.ident.slot = cl_ident, 
        sample.ident.slot = "orig.ident"
      ) 
    }
  }
}