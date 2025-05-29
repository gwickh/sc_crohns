source(file.path("sc_crohns/clustering", ".Rprofile"))

# Define values for PCA parameters
PCA_params <- list(
  xmin = 0.1,
  xmax = 10,
  disps = c(0.5, 1, 1.5),
  nfeats = c(1000, 2000, 5000)
)

xmin <- PCA_params$xmin
xmax <- PCA_params$xmax
disps <- PCA_params$disps
nfeats <- PCA_params$nfeats

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_LOC, "seurat_object.Rds"))
}
  
# Loop over variable feature selection methods
for (ymin in disps) {
  case <- paste0("mean.var.plot_disp", ymin)
  message("Processing: ", case)
  
  dr <- paste("pca", case, sep = "_")
  seurat_object <- FindVariableFeatures(
    seurat_object, 
    selection.method = "mean.var.plot", 
    mean.cutoff = c(xmin, xmax), 
    dispersion.cutoff = c(ymin, Inf)
  )
  seurat_object <- dim_reduction(
    seurat_object,
    genes.list = VariableFeatures(seurat_object),
    dr = dr,
    out_dir = file.path(OUTDIR, case)
  )
}

for (nfeat in nfeats) {
  case <- paste0("vst_top", nfeat)
  message("Processing: ", case)
  
  dr <- paste("pca", case, sep = "_")
  seurat_object <- FindVariableFeatures(
    seurat_object,
    selection.method = "vst",
    nfeatures = nfeat
  )
  seurat_object <- dim_reduction(
    seurat_object,
    genes.list = VariableFeatures(seurat_object),
    dr = dr,
    out_dir = file.path(OUTDIR, case)
  )
}

# Save Final Object
saveRDS(seurat_object, file = file.path(OUTDIR, "seurat_object.Rds"))