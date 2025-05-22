# Set static variables
SCRIPT_DIR <- "sc_crohns/clustering"
SEURAT_OBJECT_LOC <- "project-area/data/crohns_scrnaseq/clustering_output"
OUTDIR <- SEURAT_OBJECT_LOC

# Source .Rprofile from sc_crohns/clustering/
source(file.path(SCRIPT_DIR, ".Rprofile"))

# Function to run PCA
dim_reduction <- function(seurat_object, genes.list, dr = "pca", num_pcs = 50) {
  key <- paste(
    gsub('\\.|_',"",dr), 
    "_", 
    sep=""
  )
  
  seurat_object <- RunPCA(
    seurat_object, 
    features = genes.list, 
    npcs = num_pcs, 
    verbose = FALSE, 
    reduction.name = dr, 
    reduction.key = key
  )
  
  return(seurat_object)
}

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

# Set parallelisation
max_workers <- parallelly::availableCores()
cores <- min(max_workers, detectCores() - 1)
plan("multisession", workers = cores)
seurat_object <- readRDS(file.path(SEURAT_OBJECT_LOC, "merged_object.Rds"))

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
  
  out_dir <- case
  dir.create(out_dir, showWarnings = FALSE)
  
  write.table(
    VariableFeatures(seurat_object),
    file = file.path(out_dir, "features.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  seurat_object <- dim_reduction(seurat_object, genes.list = VariableFeatures(seurat_object), dr = dr)
}

# Variable Feature Selection + PCA: vst
for (nfeat in nfeats) {
  case <- paste0("vst_top", nfeat)
  message("Processing: ", case)
  
  dr <- paste("pca", case, sep = "_")
  seurat_object <- FindVariableFeatures(
    seurat_object,
    selection.method = "vst",
    nfeatures = nfeat
  )
  
  out_dir <- case
  dir.create(out_dir, showWarnings = FALSE)
  
  write.table(
    VariableFeatures(seurat_object),
    file = file.path(out_dir, "features.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  seurat_object <- dim_reduction(seurat_object, genes.list = VariableFeatures(seurat_object), dr = dr)
}

# Save Final Object
saveRDS(seurat_object, file = file.path(OUTDIR, "pca_reduced_object.Rds"))