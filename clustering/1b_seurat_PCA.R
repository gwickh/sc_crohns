# Set working directory
setwd("project-area/data/crohns_scrnaseq/crohns_samples/")

# Function to run PCA
dim_reduction <- function(object, genes.list, dr = "pca", num_pcs = 50) {
  key <- paste(
    gsub('\\.|_',"",dr), 
    "_", 
    sep=""
  )
  
  object <- RunPCA(
    object, 
    features = genes.list, 
    npcs = num_pcs, 
    verbose = FALSE, 
    reduction.name = dr, 
    reduction.key = key
  )
  
  return(object)
}

# Define values for PCA parameters
params <- list(
  cores = 1,
  xmin = 0.1,
  xmax = 10,
  disps = c(0.5, 1, 1.5),
  nfeats = c(1000, 2000, 5000)
)

cores <- params$cores
xmin <- params$xmin
xmax <- params$xmax
disps <- params$disps
nfeats <- params$nfeats

# Set parallelisation
plan("multiprocess", workers = cores)
object <- readRDS("object.Rds")

# Loop over variable feature selection methods
for (ymin in disps) {
  case <- paste0("mean.var.plot_disp", ymin)
  write(case, file="")
  dr <- paste("pca", case, sep="_")
  object <- FindVariableFeatures(object, selection.method = "mean.var.plot", mean.cutoff = c(xmin, xmax), dispersion.cutoff = c(ymin, Inf))
  out_dir <- file.path(case)
  dir.create(path = out_dir, showWarnings = FALSE)
  write.table(VariableFeatures(object), file = file.path(out_dir, "features.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  object <- dim_reduction(object, genes.list = VariableFeatures(object), dr = dr)
}

for (nfeat in nfeats) {
  case <- paste0("vst_top", nfeat)
  write(case, file="")
  dr <- paste("pca", case, sep="_")
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeat)
  out_dir <- file.path(case)
  dir.create(path = out_dir, showWarnings = FALSE)
  write.table(VariableFeatures(object), file = file.path(out_dir, "features.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  object <- dim_reduction(object, genes.list = VariableFeatures(object), dr = dr)
}

# write output
saveRDS(object, file = file.path("object.Rds"))