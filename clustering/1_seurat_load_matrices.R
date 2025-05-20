# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load packages
pkgs <- c(
  "Seurat",
  "tidyverse",
  "Matrix",
  "hdf5r"
)

local({
  for (pkg in pkgs) {
    if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
})

# Set working directory
setwd("project-area/data/crohns_scrnaseq/crohns_samples/")

## Load matrices
# Get matrix path and sample names for existing filtered matrices
matrix_paths <- list.files(
  pattern = "filtered_feature_bc_matrix\\.h5$",
  recursive = TRUE,
  full.names = TRUE
)

sample_names <- basename(dirname(dirname(matrix_paths)))

# Load matrices into Seurat object
seurat_list <- list()

# Loop through matrix path and sample name to load count matrix from .h5 file 
for (i in seq_along(matrix_paths)) {
  path <- matrix_paths[i]
  sample <- sample_names[i]
   
  counts <- Read10X_h5(path)
  
  # Extract gene expression data
  if (is.list(counts)) {
    counts <- counts[["Gene Expression"]]
  }
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample)
  
  # Add sample metadata
  seurat_obj$sample_id <- sample
  
  # Store in list
  seurat_list[[sample]] <- seurat_obj
}

# Merge all into one Seurat object
merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_names,
  project = "crohns_samples_merged"
)

## Run downstream Seurat functions
# Normalize data
merged <- NormalizeData(merged)

# Cell cycle scoring
if (!exists("cc.genes")) {
  cc.genes <- Seurat::cc.genes.updated.2019
}
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

merged <- CellCycleScoring(
  merged,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# Regress variables
VARS_TO_REGRESS <- c("S.Score", "G2M.Score")

merged <- ScaleData(
  merged,
  vars.to.regress = VARS_TO_REGRESS
)

# Write out
saveRDS(merged, file = merged_object.Rds)