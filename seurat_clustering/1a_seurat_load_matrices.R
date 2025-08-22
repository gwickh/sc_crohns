# Set static variables
MIN_CELLS <- 3
MIN_FEATURES <- 200
VARS_TO_REGRESS <- c("S.Score", "G2M.Score") # Default should be NULL

source(file.path("sc_crohns/seurat_clustering", ".Rprofile"))

if (file.exists(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))) {
  print("seurat_object already created, skipping")
  q()
}

# Load matrices into Seurat obj -----------------------------------------
# Get matrix path and sample names for existing filtered matrices

matrix_paths <- list.files(
  path = MATRIX_DIR,
  pattern = "filtered_feature_bc_matrix\\.h5$",
  recursive = TRUE,
  full.names = TRUE
)

print(matrix_paths)

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
  seurat_obj <- CreateSeuratObject(
    counts = counts, 
    project = sample,
    min.cells = MIN_CELLS,
    min.features = MIN_FEATURES
  )
  
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

# Annotate mitochondrial content
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

# Filter low-quality cells
merged <- subset(
  merged, 
  subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10
)

## Normalise and perform cell cycle scoring ------------------------
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
merged <- FindVariableFeatures(merged)
merged <- ScaleData(
  merged,
  vars.to.regress = VARS_TO_REGRESS,
  features = VariableFeatures(merged)
)

# Write out
saveRDS(merged, file = file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
saveRDS(merged, file = file.path(SEURAT_OBJECT_PATH, "seurat_object_raw.Rds"))