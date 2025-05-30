source(file.path("sc_crohns/clustering", ".Rprofile"))

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

# Function to perform dimensionality reduction by PCA ------------------
PCA_dr <- function(
  seurat_object,
  genes.list,
  dr,
  num_pcs = 50,
  out_dir
  ) {
  # Write selected features to file
  write.table(
    genes.list,
    file = file.path(paste(out_dir, "features.txt", sep = "_")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  # Run PCA
  seurat_object <- RunPCA(
    seurat_object,
    features = genes.list,
    npcs = num_pcs,
    verbose = FALSE,
    reduction.name = dr,
    reduction.key = paste(gsub('\\.|_', "", dr), "_", sep = "")
  )
  return(seurat_object)
}


# Loop over variable feature selection methods and run PCA ----------------
# Create directory to store PCA stats
PCA_OUTPUT_PATH <- file.path(SEURAT_OBJECT_PATH, "PCA_stats")
dir.create(PCA_OUTPUT_PATH, showWarnings = FALSE)

# Define values for PCA parameters
xmin = 0.1
xmax = 10
disps = c(0.5, 1, 1.5)
n_features = c(1000, 2000, 5000)
dr_li <- list()
  
# Loop over variable feature selection methods
for (disp in disps) {
  filename <- paste("mean.var.plot_disp", disp, sep = "_")
  dr_name <- paste("pca", filename, sep = "_")
  dr_li <- append(dr_li, paste("mean.var.plot_disp", disp, sep = "_"))
  
  seurat_object <- FindVariableFeatures(
    seurat_object, 
    selection.method = "mean.var.plot", 
    mean.cutoff = c(xmin, xmax), 
    dispersion.cutoff = c(disp, Inf)
  )
  # Perform PCA with PCA_dr() and write elbow plot
  seurat_object <- PCA_dr(
    seurat_object,
    genes.list = VariableFeatures(seurat_object),
    dr = dr_name,
    out_dir = file.path(PCA_OUTPUT_PATH, filename)
  )
}

for (n in n_features) {
  filename <- paste("vst_top", n, sep = "_")
  dr_name <- paste("pca", filename, sep = "_")
  dr_li <- append(dr_li, paste("vst_top", n, sep = "_"))
  
  seurat_object <- FindVariableFeatures(
    seurat_object,
    selection.method = "vst",
    nfeatures = n
  )
  # Perform PCA with PCA_dr() and write elbow plot
  seurat_object <- PCA_dr(
    seurat_object,
    genes.list = VariableFeatures(seurat_object),
    dr =  dr_name,
    out_dir = file.path(PCA_OUTPUT_PATH, filename)
  )
}

# Generate elbow and loadings visualisations ------------------------------
# Create elbow plot for each PCA performed
elbow_plots <- lapply(dr_li, function(dr) {
  ElbowPlot(
    seurat_object,
    reduction = paste("pca", dr, sep="_"),
    ndims = 50
  ) +
  ggtitle(dr) +
  theme(plot.title = element_text(size = 10))
})

combined_plot <- wrap_plots(
  elbow_plots,
  nrow = ceiling(length(dr_li) / 2)
)

ggsave(
  filename = file.path(PCA_OUTPUT_PATH, "elbow_plots.png"),
  plot = combined_plot,
  width = 18, 
  height = 4 * ceiling(length(dr_li) / 2),
  dpi = 300
)

# Create loadings plot for each PCA performed
loading_plots <- lapply(dr_li, function(dr) {
  VizDimLoadings(
    seurat_object,
    dims = 1:3,
    reduction = paste("pca", dr, sep="_")
  ) +
  ggtitle(paste0(dr, ": PCs 1-3")) +
  theme(plot.title = element_text(size = 10)) +
  ggdraw()
})

combined_loadings <- wrap_plots(
  loading_plots,
  nrow = ceiling(length(dr_li) / 2)
)
ggsave(
  filename = file.path(PCA_OUTPUT_PATH, "loading_plots_multidim.png"),
  plot = combined_loadings,
  width = 18,
  height = 4 * ceiling(length(dr_li) / 2),
  dpi = 300
)

# Save Final Object
saveRDS(seurat_object, file = file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))