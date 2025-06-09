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

# initialise empty lists
dr_li <- list()
vfeature_objects <- list()
  
# Loop over feature selection methods
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
  
  # save variable features  
  vfeature_objects[[paste(filename, "vfeatures", sep = "_")]] <- seurat_object

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
  # save variable features  
  vfeature_objects[[paste(filename, "vfeatures", sep = "_")]] <- seurat_object
  
  # Perform PCA with PCA_dr() and write elbow plot
  seurat_object <- PCA_dr(
    seurat_object,
    genes.list = VariableFeatures(seurat_object),
    dr =  dr_name,
    out_dir = file.path(PCA_OUTPUT_PATH, filename)
  )
}

# Generate elbow, PC loadings and variable features visualisations ------------
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
    reduction = paste("pca", dr, sep="_"),
    ncol = 3
  ) +
  ggtitle(paste0(dr, ": PCs 1-3")) +
  theme(plot.title = element_text(size = 10))
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

# Plot highest variable genes
vfp_list <- list()
for (name in names(vfeature_objects)) {
  hvf <- HVFInfo(vfeature_objects[[name]]) 
  if ("variance.standardized" %in% colnames(hvf)) {
    top_genes <- hvf %>%
      arrange(desc(variance.standardized)) %>%
        head(10) %>%
          rownames()
  } else if ("mvp.dispersion.scaled" %in% colnames(hvf)) {
    top_genes <- hvf %>%
      arrange(desc(mvp.dispersion.scaled)) %>%
        head(10) %>%
          rownames()
  } else {
    warning(
      paste0("No valid HVF metric found for ", name, 
        "\ncolumn names = ", paste(colnames(hvf), collapse = ", ")))
    next
  }
  
  vfplot <- VariableFeaturePlot(vfeature_objects[[name]])
  vfplot <- LabelPoints(
    plot = vfplot,
    points = top_genes, 
    repel = TRUE,
    xnudge = 0,
    ynudge = 0
  ) +
    ggtitle(paste("Variable Features:", name)) +
    theme(plot.title = element_text(size = 10))
  vfp_list[[name]] <- vfplot
}

combined_vfeatures <- wrap_plots(
  vfp_list,
  ncol = ceiling(length(dr_li) / 2)
)
  
ggsave(
  filename = file.path(PCA_OUTPUT_PATH, "v_features.png"),
  plot = combined_vfeatures,
  width = 18, height = 5, dpi = 300
)

# Save Final Object
saveRDS(seurat_object, file = file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))