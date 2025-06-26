# Load dependencies and environment
source(file.path("sc_crohns/clustering", ".Rprofile"))

FILTERING_OUTPUT_PATH <- file.path(SEURAT_OBJECT_PATH, "filtering_stats")
dir.create(FILTERING_OUTPUT_PATH, showWarnings = FALSE, recursive = TRUE)

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

# ------------------------ FUNCTION: Detect doublets --------------------------
runDoubletDetection <- function(
    seurat_obj,
    dr = "pca",
    pcs = 50,
    expDoublPerc = 2,
    epsilon = 0.0001,
    pN = 0.25,
    cluster_col,
    out_dir
) {
  all_summaries <- list()
  
  count_layers <- grep("^counts", names(seurat_obj[["RNA"]]@layers), value = TRUE)
  message("Count layers found: ", paste(count_layers, collapse = ", "))

  for (layer_name in count_layers) {
    
    counts_matrix <- seurat_obj[["RNA"]]@layers[[layer_name]]
    colnames(counts_matrix) <- colnames(seurat_obj)
    
    new_assay <- CreateAssayObject(counts = counts_matrix)
    
    # Add this assay to the Seurat object
    seurat_obj[[assay_name]] <- new_assay
    
    # Set default assay to new assay
    DefaultAssay(seurat_obj) <- assay_name
    
    # Output directory per layer
    out_dir <- file.path(out_dir, layer_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Run DoubletFind, capture returned object and get summary_df
    seurat_obj <- DoubletFind(
      object = seurat_obj,
      dr = dr,
      pcs = pcs,
      expDoublPerc = expDoublPerc,
      epsilon = epsilon,
      pN = pN,
      out_dir = out_dir,
      cluster = cluster_col
    )
    
    # Extract the summary file you wrote inside DoubletFind
    summary_file <- file.path(out_dir, "doublet_cluster_stats.tsv")
    if (file.exists(summary_file)) {
      summary_df <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      summary_df$layer <- layer_name
      summary_df$dr <- dr
      
      all_summaries[[layer_name]] <- summary_df
    } else {
      warning("Summary file not found for layer: ", layer_name)
    }
  }
  
  # Combine all summaries into one data.frame
  combined_summary <- do.call(rbind, all_summaries)
  
  # Write combined summary to disk
  combined_file <- file.path(out_dir, paste0("combined_doublet_cluster_stats_", dr, ".tsv"))
  write.table(combined_summary, combined_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Combined summary written to: ", combined_file)
  
  return(seurat_obj)
}


dr_li <- list()

for (n in n_features) {
  filename <- paste("vst_top", n, sep = "_")
  dr_name <- paste("pca", filename, sep = "_")
  dr_li <- append(dr_li, paste("vst_top", n, sep = "_"))
  
  object <- runDoubletDetection(
    seurat_obj = seurat_object, 
    dr = dr_name, 
    pcs = dim, 
    expDoublPerc = 7.8,
    cluster = "seurat_clusters",
    out_dir = file.path(FILTERING_OUTPUT_PATH, filename)
  )
}

for (disp in disps) {
  filename <- paste("mean.var.plot_disp", disp, sep = "_")
  dr_name <- paste("pca", filename, sep = "_")
  dr_li <- append(dr_li, paste("mean.var.plot_disp", disp, sep = "_"))
  
  object <- runDoubletDetection(
    seurat_obj = seurat_object, 
    dr = dr_name, 
    pcs = dim, 
    expDoublPerc = 7.8,
    cluster = "seurat_clusters",
    out_dir = file.path(FILTERING_OUTPUT_PATH, filename)
  )
}


# Save Final Object
saveRDS(seurat_object, file = file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))