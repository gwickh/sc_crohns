# Load dependencies and environment
source(file.path("sc_crohns/clustering", ".Rprofile"))

FILTERING_OUTPUT_PATH <- file.path(SEURAT_OBJECT_PATH, "filtering_stats")
dir.create(FILTERING_OUTPUT_PATH, showWarnings = FALSE, recursive = TRUE)

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

# Initialize summary storage
all_summaries <- list()

# ------------------------ FUNCTION: Detect doublets --------------------------
runDoubletDetection <- function(
    seurat_obj,
    dr = "pca",
    pcs = 50,
    expDoublPerc = 2,
    epsilon = 0.0001,
    pN = 0.25,
    clusters = "seurat_clusters",
    out_dir,
    layer_name
) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    samp_split <- SplitObject(seurat_obj, split.by = "sample_id")

    samp_split <- lapply(samp_split, run_doubletfinder_custom, multiplet_rate = 0.01)

    for (sample in samp_split) {
        # Basic safety check
        if (!"seurat_clusters" %in% colnames(sample@meta.data)) next

        annotations <- sample@meta.data[[clusters]]
        homotypic.prop <- modelHomotypic(annotations)

        min_pc <- min(pcs, length(sample@reductions[[dr]]@cell.embeddings[1,]))
        sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)
        sweep_stats <- summarizeSweep(sweep_list)
        bcmvn <- find.pK(sweep_stats)
        optimal.pk <- bcmvn %>%
          filter(BCmetric == max(BCmetric, na.rm = TRUE)) %>%
          pull(pK)

        nExp.poi <- round(0.01 * nrow(sample@meta.data))
        nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

        sample <- DoubletFinder::doubletFinder_v3(
          sample,
          PCs = 1:pcs,
          pN = pN,
          pK = as.numeric(as.character(optimal.pk[1])),
          nExp = nExp.poi.adj,
          reuse.pANN = FALSE,
          sct = FALSE
        )
    }

    # Optionally reassemble and return
    combined <- merge(samp_split[[1]], samp_split[-1])
    
    # Save summary stats
    summary_file <- file.path(out_dir, "doublet_cluster_stats.tsv")
    df <- table(combined$DF.classifications) # Adjust this to reflect your naming convention
    write.table(df, summary_file, sep = "\t", quote = FALSE, col.names = NA)
    
    # Store for global summary
    df <- as.data.frame(df)
    df$layer <- layer_name
    df$dr <- dr
    all_summaries[[layer_name]] <<- df
    
    return(combined)
}

# Set the number of PCs
pcs_to_use <- 50

# User-defined values — make sure these are declared
n_features <- c(1000, 2000, 3000)
disps <- c(0.5, 1.0, 1.5)

dr_li <- list()

# Loop over VST feature selections
for (n in n_features) {
  filename <- paste0("vst_top_", n)
  dr_name <- paste0("pca_", filename)
  dr_li <- append(dr_li, filename)
  
  seurat_object <- runDoubletDetection(
    seurat_obj = seurat_object,
    dr = dr_name,
    pcs = pcs_to_use,
    expDoublPerc = 7.8,
    clusters = "seurat_clusters",
    out_dir = file.path(FILTERING_OUTPUT_PATH, filename),
    layer_name = filename
  )
}

# Loop over dispersion values
for (disp in disps) {
  filename <- paste0("mean.var.plot_disp_", disp)
  dr_name <- paste0("pca_", filename)
  dr_li <- append(dr_li, filename)
  
  seurat_object <- runDoubletDetection(
    seurat_obj = seurat_object,
    dr = dr_name,
    pcs = pcs_to_use,
    expDoublPerc = 7.8,
    clusters = "seurat_clusters",
    out_dir = file.path(FILTERING_OUTPUT_PATH, filename),
    layer_name = filename
  )
}

# Combine and save all summaries
combined_summary <- do.call(rbind, all_summaries)
combined_file <- file.path(FILTERING_OUTPUT_PATH, "combined_doublet_cluster_stats.tsv")
write.table(combined_summary, combined_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Combined summary written to: ", combined_file)

# Save Final Object
saveRDS(seurat_object, file = file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))