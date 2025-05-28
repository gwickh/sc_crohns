# Perform dimensionality reduction by PCA
dim_reduction <- function(
    seurat_object,
    genes.list,
    dr = "pca",
    num_pcs = 50,
    out_dir = OUTDIR
) {
  # Create output directory if provided
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE)
    
    # Save variable features
    write.table(
      genes.list,
      file = file.path(out_dir, "features.txt"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  # Define reduction key
  key <- paste(gsub('\\.|_', "", dr), "_", sep = "")
  
  # Run PCA
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

# Function to identify number of informative PCs based on standard deviation 
# cutoff and cluster with Louvain at varying on resolution and kNN
ClustDr <- function(
    seurat_object, 
    dr = "pca", 
    pcs = 50, 
    pval = 1e-05, 
    non.rand.sd.frac = 0.5, 
    k, 
    res, 
    file
) {
  # Compute minimum std a PC must have to be considered “non-random” assuming
  # last 10 PCs approximate random noise to exclude SVD artifacts
  min_stdev <- (1+non.rand.sd.frac)*mean(seurat_object@reductions[[dr]]@stdev[(pcs-10):pcs])
  
  # Loops through PCs to determine nPCs, the number of informative PCs that 
  # explains > 0.05 variance with respect to the next
  for (nPCs in 0:(length(seurat_object@reductions[[dr]]@stdev)-1)) {
    if (seurat_object@reductions[[dr]]@stdev[nPCs+1] < min_stdev) 
      break
  }
  write(nPCs, file = file)
  
  # For each k build a NN graph using the top nPCs components
  for (kk in k) {
    seurat_object <- FindNeighbors(
      seurat_object, 
      k.param = kk, 
      reduction = dr, 
      dims = 1:nPCs
    )
    # For each resolution, run Louvain clustering
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
      seurat_object <- FindClusters(seurat_object, resolution = r)
      names(seurat_object@meta.data)[names(seurat_object@meta.data) == "seurat_clusters"] <- cl_ident
      seurat_object@active.ident <- as.factor(seurat_object$orig.ident)
    }
  }
  
  return(seurat_object)
}


ClustersComposition <- function(
    seurat_object, 
    out.prefix, 
    cl.ident.slot,
    sample.ident.slot = "orig.ident",
    width = 6, 
    height = 3
) {
  if (!all(c(cl.ident.slot, sample.ident.slot) %in% colnames(seurat_object@meta.data))) {
    stop("One or both of the specified metadata slots do not exist.")
  }
  
  sample_ident <- seurat_object@meta.data[[sample.ident.slot]]
  cl_ident <- seurat_object@meta.data[[cl.ident.slot]]
  
  # Count cells per (sample, cluster)
  df <- as.data.frame(table(sample = sample_ident, cluster = cl_ident))
  
  # Ensure numeric values and correct factor levels
  df$num <- as.numeric(df$Freq)
  df$cluster <- factor(df$cluster, levels = sort(unique(cl_ident)))
  df$sample <- factor(df$sample, levels = rev(sort(unique(sample_ident))))
  df$Freq <- NULL  # remove redundant column
  
  # Save raw counts
  write.table(df, file = paste0(out.prefix, ".tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Plot function
  plot_bar <- function(data, x, y, fill, filename, norm = FALSE) {
    g <- ggplot(data, aes_string(x = x, y = y, fill = fill)) + 
      theme_minimal() + 
      coord_flip() + 
      xlab("") + 
      theme(axis.text.y = element_text(size = 15))
    
    if (norm) {
      g <- g + geom_bar(stat = "identity", position = "fill") + ylab("fraction of cells")
    } else {
      g <- g + geom_bar(stat = "identity") + ylab("number of cells")
    }
    
    pdf(filename, width = width, height = height)
    print(g)
    dev.off()
  }
  
  # Generate plots
  plot_bar(df, "sample", "num", "cluster", paste0(out.prefix, "_SbyC_abs.pdf"), norm = FALSE)
  plot_bar(df, "sample", "num", "cluster", paste0(out.prefix, "_SbyC_norm.pdf"), norm = TRUE)
  plot_bar(df, "cluster", "num", "sample", paste0(out.prefix, "_CbyS_abs.pdf"), norm = FALSE)
  plot_bar(df, "cluster", "num", "sample", paste0(out.prefix, "_CbyS_norm.pdf"), norm = TRUE)
  
  return()
}


DrPlot <- function(
    seurat_object, 
    dr = "pca", 
    dr_plot = "pca", 
    pcs = 50, 
    pval = 1e-05, 
    non.rand.sd.frac = 0.5, 
    k, 
    res, 
    OUTDIR, 
    dr_type = "PCA"
) {
  ggsave_plot <- function(plot_obj, filename, width = 7.4, height = 6) {
    ggsave(filename, plot = plot_obj, width = width, height = height)
  }
  dimplot_configs <- list(
    list(group.by = "orig.ident", file_suffix = "samples"),
    list(group.by = "Phase", file_suffix = "ccScore")
  )
  for (cfg in dimplot_configs) {
    if (cfg$group.by %in% colnames(seurat_object@meta.data)) {
      p <- DimPlot(seurat_object, group.by = cfg$group.by, reduction = dr_plot)
      ggsave_plot(p, file.path(OUTDIR, paste0(dr_type, "_plot_", cfg$file_suffix, ".pdf")))
    }
    if ("nCount_RNA" %in% colnames(seurat_object@meta.data)) {
      fp <- FeaturePlot(seurat_object, features = "nCount_RNA", reduction = dr_plot, cols = c("blue", "red"))
      ggsave_plot(fp, file.path(OUTDIR, paste0(dr_type, "_plot_UMIcount.pdf")))
      
      dp <- ggplot(seurat_object@meta.data, aes(x = nCount_RNA)) +
        theme_classic() + geom_density() + ylab("density")
      ggsave_plot(dp, file.path(OUTDIR, "density_plot_UMIcount.pdf"), width = 3.5, height = 3.5)
    }
  }
  for (kk in k) {
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep = "")
      if (cl_ident %in% colnames(seurat_object@meta.data)) {
        p <- DimPlot(seurat_object, group.by = cl_ident, reduction = dr_plot)
        ggsave_plot(p, file.path(OUTDIR, paste0(dr_type, "_plot_", cl_ident, ".pdf")))
      }
    }
  }
}
run_UMAP <- function() {
  case <- paste("mean.var.plot_disp", ymin, sep = "")
  dr <- paste("pca", case, sep = "_")
  tsne_dr <- paste0("tsne_", dr)
  umap_dr <- paste0("umap_", dr)
  in_dir <- file.path(SEURAT_OBJECT_LOC, case)
  WRITEOUT <- file.path(OUTDIR, dr)
  dir.create(WRITEOUT, showWarnings = FALSE)
  nPCs <- GetNPCs(in_dir)
  if (!(tsne_dr %in% names(seurat_object@reductions)))
    seurat_object <- RunTSNE(seurat_object, reduction = dr, dims = 1:nPCs, reduction.name = tsne_dr, reduction.key = paste0("tSNE_", dr, "_"))
  if (!(umap_dr %in% names(seurat_object@reductions)))
    seurat_object <- RunUMAP(seurat_object, reduction = dr, dims = 1:nPCs, reduction.name = umap_dr, reduction.key = paste0("UMAP_", dr, "_"))
  DrPlot(seurat_object, dr = dr, dr_plot = dr, k = k, res = res, OUTDIR = WRITEOUT, dr_type = "PCA")	
  DrPlot(seurat_object, dr = dr, dr_plot = tsne_dr, k = k, res = res, OUTDIR = WRITEOUT, dr_type = "TSNE")
  DrPlot(seurat_object, dr = dr, dr_plot = umap_dr, k = k, res = res, OUTDIR = WRITEOUT, dr_type = "UMAP")
}