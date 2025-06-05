source(file.path("sc_crohns/clustering", ".Rprofile"))


# Function to calculate cluster compositions ------------------------------
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

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_LOC, "seurat_object.Rds"))
}

# get clustering statistics
disps <- seq(0.5, 1.5, 0.5)
n_features <- c(1000, 2000, 5000)
k <- seq(10, 50, 5)
res <- seq(0.5, 1.2, 0.1)

# Create list of cases and dimensionality reductions
cases <- c()
for (ymin in disps)
  cases <- c(cases, paste("mean.var.plot_disp", ymin, sep=""))
for (nfeat in nfeats)
  cases <- c(cases, paste("vst_top", nfeat, sep=""))
drs <- paste("pca", cases, sep="_")

# Create directories for each dim red
for (ll in drs) 
  dir.create(file.path(ll), showWarnings = FALSE, recursive = TRUE)

# Perform cluster composition computation
for (dr in drs) {
  print(paste("dr =", dr))
  for (kk in k) {
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
      out_prefix <- file.path(dr, cl_ident)
      ClustersComposition(
        seurat_object = seurat_object, 
        out.prefix = out_prefix, 
        cl.ident.slot = cl_ident, 
        sample.ident.slot = "orig.ident"
      ) 
    }
  }
}