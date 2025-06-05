source(file.path("sc_crohns/clustering", ".Rprofile"))

CLUSTERING_OUTPUT_PATH <- file.path(SEURAT_OBJECT_PATH, "clustering_stats")

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

# Functions to calculate cluster compositions and silhouette scores ------------
# compute cluster composition
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

# compute silhouette scores
SilhouetteSingle <- function(dist, ident) {
  id <- unique(ident)
  if (length(id) > 1) { 
    silh_data <- silhouette(as.numeric(ident), dist)
    silh_summ <- summary(silh_data)
    saveRDS(silh_summ, file = file.path(
      CLUSTERING_OUTPUT_PATH, 
      "silh_summary.Rds")
    )
    ggsave(
      filename = file.path(CLUSTERING_OUTPUT_PATH, "silh_summary.png"),
      plot = fviz_silhouette(silh_data), 
      width = 8, height = 6, dpi = 150
    )
  }
  return()
}

# get clustering statistics
disps <- seq(0.5, 1.5, 0.5)
n_features <- c(1000, 2000, 5000)
K <- seq(10, 50, 5)
res <- seq(0.5, 1.2, 0.1)

# Create list of cases and dimensionality reductions
cases <- c()
for (disp in disps) {
  cases <- c(cases, paste("mean.var.plot_disp", disp, sep="_"))
}
for (nfeat in nfeats) {
  cases <- c(cases, paste("vst_top", nfeat, sep="_"))
}
drs <- paste("pca", cases, sep="_")

# Create directories for each dim red
for (ll in drs) {
  dir.create(file.path(ll), showWarnings = FALSE, recursive = TRUE)
}
  
# Perform cluster composition computation
for (dr in drs) {
  print(paste("dr =", dr))
  for (k in K) {
    for (r in res) {
      cl_ident <- paste("clusters", dr, "k", k, "res", r, sep="_")
      out_prefix <- file.path(dr, cl_ident)
      ClustersComposition(
        seurat_object = seurat_object, 
        out.prefix = out_prefix, 
        cl.ident.slot = cl_ident, 
        sample.ident.slot = "orig.ident"
      )
      SilhouetteSingle(
        dist = dist, 
        ident = c(as.matrix(seurat_object@meta.data[cl_ident])), 
        out.prefix = out_prefix
      )
      
    }
  }
}

for (kk in k) {
  for (r in res) {
    n <- length(silh_data$clus.sizes)
    avg_silh <- silh_data$si.summary[4]
    meta <- seurat_object@meta.data[[cl_ident]]
    avg_umi <- round(mean(seurat_object@meta.data$nCount_RNA))
    avg_umi_cl <- unlist(lapply(1:n-1, function(x) round(mean(seurat_object@meta.data$nCount_RNA[meta == as.character(x)]))))
    idx_min_umi <- order(avg_umi_cl)[1]-1
    cl_size <- sum(meta == as.character(idx_min_umi))
    v <- c(case, dim, kk, r, n, avg_silh, avg_umi, min(avg_umi_cl), idx_min_umi, cl_size)
    df <- rbind(df, v)
  }
}

df <- data.frame(
  hvg.set = 1, 
  PCs = 1, 
  k = 1, 
  res = 1, 
  num.cl = 1, 
  avg.silh.width = 1, 
  avg.umi = 1, 
  min.umi.val = 1, 
  min.umi.idx = 1, 
  min.umi.cl.size = 1
)

df <- data.frame(
  "hvg.set", 
  "PCs", 
  "k", 
  "res", 
  "num.cl", 
  "avg.silh.width", 
  "avg.umi", 
  "min.umi.cl.val", 
  "min.umi.cl.idx", 
  "min.umi.cl.size"
)


