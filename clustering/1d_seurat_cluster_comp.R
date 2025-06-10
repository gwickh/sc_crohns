# Load dependencies and environment
source(file.path("sc_crohns/clustering", ".Rprofile"))

CLUSTERING_OUTPUT_PATH <- file.path(SEURAT_OBJECT_PATH, "clustering_stats")
dir.create(CLUSTERING_OUTPUT_PATH, showWarnings = FALSE, recursive = TRUE)

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

# ----------------- Function: Compute Cluster Composition -----------------
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
  
  df <- as.data.frame(table(sample = sample_ident, cluster = cl_ident))
  df$num <- as.numeric(df$Freq)
  df$cluster <- factor(df$cluster, levels = sort(unique(cl_ident)))
  df$sample <- factor(df$sample, levels = rev(sort(unique(sample_ident))))
  df$Freq <- NULL
  
  out_path <- file.path(CLUSTERING_OUTPUT_PATH, paste0(out.prefix, ".tsv"))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  write.table(
    df, 
    file = out_path, 
    row.names = FALSE, 
    sep = "\t",
    quote = FALSE
  )
  
  plot_bar <- function(data, x, y, fill, norm = FALSE) {
    g <- ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
      theme_minimal() +
      coord_flip() +
      xlab("") +
      theme(axis.text.y = element_text(size = 15))
    
    if (norm) {
      g <- g + geom_bar(stat = "identity", position = "fill") + ylab("fraction of cells")
    } else {
      g <- g + geom_bar(stat = "identity") + ylab("number of cells")
    }
    
    return(g)
  }
  
  # wrapper to generate and combine all four plots
  plot_cluster_composition <- function(data, out.prefix, width = 12, height = 8) {
    p1 <- plot_bar(data, x = sample, y = num, fill = cluster, norm = FALSE) +
      ggtitle("Sample by Cluster (absolute)")
    
    p2 <- plot_bar(data, x = sample, y = num, fill = cluster, norm = TRUE) +
      ggtitle("Sample by Cluster (normalized)")
    
    p3 <- plot_bar(data, x = cluster, y = num, fill = sample, norm = FALSE) +
      ggtitle("Cluster by Sample (absolute)")
    
    p4 <- plot_bar(data, x = cluster, y = num, fill = sample, norm = TRUE) +
      ggtitle("Cluster by Sample (normalized)")
    
    combined_plot <- (p1 | p2) / (p3 | p4)
    
    ggsave(
      filename = file.path(
        CLUSTERING_OUTPUT_PATH, 
        paste0(out.prefix, "_cluster_composition.pdf")
      ),
      plot = combined_plot,
      width = width,
      height = height
    )
  }
  plot_cluster_composition(df, out.prefix)
}
# ----------------- Function: Compute Silhouette Scores -----------------
SilhouetteSingle <- function(
  dist, 
  cl_ident,
  cluster_meta = seurat_object@meta.data[[cl_ident]]
) {
  if (length(unique(cluster_meta)) > 1) { 
    silh_data <- cluster::silhouette(
      as.numeric(as.factor(cluster_meta)), 
      dist
    )
    silh_summ <- summary(silh_data)
    
    saveRDS(
      silh_summ, 
      file = file.path(CLUSTERING_OUTPUT_PATH, paste0(cl_ident, "_silh.Rds"))
    )
    ggsave(
      filename = file.path(
        CLUSTERING_OUTPUT_PATH, 
        paste0(cl_ident, "_silh.png")
      ),
      plot = factoextra::fviz_silhouette(silh_data), 
      width = 8, height = 6, dpi = 150
    )
    return(silh_summ)
  }
  return(NULL)
}

# ----------------- Function: Summarize Clustering -----------------
SummarizeClusteringStats <- function(
    silh_data,
    seurat_object,
    cl_ident,
    case,
    dim,
    k,
    r
) {
  n_clusters <- length(silh_data$clus.sizes)
  avg_silh <- as.numeric(silh_data$avg.width)
  cluster_meta <- seurat_object@meta.data[[cl_ident]]
  avg_umi <- round(mean(seurat_object@meta.data$nCount_RNA))
  
  avg_umi_per_cluster <- sapply(
    unique(cluster_meta), 
    function(cl) round(mean(seurat_object@meta.data$nCount_RNA[cluster_meta == cl]))
  )
  names(avg_umi_per_cluster) <- as.character(unique(cluster_meta)) 
  
  idx_min_umi <- which.min(avg_umi_per_cluster)
  min_avg_umi <- avg_umi_per_cluster[idx_min_umi]
  min_cluster_size <- sum(cluster_meta == names(avg_umi_per_cluster)[idx_min_umi])
  min_umi_cluster <- names(avg_umi_per_cluster)[idx_min_umi]
  
  return(data.frame(
    case = case,
    dim = dim,
    k = k,
    r = r,
    n_clusters = n_clusters,
    avg_silhouette = avg_silh,
    avg_umi = avg_umi,
    min_avg_umi = min_avg_umi,
    min_umi_cluster = min_umi_cluster,
    min_cluster_size = min_cluster_size
  ))
}

# ----------------- Parameters -----------------
disp_cases <- paste0("mean.var.plot_disp_", disps)
vst_cases  <- paste0("vst_top_", n_features)
cases <- c(disp_cases, vst_cases)

# ----------------- Run Analysis -----------------
df <- data.frame()

for (case in cases) {
  dr <- paste0("pca_", case)
  dir.create(dr, showWarnings = FALSE, recursive = TRUE)
  dim_file <- file.path(CLUSTERING_OUTPUT_PATH, paste(dr, "num_PCs.txt", sep = "_"))
  dim <- c(as.matrix(read.table(dim_file)))
  data <- as.matrix(seurat_object@reductions[[dr]]@cell.embeddings[,1:dim])
  dist <- parDist(data, threads = cores)
  
  for (k in neighbors) {
    for (r in res) {
      cl_ident <- paste("clusters", dr, "k", k, "res", r, sep = "_")
      out_prefix <- file.path(dr, cl_ident)
      
      if (!(cl_ident %in% colnames(seurat_object@meta.data))) next
      
      ClustersComposition(
        seurat_object = seurat_object, 
        out.prefix = out_prefix, 
        cl.ident.slot = cl_ident, 
        sample.ident.slot = "orig.ident"
      )
      
      silh_data <- SilhouetteSingle(
        dist = dist, 
        cl_ident = cl_ident
      )

      if (!is.null(silh_data)) {
        summary <- SummarizeClusteringStats(
          silh_data = silh_data,
          seurat_object = seurat_object,
          cl_ident = cl_ident,
          case = case,
          dim = dim,
          k = k,
          r = r
        )
        df <- rbind(df, summary)
      }
    }
  }
}

# ----------------- Save Summary -----------------
write.table(
  df,
  file = file.path(CLUSTERING_OUTPUT_PATH, "clustering_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)