# Set static variables
SCRIPT_DIR <- "sc_crohns/clustering"
SEURAT_OBJECT_LOC <- "project-area/data/crohns_scrnaseq/clustering_output"
OUTDIR <- SEURAT_OBJECT_LOC

# Source .Rprofile from sc_crohns/clustering/
source(file.path(SCRIPT_DIR, ".Rprofile"))

# Function to compute cluster composition 
ClustersComposition <- function(
    seurat_object, 
    out.prefix, 
    cl.ident.slot = "seurat_clusters", 
    sample.ident.slot = "orig.ident"
  ) {
  
  sample_ident <- seurat_object@meta.data[[sample.ident.slot]]
  samples <- unique(sample_ident)
  
  cl_ident <- seurat_object@meta.data[[cl.ident.slot]]
  clusters <- unique(cl_ident)
  
  df <- matrix(NA, nrow = length(samples)*length(clusters), ncol = 3)
  df <- as.data.frame(df)
  colnames(df) <- c("sample", "cluster", "num")
  df$num <- 0 
  
  n <- 1
  for (s in rev(samples)) {
    for (c in clusters) {
      count <- sum(sample_ident == s & cl_ident == c)
      df[n,] <- c(s, c, count)
      n <- n + 1
    }
  }
  
  df$num <- as.numeric(df$num)
  df$cluster <- factor(df$cluster, levels = sort(clusters))
  write.table(df, file = paste0(out.prefix, ".tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  g <- ggplot(df, aes(x=sample, y=num, fill=cluster)) + theme_minimal() + coord_flip() + xlab("") + theme(axis.text.y=element_text(size=15))
  g_abs <- g + geom_bar(stat="identity") + ylab("number of cells")
  pdf(paste0(out.prefix, "_SbyCabs.pdf"), width = 6, height = 3)
  print(g_abs)
  dev.off()
  
  g_norm <- g + geom_bar(stat="identity", position="fill") + ylab("fraction of cells")
  pdf(paste0(out.prefix, "_SbyCnorm.pdf"), width = 6, height = 3)
  print(g_norm)
  dev.off()
  
  g <- ggplot(df, aes(x=cluster, y=num, fill=sample)) + theme_minimal() + coord_flip() + xlab("") + theme(axis.text.y=element_text(size=15))
  g_abs <- g + geom_bar(stat="identity") + ylab("number of cells")
  pdf(paste0(out.prefix, "_CbySabs.pdf"), width = 6, height = 3)
  print(g_abs)
  dev.off()
  
  g_norm <- g + geom_bar(stat="identity", position="fill") + ylab("fraction of cells")
  pdf(paste0(out.prefix, "_CbySnorm.pdf"), width = 6, height = 3)
  print(g_norm)
  dev.off()
  
  return()
}

clustering_params <- list(
  disps = 0.5,
  nfeats = 5000,
  k = 30,
  res = c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2)
)

disps <- clustering_params$disps
nfeats <- clustering_params$nfeats
k <- clustering_params$k
res <- clustering_params$res

options(future.globals.maxSize = 60 * 1024^3)
max_workers <- parallelly::availableCores()
cores <- min(max_workers, detectCores() - 1)
plan("multisession", workers = cores, future.seed = TRUE)
seurat_object <- readRDS(file.path(SEURAT_OBJECT_LOC, "clustered_object.Rds"))

# Create list of cases and dimensionality reductions
cases <- c()
for (ymin in disps)
  cases <- c(cases, paste("mean.var.plot_disp", ymin, sep=""))
for (nfeat in nfeats)
  cases <- c(cases, paste("vst_top", nfeat, sep=""))
drs <- paste("pca", cases, sep="_")

# Create directories for each dim red
for (ll in drs) 
  dir.create(file.path(OUTDIR, ll), showWarnings = FALSE, recursive = TRUE)

# Perform cluster composition computation
for (dr in drs) {
  print(paste("dr =", dr))
  for (kk in k) {
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
      out_prefix <- file.path(OUTDIR, dr, cl_ident)
      ClustersComposition(
        seurat_object = seurat_object, 
        out.prefix = out_prefix, 
        cl.ident.slot = cl_ident, 
        sample.ident.slot = "sample.name"
      ) 
    }
  }
}