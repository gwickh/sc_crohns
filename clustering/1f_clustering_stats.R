source(file.path("sc_crohns/clustering", ".Rprofile"))

params <- as.matrix(read.table(PARAMS, sep="="))
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))
k <- as.numeric(unlist(strsplit(params[params[,1]=="k", 2], split=",")))
res <- as.numeric(unlist(strsplit(params[params[,1]=="res", 2], split=",")))

object <- readRDS(OBJ)

df <- data.frame(hvg.set = NULL, PCs = NULL, k = NULL, res = NULL, num.cl = NULL, avg.silh.width = NULL, avg.umi = NULL, min.umi.val = NULL, min.umi.idx = NULL, min.umi.cl.size = NULL)

cases <- c()
for (ymin in disps)
  cases <- c(cases, paste("mean.var.plot_disp", ymin, sep=""))
for (nfeat in nfeats)
  cases <- c(cases, paste("vst_top", nfeat, sep=""))

for (case in cases) {
  dr <- paste("pca", case, sep="_")
  dim_file <- file.path(PC_DIR, case, "num_PCs.txt")
  dim <- read.table(dim_file)[,1]
  for (kk in k) {
    for (r in res) {
      cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
      silh_prefix <- file.path(SILH_DIR, dr, cl_ident)
      silh_path <- paste0(silh_prefix, ".Rds")
      silh <- readRDS(silh_path)
      n <- length(silh$clus.sizes)
      avg_silh <- silh$si.summary[4]
      meta <- object@meta.data[[cl_ident]]
      avg_umi <- round(mean(object@meta.data$nCount_RNA))
      avg_umi_cl <- unlist(lapply(1:n-1, function(x) round(mean(object@meta.data$nCount_RNA[meta == as.character(x)]))))
      idx_min_umi <- order(avg_umi_cl)[1]-1
      cl_size <- sum(meta == as.character(idx_min_umi))
      v <- c(case, dim, kk, r, n, avg_silh, avg_umi, min(avg_umi_cl), idx_min_umi, cl_size)
      df <- rbind(df, v)
    }
  }
}

colnames(df) <- data.frame("hvg.set", "PCs", "k", "res", "num.cl", "avg.silh.width", "avg.umi", "min.umi.cl.val", "min.umi.cl.idx", "min.umi.cl.size")

write.table(df, file = OUT_TSV, sep = "\t", row.names = FALSE, quote = FALSE)
