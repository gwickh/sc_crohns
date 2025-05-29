DoubletFind <- function(object, dr = "pca", pcs = 50, expDoublPerc = 2, out_dir) {
  
  names(object@reductions)[names(object@reductions) == dr] <- "pca"
  
  pdf(file.path(out_dir, "doubletFinder.pdf"))
  
  sweep.res.list <- paramSweep_v3(object, PCs = 1:pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal_pK <- bcmvn$pK[which(bcmvn$BCmetric > max(bcmvn$BCmetric) - EPSILON)]
  optimal_pK <- as.numeric(as.character(optimal_pK))
  
  expDoubl <- round(ncol(object)*expDoublPerc/100)
  object <- doubletFinder_v3(object, PCs = 1:pcs, pN = pN, pK = optimal_pK, nExp = expDoubl, reuse.pANN = FALSE, sct = FALSE)
  
  dev.off()
  
  id <- paste(pN, optimal_pK, expDoubl, sep = "_")
  score_slot <- paste0("pANN_", id)
  label_slot <- paste0("DF.classifications_", id)
  
  colnames(object@meta.data)[colnames(object@meta.data) == score_slot] <- paste0("doublet_score_", dr)
  colnames(object@meta.data)[colnames(object@meta.data) == label_slot] <- paste0("doublet_label_", dr)
  names(object@reductions)[names(object@reductions) == "pca"] <- dr
  
  return(object)
}

dir.create(OUT_DIR, showWarnings = FALSE)

object <- readRDS(OBJECT)

feature_method="vst"
for (nfeat in nfeats) {
  
  case <- paste("vst_top", nfeat, sep = "")
  dr <- paste("pca", case, sep = "_")
  in_dir <- file.path(IN_DIR, case)
  out_dir <- file.path(OUT_DIR, dr)
  dir.create(out_dir, showWarnings = FALSE)
  
  ExtractPCs(object, dr = dr, file = file.path(in_dir, "num_PCs.txt"))
  nPCs <- GetNPCs(in_dir)
  object <- DoubletFind(object, dr = dr, pcs = nPCs, expDoublPerc = DOUBLET_PERC, out_dir = out_dir)
  
  tsne_dr <- paste0("tsne_", dr)
  umap_dr <- paste0("umap_", dr)
  if (!(tsne_dr %in% names(object@reductions)))
    object <- RunTSNE(object, reduction = dr, dims = 1:nPCs, reduction.name = tsne_dr, reduction.key = tsne_dr)
  if (!(umap_dr %in% names(object@reductions)))
    object <- RunUMAP(object, reduction = dr, dims = 1:nPCs, reduction.name = umap_dr, reduction.key = umap_dr)
  DrPlot(object, dr = dr, dr_plot = tsne_dr, out_dir = out_dir, dr_type = "TSNE")
  DrPlot(object, dr = dr, dr_plot = umap_dr, out_dir = out_dir, dr_type = "UMAP")
}

saveRDS(object, file=OBJECT)

object <- readRDS(CL_OBJ)
cl_all <- object@meta.data[[SLOT]]
names(cl_all) <- colnames(object)

object <- readRDS(DOUBL_OBJ)
cells <- colnames(object)

dr <- gsub("_k.*", "", gsub("clusters_pca_", "", SLOT))
dr <- paste0("pca_", dr)
doublet_label_slot <- paste0("doublet_label_", dr)
doublet_score_slot <- paste0("doublet_score_", dr)

df <- object@meta.data[,match(c(doublet_label_slot, doublet_score_slot),colnames(object@meta.data))]
rownames(df) <- colnames(object)
colnames(df) <- c("doublet_label", "doublet_score")

OUT_TSV <- paste0(OUT_PREFIX, "_res.tsv")
write.table(df, file = OUT_TSV, sep = "\t", quote = FALSE)

cl <- cl_all[match(cells, names(cl_all))]
cl <- factor(cl, levels = sort(unique(cl)))

tab <- table(data.frame(clusters = cl, doublets = df[,1]))
doublet.frac <- tab[,"Doublet"] / rowSums(tab)
df <- data.frame(clusters = rownames(tab), doublets = tab[,"Doublet"], doublet.frac = doublet.frac)

OUT_TSV <- paste0(OUT_PREFIX, "_stats.tsv")
write.table(df, file = OUT_TSV, sep = "\t", row.names = FALSE, quote = FALSE)