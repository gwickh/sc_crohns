MODE <- c(
  "clusters_pca_mean.var.plot_disp_0.5_k_10_res_1.2" 
)
CL <- c("44")
OUT <- cells.tsv
OUT_OBJ <- "obj"
backup_object <- seurat_object

cls <- unlist(strsplit(CL, split=","))
dr <- gsub("_k.*", "", gsub("clusters_", "umap_", MODE))

object <- readRDS(seurat_object)
dir <- dirname(OUT)

idx <- which(!(seurat_object@meta.data[[MODE]] %in% cls))
cells <- colnames(seurat_object)[idx]

write.table(cells, file = CELLS, col.names = FALSE, row.names = FALSE, quote = FALSE)
seurat_object <- subset(seurat_object, cells = cells)
seurat_object <- ScaleData(seurat_object)


seurat_object <- PercentageFeatureSet(seurat_object, "^MT-", col.name = "percent_mito")
plot <- FeaturePlot(seurat_object, reduction = dr, features = "percent_mito")
pdf("percent_mito_plot.pdf")
print(plot)
dev.off()

plot <- FeaturePlot(seurat_object,  reduction = dr, features = "nCount_RNA")    # Total UMI count
pdf("nCount_RNA.pdf")
print(plot)
dev.off()

plot <- FeaturePlot(seurat_object,  reduction = dr, features = "nFeature_RNA")
pdf("nFeature_RNA_plot.pdf")
print(plot)
dev.off()

seurat_object <- subset(seurat_object, subset = percent_mito < 10)

saveRDS(seurat_object, file = OUT_OBJ)