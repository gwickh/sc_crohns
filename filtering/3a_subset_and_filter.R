OBJ <- args[1]
MODE <- args[2]
CL <- args[3]
CELLS <- args[4]
OUT_OBJ <- args[5]

cls <- unlist(strsplit(CL, split=","))
dr <- gsub("_k.*", "", gsub("clusters_", "umap_", MODE))

object <- readRDS(OBJ)
dir <- dirname(CELLS)

idx <- which(!(object@meta.data[[MODE]] %in% cls))
cells <- colnames(object)[idx]

write.table(cells, file = CELLS, col.names = FALSE, row.names = FALSE, quote = FALSE)

object <- PercentageFeatureSet(object, "^MT-", col.name = "percent_mito")
pdf(file.path(dir, "umap_orig_percent_mito.pdf"))
FeaturePlot(object, reduction = dr, features = "percent_mito")
dev.off()

pdf(file.path(dir, "umap_filt_percent_mito.pdf"))
FeaturePlot(object, reduction = dr, cells = cells, features = "percent_mito")



object <- readRDS(OBJ)
cells <- read.table(CELLS)[,1]

object <- subset(object, cells = cells)
object <- ScaleData(object)

saveRDS(object, file = OUT_OBJ)