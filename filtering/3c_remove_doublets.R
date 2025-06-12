IN_OBJ <- seurat_object
CELLS <- cells
OUT_OBJ <- seurat_object
DR <- dr
CL <- cluster

if (length(args) > 5) { # optionally apply by-cluster filtering
  cls <- unlist(strsplit(CL, split=","))
  object <- readRDS(IN_OBJ)
  idx <- which(!(object@meta.data[[CL]] %in% cls))
  cells <- colnames(object)[idx]
} 

object <- readRDS(OBJ2)
doublet_score <- paste0("doublet_score_", DR)
doublet_label <- paste0("doublet_label_", DR)

if (length(args) < 6) {
  cells <- colnames(object)
}

doublet_scores <- object@meta.data[[doublet_score]]
total_doublets <- sum(object@meta.data[[doublet_label]] == "Doublet")
doublets_filt_by_cluster <- sum(!(colnames(object) %in% cells))
doublets_to_filter <- total_doublets - doublets_filt_by_cluster

idx_filtered <- which(colnames(object) %in% cells)
threshold <- sort(doublet_scores[idx_filtered], decreasing = TRUE)[doublets_to_filter]
idx_filtered_2 <- which(doublet_scores[idx_filtered] < threshold)

cells <- colnames(object)[idx_filtered[idx_filtered_2]]

write.table(cells, file = CELLS, col.names = FALSE, row.names = FALSE, quote = FALSE)

object <- readRDS(IN_OBJ)
cells <- read.table(CELLS)[,1]

object <- subset(object, cells = cells)
object <- ScaleData(object)

saveRDS(object, file = OBJ)

