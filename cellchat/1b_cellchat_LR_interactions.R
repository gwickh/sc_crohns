library(CellChat)
library(Seurat)
library(Matrix)
library(patchwork)
library(igraph)
library(future)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)

# Vars
path <- path.expand(
    "~/project-area/data/crohns_scrnaseq/10c_14n_analysis/cellchat_output"
)

# load files
cellchat.crohns <- readRDS(
  file.path(path, "cellchat_object_category_Crohn_s_disease.Rds")
)
cellchat.normal <- readRDS(
  file.path(path, "cellchat_object_category_Normal.Rds")
)

object.list <- list(Crohns = cellchat.crohns, Normal = cellchat.normal)

cellchat_merged <- mergeCellChat(
  object.list, add.names = names(object.list), cell.prefix = TRUE
)

crohns_celltypes = levels(cellchat.crohns@idents)
normal_celltypes = levels(cellchat.normal@idents)

plots <- list()
for (pw in levels(cellchat_merged@idents$joint)) {
  gg1 <- netVisual_bubble(
    object = cellchat_merged,
    sources.use = pw,
    targets.use = NULL,
    sort.by.source = FALSE,
    comparison = c(1, 2),
    title.name = paste("Signalling from", pw),
    angle.x = 45,
    show.legend = TRUE
  )
  plots[[pw]] <- gg1
  h <- max(6, 1.25 + length(levels(gg1@data$interaction_name_2)) * 0.2)
  pdf(file.path(path, paste0(pw, "_LRs_category.pdf")), width = 6, height = h)
  print(gg1)
  dev.off()
}



netVisual_chord_gene(
  object = cellchat.crohns,
  slot.name = "netP",
  lab.cex = 0.5,
  legend.pos.x = 30,
  legend.pos.y = 40,
)


# load files
cellchat.crohns <- readRDS(
  file.path(path, "cellchat_object_final_curated_labels_Crohn_s_disease.Rds")
)
cellchat.normal <- readRDS(
  file.path(path, "cellchat_object_final_curated_labels_Normal.Rds")
)

object.list <- list(Crohns = cellchat.crohns, Normal = cellchat.normal)

cellchat_merged <- mergeCellChat(
  object.list, add.names = names(object.list), cell.prefix = TRUE
)

crohns_celltypes = levels(cellchat.crohns@idents)
normal_celltypes = levels(cellchat.normal@idents)

plots <- list()
for (pw in levels(cellchat_merged@idents$joint)) {
  gg1 <- netVisual_bubble(
    object = cellchat_merged,
    sources.use = pw,
    targets.use = NULL,
    sort.by.source = FALSE,
    comparison = c(1, 2),
    title.name = paste("Signalling from", pw),
    angle.x = 45,
    show.legend = TRUE
  )
  plots[[pw]] <- gg1
  h <- max(6, 1.25 + length(levels(gg1@data$interaction_name_2)) * 0.2)
  pdf(file.path(path, paste0(pw, "_LRs_curated.pdf")), width = 16, height = h)
  print(gg1)
  dev.off()
}



netVisual_chord_gene(
  object = cellchat.crohns,
  slot.name = "netP",
  lab.cex = 0.5,
  legend.pos.x = 30,
  legend.pos.y = 40,
)

