library(CellChat)
library(Seurat)
library(Matrix)
library(patchwork)
library(igraph)
library(future)

options(stringsAsFactors = FALSE)

# Vars
path <- "dev/project-area/data/crohns_scrnaseq/cellchat_output"

# load files
cellchat.crohns <- readRDS(file.path(path, "cellchat_object_category_Crohn's Disease.Rds"))
cellchat.normal <- readRDS(file.path(path, "cellchat_object_category_Normal.Rds"))

crohns_celltypes = levels(cellchat.crohns@idents)
normal_celltypes = levels(cellchat.normal@idents)

object.list <- list(Crohns = cellchat.crohns, Normal = cellchat.normal)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# total interactions barplot
gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(2,1))
gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(2,1), measure = "weight")
gg1 + gg2

# heatmap of interactions
# crohns interactions
gg1 <- netVisual_heatmap(cellchat.crohns)
gg2 <- netVisual_heatmap(cellchat.crohns, measure = "weight")
draw(gg1+gg2)

decorate_heatmap_body("Number of interactions", {
  grid.text("Destination (Receiver)", x = unit(19, "lines"), y = unit(-7, "lines"), gp = gpar(fontsize = 10))
})

# normal interactions
gg1 <- netVisual_heatmap(cellchat.normal)
gg2 <- netVisual_heatmap(cellchat.normal, measure = "weight")
draw(gg1+gg2)

decorate_heatmap_body("Number of interactions", {
  grid.text("Destination (Receiver)", x = unit(19, "lines"), y = unit(-7, "lines"), gp = gpar(fontsize = 10))
})

# differential interactions
gg1 <- netVisual_heatmap(cellchat_merged, comparison = c(2, 1))
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight", comparison = c(2, 1))
draw(gg1+gg2)

decorate_heatmap_body("Relative values", {
  grid.text("Destination (Receiver)", x = unit(-1, "lines"), y = unit(-7, "lines"), gp = gpar(fontsize = 10))
})


