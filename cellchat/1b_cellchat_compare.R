library(CellChat)
library(Seurat)
library(Matrix)
library(patchwork)
library(igraph)
library(future)

options(stringsAsFactors = FALSE)

# Vars
path <- "project-area/data/crohns_scrnaseq/cellchat_output"

# load files
cellchat.crohns <- readRDS(file.path(path, "cellchat_object_category_Crohn's Disease.Rds"))
cellchat.normal <- readRDS(file.path(path, "cellchat_object_category_Normal.Rds"))

crohns_celltypes = levels(cellchat.crohns@idents)
normal_celltypes = levels(cellchat.normal@idents)

object.list <- list(Crohns = cellchat.crohns, Normal = cellchat.normal)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# total interactions barplot
gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# heatmap of interactions
gg1 <- netVisual_heatmap(cellchat_merged)
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
gg1 + gg2