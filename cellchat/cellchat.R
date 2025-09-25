library(Seurat)
library(Matrix)
library(CellChat)
library(patchwork)
library(igraph)
options(stringsAsFactors = FALSE)

# Vars
path <- "project-area/data/crohns_scrnaseq/scvi_tools_output/"
h5seurat_path <- file.path(path, "Integrated_05_label/query_concat_curated.h5seurat")
label_column <- c("curated", "category")

h5seurat_file <- readH5Seurat(h5seurat_path)

data.input <- h5seurat_file[["RNA"]]$data # normalized data matrix
labels <- Idents(h5seurat_file)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellChat <- createCellChat(object = h5seurat_file, group.by = "ident", assay = "RNA")

# # Run CellChat
# data.use <- GetAssayData(seu, assay = "RNA", slot = "data")
# meta <- data.frame(group = Idents(seu))
# rownames(meta) <- colnames(seu)

# cellchat <- createCellChat(
#   object = data.use,
#   meta = meta,
#   group.by = label_column
# )

# data(CellChatDB.human)
# cellchat@DB <- CellChatDB.human
# data(PPI.human)
# PPI.use <- PPI.human


# cellchat <- subsetData(cellchat)
# future::plan("multisession", workers = 4)

# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.use)
# cellchat <- computeCommunProb(cellchat)
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# cellchat <- computeCommunProbPathway(cellchat)
# cellchat <- aggregateNet(cellchat)
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# # save cellchat object and visualisations
# saveRDS(cellchat, file.path(outdir, "cellchat_object.rds"))

# pdf(file.path(outdir, "bubble_top_pathways.pdf"))
# netVisual_bubble(cellchat, remove.isolate = TRUE)
# dev.off()

# pdf(file.path(outdir, "heatmap_overall.pdf"))
# netVisual_heatmap(cellchat, measure = "weight")
# dev.off()
