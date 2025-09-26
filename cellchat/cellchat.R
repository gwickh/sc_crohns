library(Seurat)
library(Matrix)
library(CellChat)
library(patchwork)
library(igraph)
library(future)

options(future.globals.maxSize = 16 * 1024^3) 
options(stringsAsFactors = FALSE)

# Vars
path <- "project-area/data/crohns_scrnaseq/scvi_tools_output/"
rds_path <- file.path(path, "Integrated_05_label/query_concat_curated.Rds")
label_column <- c("curated", "category")
outdir <- "project-area/data/crohns_scrnaseq/cellchat_output/"

# load file
seurat_file <- readRDS(rds_path)

# check data structure
print("class")
class(seurat_file) 

print("Assays")
Assays(seurat_file)

print("DefaultAssay")
DefaultAssay(seurat_file) 

print("slotNames")
slotNames(seurat_file[["originalexp"]])

# rename assays
seurat_file <- RenameAssays(seurat_file, originalexp = "RNA")
DefaultAssay(seurat_file) <- "RNA"
colnames(seurat_file@meta.data)[colnames(seurat_file@meta.data) == "sample_id"] <- "samples"

# create cellchat object
cellChat <- createCellChat(object = seurat_file, group.by = "category", assay = "RNA")

# set DB
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellChat@DB <- CellChatDB.use
cellChat <- subsetData(cellChat)

# Identify overexpressed genes and interactions
future::plan("multisession", workers = 4)

cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
cellChat <- smoothData(cellChat, adj = PPI.human)

# Compute the communication probability
cellChat <- computeCommunProb(cellChat, raw.use = FALSE, type = "triMean")
cellChat <- filterCommunication(cellChat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "net")

# save cellchat object and visualisations
saveRDS(cellChat, file.path(outdir, "cellchat_object.rds"))




