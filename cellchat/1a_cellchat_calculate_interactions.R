library(Seurat)
library(Matrix)
library(CellChat)
library(patchwork)
library(igraph)
library(future)

options(stringsAsFactors = FALSE)
# options(future.globals.maxSize = 16 * 1024^3) 
# future::plan("multisession", workers = 4)

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

diagnosis_type <- unique(seurat_file@meta.data$Diagnosis)

# create cellchat object
run_cellchat = function(celltypes) {
    for (samp in diagnosis_type) {
        seurat_file_subset <- subset(seurat_file, subset = Diagnosis == samp)
        cellChat <- createCellChat(object = seurat_file_subset, group.by = celltypes, assay = "RNA")

        # set DB
        CellChatDB <- CellChatDB.human 
        CellChatDB.use <- subsetDB(CellChatDB)
        cellChat@DB <- CellChatDB.use
        cellChat <- subsetData(cellChat)

        # Identify overexpressed genes and interactions
        cellChat <- identifyOverExpressedGenes(cellChat)
        cellChat <- identifyOverExpressedInteractions(cellChat)

        # Compute the communication probability
        cellChat <- computeCommunProb(cellChat, trim = NULL, type = "triMean")
        cellChat <- filterCommunication(cellChat, min.cells = 10)

        # Infer the cell-cell communication at a signaling pathway level
        cellChat <- computeCommunProbPathway(cellChat)
        cellChat <- aggregateNet(cellChat)
        cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")
        cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "net")

        # save cellchat object and visualisations
        saveRDS(cellChat, file.path(outdir, paste0("cellchat_object_", celltypes, "_", samp, ".Rds")))
    }
}

run_cellchat(celltypes = "category")
run_cellchat(celltypes = "curated")

