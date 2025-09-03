library(Seurat)
library(Matrix)
library(CellChat)
library(patchwork)
library(igraph)
library(SeuratDisk)

# Vars
path <- "project-area/data/crohns_scrnaseq/scvi_tools_output/"
h5ad_file <- file.path(path, "Integrated_05_label/query_concat_curated.h5ad")
label_column <- c("curated", "category")
species <- "human"
outdir <- file.path(path, "../cellchat_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Convert h5ad to Seurat object
tmp_h5s <- file.path(tempdir(), "tmp.h5seurat")
SeuratDisk::Convert(
  h5ad_file,
  dest = "h5seurat",
  overwrite = TRUE,
  filename = tmp_h5s
)
seu <- SeuratDisk::LoadH5Seurat(tmp_h5s)

DefaultAssay(seu) <- "RNA"
Idents(seu) <- factor(seu@meta.data[[label_column]])

# Run CellChat
data.use <- GetAssayData(seu, assay = "RNA", slot = "data")
meta <- data.frame(group = Idents(seu))
rownames(meta) <- colnames(seu)

cellchat <- createCellChat(
  object = data.use,
  meta = meta,
  group.by = label_column
)

data(CellChatDB.human)
cellchat@DB <- CellChatDB.human
data(PPI.human)
PPI.use <- PPI.human


cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.use)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# save cellchat object and visualisations
saveRDS(cellchat, file.path(outdir, "cellchat_object.rds"))

pdf(file.path(outdir, "bubble_top_pathways.pdf"))
netVisual_bubble(cellchat, remove.isolate = TRUE)
dev.off()

pdf(file.path(outdir, "heatmap_overall.pdf"))
netVisual_heatmap(cellchat, measure = "weight")
dev.off()
