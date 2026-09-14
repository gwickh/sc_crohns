library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)

setwd("~/project-area/data/crohns_scrnaseq/10c_14n_analysis/")

sce <- readH5AD("adata_final.h5ad", version = "0.12.3", reader = "R")

seu <- as.Seurat(sce, counts = "X", data = NULL)

seu <- NormalizeData(seu)

saveRDS(seu, file = "cellchat_output/adata_low_conf_hc_celltype_anchor_profiles.Rds")