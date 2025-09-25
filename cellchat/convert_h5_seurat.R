library(zellkonverter)
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")

setwd("~/NetworkDrives/project-area/data/crohns_scrnaseq/scvi_tools_output/Integrated_05_label")

sce <- readH5AD("query_concat_curated.h5ad")
seu <- as.Seurat(sce, counts = "counts", data = "X") 

saveRDS(seu, file = "query_concat_curated.Rds")