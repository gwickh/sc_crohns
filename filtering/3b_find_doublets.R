library(Seurat)
library(DoubletFinder)
library(dplyr)

# Initialize results
results_list <- list()

# Loop through samples
for (s in c("C_03","C_08","C_13","N_04", "N_06", "N_07", "N_12")) {
  while (!is.null(dev.list())) dev.off()
  graphics.off()
  message("Processing sample: ", s)
  
  # Subset
  seu <- subset(seurat_obj, subset = seurat_obj@meta.data$sample_id == s)

  # Parameter sweep to find optimal pK
  sweep.res <- paramSweep(seu, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  best_pK <- as.numeric(as.character(best_pK))
  
  # Estimate expected doublet rate
  exp_rate <- ifelse(s == "C_08", 0.10, 0.076)
  nExp_poi <- round(ncol(seu) * exp_rate)
  
  # Run DoubletFinder
  seu <- doubletFinder(
    seu, 
    PCs = 1:10, 
    pN = 0.25, 
    pK = best_pK, 
    nExp = nExp_poi, 
    reuse.pANN = NULL, 
    sct = FALSE
  )
  
  # Extract doublet classification column (last metadata column)
  df_col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)
  
  # Subset singlets
  seu_singlets <- subset(seu, subset = !!sym(df_col) == "Singlet")
  
  seu$DoubletFinder <- seu@meta.data[[df_col]]

  results_list[[s]] <- seu
}


seurat_doublet <- merge(results_list[[1]], y = results_list[-1])
table(seurat_doublet$DoubletFinder)