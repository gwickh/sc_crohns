NO_PAIRS <- FALSE
MIN_PERC <- 0.1
LOGFC <- 0.25
LOGFC_FILT <- 0.5
ADJ_PVAL_FILT <- 0.05
LATENT_VARS <- NULL
TOP_DEG <- 50
GENES_EXCLUDE <- NULL
CORES <- 8

DEA_OUT_DIR <- file.path(SEURAT_OBJECT_PATH, "DEA_output")
dir.create(DEA_OUT_DIR, showWarnings = FALSE)

if (!exists("seurat_object")) {
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_LOC, "seurat_object.Rds"))
}

source(file.path("sc_crohns/clustering", ".Rprofile"))
source(file.path("sc_crohns/dea", "2_dea_functions.R"))

Idents(seurat_object) <- seurat_object@meta.data[["seurat_clusters"]]

genes <- rownames(LayerData(seurat_object, layer = "counts"))
if (!is.null(GENES_EXCLUDE)) {
  genes <- genes[!(genes %in% GENES_EXCLUDE)]
}

if (!NO_PAIRS) {
  ClusterGeneMarkersByPairs(
    object = seurat_object,
    out.dir = DEA_OUT_DIR,
    id = CL_MODE,
    test.use = TEST,
    min.pct = MIN_PERC,
    logFC = LOGFC,
    latent.vars = LATENT_VARS
  )
}

ClusterGeneMarkersVsAll(
  object = seurat_object,
  out.dir = DEA_OUT_DIR,
  id = CL_MODE,
  test.use = TEST,
  min.pct = MIN_PERC,
  logFC = LOGFC,
  latent.vars = LATENT_VARS
)

# ==== Filtering + Plots ====
if (!NO_PAIRS) {
  FilterClusterGeneMarkersPairs(
    object = seurat_object,
    out.dir = DEA_OUT_DIR,
    id = CL_MODE,
    test.use = TEST,
    logFC.filt = LOGFC_FILT,
    adjpval.filt = ADJ_PVAL_FILT,
    num = TOP_DEG,
    genes.use = genes,
    heatmap = TRUE,
    volcano = TRUE
  )
}

FilterClusterGeneMarkersAll(
  object = seurat_object,
  out.dir = DEA_OUT_DIR,
  id = CL_MODE,
  test.use = TEST,
  logFC.filt = LOGFC_FILT,
  adjpval.filt = ADJ_PVAL_FILT,
  num = TOP_DEG,
  genes.use = genes,
  heatmap = TRUE,
  volcano = TRUE
)