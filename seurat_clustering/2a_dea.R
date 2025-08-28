PAIRED <- FALSE
MIN_PERC <- 0.1
LOGFC <- 0.25
LOGFC_FILT <- 0.5
ADJ_PVAL_FILT <- 0.05
LATENT_VARS <- NULL
TOP_DEG <- 50
GENES_EXCLUDE <- NULL
CORES <- 8
CL_MODE <- "clusters_pca_vst_top_5000_k_30_res_1.1"
TEST <- "wilcox"

source(file.path("sc_crohns/seurat_clustering", ".Rprofile"))
source(file.path("sc_crohns/seurat_clustering", "2_dea_functions.R"))

DEA_OUT_DIR <- file.path(SEURAT_OBJECT_PATH, "DEA_output")
dir.create(DEA_OUT_DIR, showWarnings = FALSE)

if (!exists("seurat_object")) {
  print("Loading seurat object")  # print to stdout
  seurat_object <- readRDS(file.path(SEURAT_OBJECT_PATH, "seurat_object.Rds"))
}

# Check if CL_MODE exists in seurat_object@meta.data
print(CL_MODE %in% colnames(seurat_object@meta.data))

print("Merging data layers")  # print to stdout
seurat_object <- MergeDataLayersSparse(seurat_object = seurat_object)
Idents(seurat_object) <- seurat_object@meta.data[[CL_MODE]]

print("Identifying genes")  # print to stdout
genes <- rownames(LayerData(seurat_object, layer = "merged.data"))
if (!is.null(GENES_EXCLUDE)) {
  print("Excluding genes")  # print to stdout
  genes <- genes[!(genes %in% GENES_EXCLUDE)]
}

if (PAIRED) {
  print("Performing pairwise DEA")  # print to stdout
  ClusterGeneMarkersByPairs(
    seurat_object = seurat_object,
    out.dir = DEA_OUT_DIR,
    id = CL_MODE,
    test.use = TEST,
    min.pct = MIN_PERC,
    logFC = LOGFC,
    latent.vars = LATENT_VARS
  )
}

print("Performing one vs all DEA")  # print to stdout
ClusterGeneMarkersVsAll(
  seurat_object = seurat_object,
  out.dir = DEA_OUT_DIR,
  id = CL_MODE,
  test.use = TEST,
  min.pct = MIN_PERC,
  logFC = LOGFC,
  latent.vars = LATENT_VARS
)

# ==== Filtering + Plots ====
if (PAIRED) {
  print("Filtering pairwise DEA")  # print to stdout
  FilterClusterGeneMarkersPairs(
    seurat_object = seurat_object,
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

print("Filtering one vs all DEA")  # print to stdout
FilterClusterGeneMarkersAll(
  seurat_object = seurat_object,
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