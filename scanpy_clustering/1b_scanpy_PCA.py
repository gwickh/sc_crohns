import os

import anndata as an
import scanpy as sc
from utils.scanpy_pca_utils import (
    pca_dispersion,
    pca_varfeatures,
    plot_elbow_plots,
    plot_pca_loadings,
    plot_variable_feature_plots,
)

SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy"

if os.path.exists(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_reduced.h5ad")):
    print("adata_merged_reduced already created, skipping")
    exit()

# Load existing AnnData object
adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))

# Create output directory
PCA_OUTPUT_PATH = os.path.join(SCANPY_OBJECT_PATH, "PCA_stats")
os.makedirs(PCA_OUTPUT_PATH) if not os.path.exists(PCA_OUTPUT_PATH) else None

# Cell cycle scoring
s_genes = [
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM7",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "CENPU",
    "HELLS",
    "RFC2",
    "POLR1B",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "MRPL36",
    "E2F8",
]

g2m_genes = [
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "PIMREG",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "CDCA3",
    "JPT1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
]

# Parameters for variable gene selection
disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
xmin = 0.1
xmax = 10

# Loop over dispersion cutoffs using mean.var.plot equivalent
pca_dispersion(adata, PCA_OUTPUT_PATH, disps, xmin, xmax, s_genes, g2m_genes)

# Loop over nfeatures with vst method
pca_varfeatures(
    adata,
    PCA_OUTPUT_PATH,
    n_features,
    s_genes,
    g2m_genes,
)

adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_reduced.h5ad"))

# Plot elbow plots
plot_elbow_plots(adata, PCA_OUTPUT_PATH)

# Plot loadings (PC1-3)
n_pcs_to_plot = 5
top_n = 10

plot_pca_loadings(adata, PCA_OUTPUT_PATH)

# Variable feature plots
plot_variable_feature_plots(adata, PCA_OUTPUT_PATH)
