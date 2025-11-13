import os
from glob import glob
import anndata as ad
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set paths
SCRIPT_DIR = "sc_crohns/scanpy_clustering"
MATRIX_DIR = "project-area/data/crohns_scrnaseq/10c_14n_analysis/crohns_samples"

OUTPATH = os.path.join("project-area/data/crohns_scrnaseq/10c_14n_analysis", "scanpy", "qc_stats")
os.makedirs(OUTPATH, exist_ok=True) if not os.path.exists(OUTPATH) else None

# Check if adata_merged already exists, if true then end script
if os.path.exists(os.path.join(OUTPATH, "adata_merged.h5ad")):
    print("adata_merged already created, skipping")
    exit()

# Load matrix paths and sample names
matrix_paths = sorted(
    [p for p in glob(os.path.join(MATRIX_DIR, "**"), recursive=True) if not os.path.isdir(p)]
)
sample_names = [os.path.basename(p) for p in matrix_paths]

# Create AnnData objects
adata_list = []
for path, sample in zip(matrix_paths, sample_names):
    if ".h5ad" in sample:
        ad = sc.read_h5ad(path)
        ad.obs["platform"] = "Parse"  
    elif ".h5" in sample:
        ad = sc.read_10x_h5(path)
        ad.var["platform"] = "10X_Chromium"  
        ad.var["gene_name"] = ad.var.index
    else:
        raise ValueError(f"Unknown file type: {sample}")
    ad.obs["sample_id"] = sample                                # add sample_id column
    ad.var_names_make_unique()                                  # make var names unique
    ad.obs_names = [f"{sample}_{bc}" for bc in ad.obs_names]    # make barcodes unique
    ad.raw = ad                                                 # store raw counts
    adata_list.append(ad)

# Compute per-sample QC metrics and mitochondrial and ribosomal genes
qc_list = []
for sample in adata_list:
    # Identify mitochondrial and ribosomal genes
    sample.var["ribo"] = sample.var["gene_name"].str.startswith(("RPS", "RPL"))
    if sample.var["ribo"].sum() == 0:
        raise ValueError(f"No ribosomal genes found in {sample.obs['sample_id'][0]}")
    
    sample.var["mt"] = sample.var["gene_name"].str.startswith("MT-")
    if sample.var["mt"].sum() == 0:
        raise ValueError(f"No mitochondrial genes found in {sample.obs['sample_id'][0]}")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        sample,
        qc_vars=["mt", "ribo"],
        inplace=True,
        log1p=True
    )

# Visualize QC metrics across samples
def qc_violin_plot(adata_list,label, metric):
    qc_list = []
    for sample in adata_list:
        qc_list.append(pd.DataFrame({
            "Sample ID":    sample.obs["sample_id"].iloc[0],
            label:          sample.obs[metric]
        }))
    df = pd.concat(qc_list)
    plt.figure(figsize=(12, 5)) 
    sns.violinplot(data=df, x="Sample ID", y=label, cut=0, inner="quart")
    plt.xticks(rotation=90)

    plt.savefig(os.path.join(OUTPATH, f"qc_plot_{metric}_raw.png"), dpi=300, bbox_inches="tight")

qc_dict = {
    "Number of genes": "n_genes_by_counts", 
    "Total read count": "total_counts",
    "Log10 total read count": "log1p_total_counts",
    "Percent mitochondrial reads": "pct_counts_mt",
    "Percent ribosomal reads": "pct_counts_ribo",
}

for label, metric in qc_dict.items():
    qc_violin_plot(adata_list, label, metric)


# # Set static variables
# MIN_CELLS = 3
# MIN_FEATURES = 200
# VARS_TO_REGRESS = ["S_score", "G2M_score"]

# # Filter low-quality cells
# adata = adata[(adata.obs.n_genes_by_counts > 500) & 
#               (adata.obs.n_genes_by_counts < 6000) & 
#               (adata.obs.pct_counts_mt < 10), :]

# # Normalize
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# # Cell cycle scoring
# # Load cell cycle gene lists (use Seurat's list)
# s_genes = [
#     "MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG", 
#     "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", 
#     "RFC2", "POLR1B", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", 
#     "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", 
#     "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", 
#     "MRPL36", "E2F8"
# ]

# g2m_genes = [
#     "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", 
#     "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", 
#     "TACC3", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", 
#     "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", 
#     "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", 
#     "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", 
#     "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", 
#     "G2E3", "GAS2L3", "CBX5", "CENPA"
# ]
# sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
# adata.obs['Phase'] = adata.obs['phase']

# sc.pp.regress_out(adata, keys=VARS_TO_REGRESS)
# sc.pp.scale(adata, max_value=10)

# # Save the object
# adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))
