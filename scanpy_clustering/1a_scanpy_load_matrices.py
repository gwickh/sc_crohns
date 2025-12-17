import os

import anndata as ad
import scanpy as sc
from utils.scanpy_doublet_utils import calculate_doublet_threshold, run_scrublet
from utils.scanpy_qc_utils import (
    compute_qc_metrics,
    filter_low_count_cells,
    load_count_matrices,
    mad_filter,
    obtain_qc_stats,
    qc_plots,
)

os.chdir("/users/yep25yan/dev")

# Set paths
MATRIX_DIR = "project-area/data/crohns_scrnaseq/10c_14n_analysis/crohns_samples"

OUTPATH = os.path.join(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis", "scanpy", "qc_stats"
)

os.makedirs(OUTPATH) if not os.path.exists(OUTPATH) else None

# Check if adata_merged already exists, if true then end script
if os.path.exists(os.path.join(OUTPATH, "adata_merged.h5ad")):
    print("adata_merged already created, skipping")
    exit()

# Load matrix paths and sample names
adata_list = load_count_matrices(MATRIX_DIR)

adata_list_raw = [adata.copy() for adata in adata_list]

# Run preprocessing steps
qc_dict = {
    "Number of genes": "n_genes_by_counts",
    "Log10 total read count": "log1p_total_counts",
    "Percent mitochondrial reads": "pct_counts_mt",
    "Percent ribosomal reads": "pct_counts_ribo",
}

adata_list = filter_low_count_cells(adata_list)

adata_list_raw = compute_qc_metrics(adata_list_raw)
adata_list = compute_qc_metrics(adata_list)

adata_list_filtered = mad_filter(adata_list, qc_dict)

adata_list_filtered = run_scrublet(adata_list_filtered)

adata_list_filtered = calculate_doublet_threshold(
    adata_list_filtered, outpath=OUTPATH, transformation=["probit"], bins=50
)

obtain_qc_stats(adata_list, adata_list_filtered, qc_dict, OUTPATH)

stages_dict = {
    "Raw": adata_list_raw,
    "Background filtering": adata_list,
    "Platform-level MAD filtering": adata_list_filtered,
}

qc_plots(stages_dict, qc_dict, OUTPATH)

# # Merge filtered AnnData objects
adata = ad.concat(
    adata_list_filtered,
    label="sample_id",
    keys=[a.obs["sample_id"].iloc[0] for a in adata_list_filtered],
    join="outer",
    merge="unique",
)

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# Save merged AnnData object
adata.write_h5ad(os.path.join(OUTPATH, "adata_merged.h5ad"))
