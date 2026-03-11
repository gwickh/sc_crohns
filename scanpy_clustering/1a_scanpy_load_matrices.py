import os

import anndata as ad
import pandas as pd
import scanpy as sc
from utils.scanpy_doublet_utils import calculate_doublet_threshold, run_scrublet
from utils.scanpy_qc_utils import (
    compute_qc_metrics,
    concatenate_adata,
    filter_low_count_cells,
    load_count_matrices,
    mad_filter,
    obtain_qc_stats,
    qc_plots,
)

# set pandas string handling to use builtin str type, not pyarrow to avoid anndata IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

# Set paths
MATRIX_DIR = "project-area/data/crohns_scrnaseq/10c_14n_analysis/crohns_samples"

OUTPATH = os.path.join(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis", "scanpy", "qc_stats"
)

os.makedirs(OUTPATH, exist_ok=True)

# Run preprocessing steps
qc_dict = {
    "Number of genes": "n_genes_by_counts",
    "Log10 total read count": "log1p_total_counts",
    "Percent mitochondrial reads": "pct_counts_mt",
    "Percent ribosomal reads": "pct_counts_ribo",
}


def main() -> None:
    # Check if adata_merged already exists, if true then end script
    if os.path.exists(
        os.path.join(os.path.dirname(OUTPATH), "adata_merged.h5ad")
    ) and os.path.exists(os.path.join((OUTPATH), "qc_stats_mad_filtering.csv")):
        print("adata_merged already created, skipping")
        return

    # Load matrix paths and sample names
    adata_list = load_count_matrices(MATRIX_DIR)

    adata_list_raw = [adata.copy() for adata in adata_list]
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

    # Merge filtered AnnData objects
    adata = concatenate_adata(adata_list_filtered)

    print(adata.var_names[:10])
    print(adata.var.index[:10])
    print(adata.var["ensembl_id"].head(10))

    # Store raw counts in layers before normalization
    adata.layers["counts"] = adata.X.copy()

    # Normalize and log-transform the data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Save merged AnnData object
    adata.write_h5ad(os.path.join(os.path.dirname(OUTPATH), "adata_merged.h5ad"))


if __name__ == "__main__":
    main()
