import os
from glob import glob

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation


# Create AnnData objects
def load_count_matrices(path) -> list:
    """
    Load count matrices from given path into AnnData objects.
    """
    matrix_paths = sorted(
        [
            p
            for p in glob(os.path.join(path, "**"), recursive=True)
            if not os.path.isdir(p) and os.path.basename(p).endswith((".h5", ".h5ad"))
        ]
    )
    sample_names = [os.path.basename(p) for p in matrix_paths]

    adata_list = []
    for path, sample in zip(matrix_paths, sample_names):
        if ".h5ad" in sample:
            adata = sc.read_h5ad(path)
            adata.obs["platform"] = "Parse"
        elif ".h5" in sample:
            adata = sc.read_10x_h5(path)
            adata.obs["platform"] = "10X_Chromium"
            adata.var["gene_name"] = adata.var.index
        else:
            raise ValueError(f"Unknown file type: {sample}")
        adata.obs["sample_id"] = sample  # add sample_id column
        adata.var_names_make_unique()  # make var names unique
        adata.obs_names = [
            f"{sample}_{bc}" for bc in adata.obs_names
        ]  # make barcodes unique
        adata.raw = adata  # store raw counts
        adata_list.append(adata)
    return adata_list


# Filter cells with low counts
def filter_low_count_cells(adata_list) -> list:
    """
    Filter cells with low total counts.
    """
    for sample in adata_list:
        sc.pp.filter_cells(sample, min_counts=200)
        sc.pp.filter_cells(sample, min_genes=100)
    return adata_list


# Compute per-sample QC metrics and mitochondrial and ribosomal genes
def compute_qc_metrics(adata_list) -> list:
    """
    Compute QC metrics for each sample in adata_list.
    """
    for sample in adata_list:
        # Identify mitochondrial and ribosomal genes
        sample.var["ribo"] = sample.var["gene_name"].str.startswith(("RPS", "RPL"))
        if sample.var["ribo"].sum() == 0:
            raise ValueError(
                f"No ribosomal genes found in {sample.obs['sample_id'][0]}"
            )

        sample.var["mt"] = sample.var["gene_name"].str.startswith("MT-")
        if sample.var["mt"].sum() == 0:
            raise ValueError(
                f"No mitochondrial genes found in {sample.obs['sample_id'][0]}"
            )

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(
            sample, qc_vars=["mt", "ribo"], inplace=True, log1p=True
        )
    return adata_list


# Plot QC metrics
def qc_plots(stages_dict, qc_dict, outpath) -> None:
    """
    Plot QC metrics across samples
    """
    stages_list = []

    for stage_name, adata_list in stages_dict.items():
        for sample in adata_list:
            sample_id = sample.obs["sample_id"].iloc[0]
            diagnosis = "Crohn's Disease" if "Crohns" in sample_id else "Normal"

            for label, metric in qc_dict.items():
                values = sample.obs[metric].values

                stages_list.append(
                    pd.DataFrame(
                        {
                            "Sample ID": sample_id,
                            "Diagnosis": diagnosis,
                            "Stage": stage_name,
                            "Metric": label,
                            "Value": values,
                            "Log10(cell count)": np.log10(sample.n_obs),
                        }
                    )
                )

    stages_df = pd.concat(stages_list, ignore_index=True)

    for label, metric in qc_dict.items():
        plt.figure(figsize=(22, 6))

        sns.violinplot(
            data=stages_df[stages_df["Metric"] == label],
            x="Sample ID",
            y="Value",
            hue="Stage",
            cut=0,
            inner="quart",
            linewidth=1,
            density_norm="width",
        )
        plt.ylabel(label)
        plt.title(f"{label} across preprocessing stages")
        plt.xticks(rotation=90)
        plt.legend(title="Stage", bbox_to_anchor=(1.05, 1), loc="upper left")

        if outpath:
            fname = f"qc_compare_{label.replace(' ', '_')}.png"
            plt.savefig(os.path.join(outpath, fname), dpi=300, bbox_inches="tight")

        plt.show()

    # Plot log10 cell counts per sample
    cellcount_df = stages_df[
        ["Sample ID", "Stage", "Log10(cell count)"]
    ].drop_duplicates()
    plt.figure(figsize=(22, 6))
    sns.barplot(
        data=cellcount_df,
        x="Sample ID",
        y="Log10(cell count)",
        hue="Stage",
    )
    plt.xticks(rotation=90)
    plt.title("Log10 cell counts")
    plt.savefig(
        os.path.join(outpath, "qc_plot_log_cell_counts_per_sample.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()


# Utility to create MAD mask for a given metric
def create_mad_mask(adata_list, sample, metric, nmads=3) -> np.ndarray:
    platform = sample.obs["platform"].iloc[0]
    platform_values = []

    for s in adata_list:
        if s.obs["platform"].iloc[0] == platform:
            platform_values.append(s.obs[metric].values)

    platform_values = np.concatenate(platform_values)
    sample_values = sample.obs[metric].values

    med = np.median(platform_values)
    mad = median_abs_deviation(platform_values, scale=0.6745)

    # guard for collapse of bounds and dropping all cells when mad is near zero
    if mad < 1e-8:
        mask = np.ones(sample.n_obs, dtype=bool)
    else:
        mask = (sample_values > med - nmads * mad) & (sample_values < med + nmads * mad)
    return mask


# Apply MAD filtering across samples
def mad_filter(adata_list, qc_dict, nmads=3) -> list:
    """
    Compute MAD boolean masks and filter per sample
    """
    adata_list_filtered = []

    for sample in adata_list:
        combined_mask = np.ones(
            sample.n_obs, dtype=bool
        )  # init combined mask to all True

        for label, metric in qc_dict.items():
            mask = create_mad_mask(adata_list, sample, metric)
            combined_mask &= mask  # logical AND to combine masks across metrics

        adata_filtered = sample[combined_mask, :].copy()

        adata_filtered.layers["raw_counts"] = adata_filtered.X.copy()

        adata_list_filtered.append(adata_filtered)

    return adata_list_filtered


# Get QC stats and save to CSV
def obtain_qc_stats(adata_list, adata_list_filtered, qc_dict, outpath) -> None:
    """
    Save QC stats to CSV.
    """
    adata_stats = []

    for sample, filtered in zip(adata_list, adata_list_filtered):
        sample_id = sample.obs["sample_id"].iloc[0]
        diagnosis = "Crohn's Disease" if "Crohns" in sample_id else "Normal"

        row = {
            "Sample ID": sample_id,
            "Diagnosis": diagnosis,
            "Initial number of cells": sample.n_obs,
            "Final number of cells": filtered.n_obs,
            "Proportion retained": filtered.n_obs / sample.n_obs,
        }
        for label, metric in qc_dict.items():
            values = sample.obs[metric].values
            med = np.median(values)
            mad = median_abs_deviation(values, scale=0.6745)
            mask = create_mad_mask(adata_list, sample, metric)

            row.update(
                {
                    f"{label} median": med,
                    f"{label} MAD": mad,
                    f"{label} passed": mask.sum(),
                    f"{label} passed fraction": mask.sum() / sample.n_obs,
                }
            )

        adata_stats.append(row)

    adata_df = pd.DataFrame(adata_stats)
    adata_df.to_csv(os.path.join(outpath, "qc_stats_mad_filtering.csv"), index=False)
