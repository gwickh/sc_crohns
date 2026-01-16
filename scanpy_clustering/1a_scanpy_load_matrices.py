# from utils.scanpy_doublet_utils import calculate_doublet_threshold, run_scrublet
# from utils.scanpy_qc_utils import (
#     compute_qc_metrics,
#     filter_low_count_cells,
#     load_count_matrices,
#     mad_filter,
#     obtain_qc_stats,
#     qc_plots,
# )
import os
from glob import glob

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation, norm
from sklearn.mixture import GaussianMixture


# Calculate doublet scores using Scrublet
def run_scrublet(adata_list_filtered, expected_doublet_rate=0.06) -> list:
    """
    Run Scrublet on each sample
    """
    for sample in adata_list_filtered:
        sc.pp.scrublet(sample, expected_doublet_rate=expected_doublet_rate)
    return adata_list_filtered


# prepare scores for GMM fitting
def score_preprocessing(sample, transformation="raw") -> np.ndarray:
    """
    Utility to preprocess scrublet scores for GMM fitting.
    """

    scores = sample.obs["doublet_score"].values
    # transform scores to fit GMM better
    eps = 1e-10
    p = np.clip(scores, eps, 1 - eps)

    if transformation == "log":
        transformed_scores = np.log1p(scores + eps)

    elif transformation == "logit":
        transformed_scores = np.log(p / (1 - p))

    elif transformation == "probit":
        transformed_scores = norm.ppf(scores + eps)

    elif transformation == "raw":
        transformed_scores = scores

    else:
        raise ValueError(
            f"Unknown transformation: {transformation},\
            choose from 'log', 'logit', 'probit', or 'raw'."
        )

    return transformed_scores


# Fit GMM and obtain threshold
def fit_gmm(transformed_scores, tail_fraction=0.05) -> GaussianMixture:
    """
    Utility to fit a 2-component GMM to scrublet scores.
    """
    # Reshape to column vector for sklearn
    transformed_scores = transformed_scores.reshape(-1, 1)

    # Initialize means based on tail points
    q_tail = np.quantile(transformed_scores, 1 - tail_fraction)
    tail_points = transformed_scores[transformed_scores >= q_tail]

    high_init = np.median(tail_points)
    low_init = np.quantile(transformed_scores, 0.2)

    means_init = np.array([[low_init], [high_init]])
    means_init = np.sort(means_init, axis=0)

    # Fit 2-component GMM
    gmm = GaussianMixture(
        n_components=2,
        covariance_type="full",
        random_state=0,
        reg_covar=1e-6,
        n_init=1,
        init_params="kmeans",
        means_init=means_init,
    ).fit(transformed_scores)

    return gmm


def get_gmm_params(gmm) -> pd.DataFrame:
    """
    Create dataframe with Gaussian parameters and decision threshold.
    """

    # Extract means, variances, weights
    means = gmm.means_.flatten()
    vars_ = gmm.covariances_.flatten()
    stds = np.sqrt(vars_)
    weights = gmm.weights_.flatten()

    # Sort components by mean
    idx = np.argsort(means)

    # Create DataFrame of GMM parameters
    gmm_params = pd.DataFrame(
        {
            "mean": means[idx],
            "std": stds[idx],
            "weight": weights[idx],
        },
        index=["singlet", "doublet"],
    )

    return gmm_params


def backtransform_threshold(
    gaussian_intersection, transformation, transformed_scores, x, comp1_pdf, comp2_pdf
) -> None:
    """
    Backtransform threshold to original score space.
    """
    if transformation == "log":
        transformed_scores = np.expm1(transformed_scores)
        gaussian_intersection = np.expm1(gaussian_intersection)

        comp1_pdf = comp1_pdf * 1 / (1 + np.expm1(x))
        comp2_pdf = comp2_pdf * 1 / (1 + np.expm1(x))

        x = np.expm1(x)

    elif transformation == "logit":
        transformed_scores = 1 / (1 + np.exp(-transformed_scores))
        gaussian_intersection = 1 / (1 + np.exp(-gaussian_intersection))

        term = 1 + np.exp(-x)
        comp1_pdf = comp1_pdf * 1 / (1 / term) * (1 - 1 / term)
        comp2_pdf = comp2_pdf * 1 / (1 / term) * (1 - 1 / term)

        x = 1 / (1 + np.exp(-x))

    elif transformation == "probit":
        transformed_scores = norm.cdf(transformed_scores)
        gaussian_intersection = norm.cdf(gaussian_intersection)

        comp1_pdf = comp1_pdf * 1 / norm.pdf(x)
        comp2_pdf = comp2_pdf * 1 / norm.pdf(x)

        x = norm.cdf(x)


def compute_threshold(
    sample,
    gmm_params,
    transformed_scores,
    transformation=None,
    backtransform=False,
    bins=50,
) -> tuple:
    """
    Utility to compute GMM threshold based on intersection of components.
    """
    sample_id = sample.obs["sample_id"].iloc[0]

    # X grid
    x = np.linspace(transformed_scores.min(), transformed_scores.max(), 500)

    # Component densities
    means = gmm_params["mean"].values
    stds = gmm_params["std"].values
    weights = gmm_params["weight"].values

    comp1_pdf = weights[0] * norm.pdf(x, means[0], stds[0])
    comp2_pdf = weights[1] * norm.pdf(x, means[1], stds[1])

    # Calculate threshold at intersection of components
    intersection_points = np.where(np.diff(np.sign(comp1_pdf - comp2_pdf)) != 0)[0]
    print(intersection_points)

    if intersection_points.size == 0:
        print(
            f"""No intersection found between GMM components for {sample_id}
             using {transformation} transformation."""
        )
        return None, transformed_scores, x, comp1_pdf, comp2_pdf

    elif intersection_points.size > 1:
        print(
            f"""Warning: Multiple intersections found for {sample_id} using the
             {transformation} transformation. Using the intersection closest to 0
             as threshold."""
        )

    gaussian_intersection = x[
        intersection_points[np.argmin(np.abs(x[intersection_points]))]
    ]

    # Backtransform to original space
    if backtransform is True:
        backtransform_threshold(
            gaussian_intersection,
            transformation,
            transformed_scores,
            x,
            comp1_pdf,
            comp2_pdf,
        )

    # Components in count scale
    binwidth = (transformed_scores.max() - transformed_scores.min()) / bins

    comp1 = comp1_pdf * len(transformed_scores) * binwidth
    comp2 = comp2_pdf * len(transformed_scores) * binwidth

    return gaussian_intersection, transformed_scores, x, comp1, comp2


def call_doublet(sample, gaussian_intersection) -> sc.AnnData:
    """
    Utility to call doublets based on GMM threshold.
    """
    if gaussian_intersection is None:
        raise ValueError("No threshold found; cannot call doublets.")
    else:
        sample_id = sample.obs["sample_id"].iloc[0]
        threshold = float(gaussian_intersection)

        sample.obs["predicted_doublet"] = np.where(
            sample.obs["doublet_score"] >= threshold,
            "doublet",
            "singlet",
        )
        print(
            f"Doublet threshold for {sample_id}: {threshold:.3f}",
            f"Predicted doublet fraction for {sample_id}: {(sample.obs['predicted_doublet'] == 'doublet').sum()}"
            f"Predicted singlet fraction for {sample_id}: {(sample.obs['predicted_doublet'] == 'singlet').sum()}",
        )

    return sample


def plot_gmm(
    sample,
    gaussian_intersection,
    x,
    comp1,
    comp2,
    gmm_params,
    outpath,
    transformed_scores,
    transformation=None,
    backtransform=False,
    bins=50,
) -> None:
    """
    Utility to plot histogram + GMM components + threshold.
    """
    sample_id = sample.obs["sample_id"].iloc[0]
    means = gmm_params["mean"].values
    stds = gmm_params["std"].values

    # Plot
    plt.figure(figsize=(10, 6))
    plt.hist(
        transformed_scores,
        bins=bins,
        density=False,
        alpha=0.5,
        color="gray",
        label="Observed",
    )

    plt.plot(
        x,
        comp1,
        label=f"Singlet component (μ={means[0]:.2f}, σ={stds[0]:.2f})",
        color="blue",
    )

    plt.plot(
        x,
        comp2,
        label=f"Doublet component (μ={means[1]:.2f}, σ={stds[1]:.2f})",
        color="red",
    )

    plt.plot(
        x, comp1 + comp2, label="Gaussian Mixture model", color="black", linestyle="--"
    )

    # Add threshold line if intersection exists
    if len(np.atleast_1d(gaussian_intersection)) == 1:
        threshold = float(gaussian_intersection)

        plt.axvline(
            threshold,
            color="black",
            linestyle=":",
            linewidth=2,
            label=f"Threshold = {threshold:.3f}",
        )

    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()

    if backtransform is True:
        plt.title(f"2-component GMM on backtransformed Scrublet scores of {sample_id}")
        plt.xlabel("Scrublet doublet score")
        plt.savefig(
            os.path.join(
                outpath, f"{sample_id}_{transformation}_doublets_GM_backtransformed.png"
            ),
            dpi=300,
            bbox_inches="tight",
        )
    else:
        plt.title(f"2-component GMM on {transformation} Scrublet scores of {sample_id}")
        plt.xlabel(f"{transformation}-transformed Scrublet doublet score")
        plt.savefig(
            os.path.join(outpath, f"{sample_id}_{transformation}_doublets_GMM.png"),
            dpi=300,
            bbox_inches="tight",
        )
    plt.close()


def plot_doublet_proportions(adata_list_filtered, thresholds, outpath) -> None:
    """
    Utility to plot doublet proportions based on GMM threshold.
    """

    results = []

    for adata in adata_list_filtered:
        sample_id = adata.obs["sample_id"].iloc[0]

        if sample_id not in thresholds:
            print(f"No threshold for {sample_id}; skipping.")
            continue

        threshold = thresholds[sample_id]

        scores = adata.obs["doublet_score"].values
        frac = np.mean(scores >= threshold)

        results.append({"Sample ID": sample_id, "Doublet fraction": frac})

    df = pd.DataFrame(results)

    # Plot
    plt.figure(figsize=(12, 6))
    sns.barplot(data=df, x="Sample ID", y="Doublet fraction", color="red")
    plt.xticks(rotation=90)
    plt.ylim(0, df["Doublet fraction"].max() * 1.1)
    plt.title("Doublet Fraction per Sample")
    plt.ylabel("Fraction of cells ≥ threshold")

    plt.tight_layout()
    plt.savefig(
        os.path.join(outpath, "doublet_fraction_barplot.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def calculate_doublet_threshold(
    adata_list_filtered, outpath, transformation, bins=50
) -> list:
    """
    Calculate doublet score threshold based on GMM parameters and plot results
    """

    thresholds = {}

    for sample in adata_list_filtered:
        for transform in transformation:
            transformed_scores = score_preprocessing(sample, transformation=transform)
            gmm = fit_gmm(transformed_scores)
            gmm_params = get_gmm_params(gmm)

            # Compute threshold in transformed space
            t_gaussian_intersection, t_transformed_scores, t_x, t_comp1, t_comp2 = (
                compute_threshold(
                    sample,
                    gmm_params,
                    transformed_scores,
                    transformation=transform,
                    bins=bins,
                    backtransform=False,
                )
            )

            plot_gmm(
                sample,
                t_gaussian_intersection,
                t_x,
                t_comp1,
                t_comp2,
                gmm_params,
                outpath,
                t_transformed_scores,
                backtransform=False,
                transformation=transform,
                bins=bins,
            )

            # Compute threshold in original space
            (
                bt_gaussian_intersection,
                bt_transformed_scores,
                bt_x,
                bt_comp1,
                bt_comp2,
            ) = compute_threshold(
                sample,
                gmm_params,
                transformed_scores,
                transformation=transform,
                bins=bins,
                backtransform=True,
            )

            plot_gmm(
                sample,
                bt_gaussian_intersection,
                bt_x,
                bt_comp1,
                bt_comp2,
                gmm_params,
                outpath,
                bt_transformed_scores,
                backtransform=True,
                transformation=transform,
                bins=bins,
            )

            sample = call_doublet(sample, bt_gaussian_intersection)
            thresholds[sample.obs["sample_id"].iloc[0]] = bt_gaussian_intersection

    plot_doublet_proportions(adata_list_filtered, thresholds, outpath)

    return adata_list_filtered


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
def mad_filter(adata_list, qc_dict) -> list:
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
adata.write_h5ad(os.path.join(os.path.dirname(OUTPATH), "adata_merged.h5ad"))
