import os
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.mixture import GaussianMixture
from scipy.stats import norm


# Calculate doublet scores using Scrublet
def run_scrublet(adata_list_filtered, expected_doublet_rate=0.06) -> list:
    """
    Run Scrublet on each sample
    """
    for sample in adata_list_filtered:
        sc.pp.scrublet(sample, expected_doublet_rate=expected_doublet_rate)
    return adata_list_filtered


# Fit GMM and obtain threshold
def fit_gmm(sample) -> GaussianMixture:
    """
    Utility to fit a 2-component GMM to scrublet scores.
    """
    scores = sample.obs["doublet_score"].values.reshape(-1, 1)

    # Fit 2-component GMM
    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=0).fit(
        scores
    )

    return gmm


def obtain_gmm_threshold(gmm) -> pd.DataFrame:
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


def plot_gmm(sample, gmm_params, outpath, bins=50):
    """
    Utility to plot histogram + GMM components + threshold.
    """
    scores = sample.obs["doublet_score"].values
    sample_id = sample.obs["sample_id"].iloc[0]

    # X grid
    x = np.linspace(scores.min(), scores.max(), 500)

    # Component densities
    means = gmm_params["mean"].values
    stds = gmm_params["std"].values
    weights = gmm_params["weight"].values

    comp1 = weights[0] * norm.pdf(x, means[0], stds[0])
    comp2 = weights[1] * norm.pdf(x, means[1], stds[1])

    #
    gaussian_intersection = np.where(np.diff(np.sign(comp1 - comp2)) != 0)[0]
    threshold = float(x[gaussian_intersection])

    # Plot
    plt.figure(figsize=(10, 6))
    plt.hist(scores, bins=bins, density=True, alpha=0.5, color="gray", label="Observed")

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

    plt.axvline(
        threshold,
        color="black",
        linestyle=":",
        linewidth=2,
        label=f"Threshold = {threshold:.3f}",
    )

    plt.title(f"GMM fit on Scrublet scores — {sample_id}")
    plt.xlabel("Scrublet doublet score")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        os.path.join(outpath, f"{sample_id}_doublet_GMM.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def calculate_doublet_threshold(adata_list_filtered, outpath, bins=50) -> None:
    """
    Calculate doublet score threshold based on GMM parameters and plot results
    """
    for sample in adata_list_filtered:
        gmm = fit_gmm(sample)
        gmm_params = obtain_gmm_threshold(gmm)
        plot_gmm(sample, gmm_params, outpath, bins)
