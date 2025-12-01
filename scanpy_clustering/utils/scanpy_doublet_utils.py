import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from scipy.stats import boxcox, norm
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

    elif transformation == "boxcox":
        transformed_scores = np.asarray(boxcox(scores + eps)[0])

    elif transformation == "raw":
        transformed_scores = scores

    else:
        raise ValueError(
            f"Unknown transformation: {transformation},\
            choose from 'log', 'logit', 'probit', 'boxcox', or 'raw'."
        )

    return transformed_scores


# Fit GMM and obtain threshold
def fit_gmm(sample, transformed_scores) -> GaussianMixture:
    """
    Utility to fit a 2-component GMM to scrublet scores.
    """
    # Reshape to column vector for sklearn
    transformed_scores = transformed_scores.reshape(-1, 1)

    # Fit 2-component GMM
    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=0).fit(
        transformed_scores
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


def calculate_point_coords(
    sample,
    gmm_params,
    transformed_scores,
    transformation=None,
    backtransform=False,
    bins=50,
):
    """
    Utility to plot histogram + GMM components + threshold.
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

    if intersection_points.size == 0:
        print(
            f"""No intersection found between GMM components for {sample_id}
             using {transformation} transformation."""
        )

    elif intersection_points.size > 1:
        print(
            f"""Warning: Multiple intersections found for {sample_id} using the
             {transformation} transformation. Using the intersection closest to 0
             as threshold."""
        )

    gaussian_intersection = intersection_points[
        np.argmin(np.abs(x[intersection_points]))
    ]

    # Backtransform to original space
    if backtransform is True:
        if transformation == "log":
            transformed_scores = np.expm1(transformed_scores)
            gaussian_intersection = np.expm1(gaussian_intersection)

            comp1_pdf = comp1_pdf * (np.expm1(x))
            comp2_pdf = comp2_pdf * (np.expm1(x))

        elif transformation == "logit":
            transformed_scores = 1 / (1 + np.exp(-transformed_scores))
            gaussian_intersection = 1 / (1 + np.exp(-gaussian_intersection))

            comp1_pdf = comp1_pdf * (x * (1 - x))
            comp2_pdf = comp2_pdf * (x * (1 - x))

        elif transformation == "probit":
            transformed_scores = norm.cdf(transformed_scores)
            gaussian_intersection = norm.cdf(gaussian_intersection)

            comp1_pdf = comp1_pdf * 1 / norm.pdf(x)
            comp2_pdf = comp2_pdf * 1 / norm.pdf(x)

            x = norm.cdf(x)


    # Components in count scale
    bindwidth = (transformed_scores.max() - transformed_scores.min()) / bins

    comp1 = comp1_pdf * len(transformed_scores) * bindwidth
    comp2 = comp2_pdf * len(transformed_scores) * bindwidth

    return gaussian_intersection, transformed_scores, x, comp1, comp2


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
        threshold = float(x[gaussian_intersection])

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
                outpath, 
                f"{sample_id}_{transformation}_doublets_GM_backtransformed.png"
            ),
            dpi=300,
            bbox_inches="tight",
        )
    else:
        plt.title(f"2-component GMM on {transformation} Scrublet scores of {sample_id}")
        plt.xlabel(f"{transformation}-transformed Scrublet doublet score")
        plt.savefig(
            os.path.join(
                outpath, 
                f"{sample_id}_{transformation}_doublets_GMM.png"
            ),
            dpi=300,
            bbox_inches="tight",
    )
    plt.close()


def calculate_doublet_threshold(
    adata_list_filtered, outpath, transformation, bins=50
) -> None:
    """
    Calculate doublet score threshold based on GMM parameters and plot results
    """
    for sample in adata_list_filtered:
        for transform in transformation:
            transformed_scores = score_preprocessing(sample, transformation=transform)
            gmm = fit_gmm(sample, transformed_scores)
            gmm_params = obtain_gmm_threshold(gmm)

            for b in [True, False]:
                gaussian_intersection, transformed_scores, x, comp1, comp2 = calculate_point_coords(
                    sample,
                    gmm_params,
                    transformed_scores,
                    transformation=transform,
                    bins=bins,
                    backtransform=b,
                )
                plot_gmm(
                    sample,
                    gaussian_intersection,
                    x,
                    comp1,
                    comp2,
                    gmm_params,
                    outpath,
                    transformed_scores,
                    backtransform=b,
                    transformation=transform,
                    bins=bins,
                )
