#!/usr/bin/env python3
"""Investigate distribution of low confidence cells."""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

adata_path = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/",
)


def plot_low_conf_umap(adata, conf):
    """Plot full-object UMAP and low-confidence-only UMAP."""
    adata = adata[adata.obs["label_spreading_confidence"] > 0.3]

    adata.obs["assignment_status"] = np.where(
        adata.obs["label_spreading_confidence"] < conf,
        "low_conf",
        "high_conf",
    )

    # low-confidence cells in the full embedding space
    sc.pp.neighbors(
        adata,
        use_rep="X_embeddings",
        n_neighbors=15,
    )

    sc.tl.umap(
        adata,
        min_dist=0.3,
    )

    low_conf_in_full_umap = adata[adata.obs["assignment_status"].eq("low_conf")].copy()

    fig = sc.pl.embedding(
        low_conf_in_full_umap,
        basis="X_umap",
        color=[
            "assignment_status",
            "label_spreading_confidence",
        ],
        ncols=2,
        frameon=True,
        return_fig=True,
        show=False,
    )

    fig.set_size_inches(12, 6)
    fig.savefig(
        adata_path / "umap_low_conf_cells_full_0p8.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )
    plt.close(fig)


def cluster_low_conf_by_high_conf_connectivity(
    adata,
    output_dir,
    embedding_key="X_embeddings",
    confidence_key="label_spreading_confidence",
    prediction_key="label_spreading_prediction",
    platform_key="platform",
    query_platform="Parse",
    confidence_threshold=0.8,
    full_n_neighbors=75,
    high_resolution=0.4,
    low_profile_n_neighbors=15,
    low_profile_resolution=0.2,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # Define low- and high-confidence cells
    # ------------------------------------------------------------

    low_mask = (
        (
            adata.obs[platform_key].astype(str).eq(query_platform)
            & adata.obs[confidence_key].astype(float).lt(confidence_threshold)
        )
        .fillna(False)
        .to_numpy()
    )

    high_mask = ~low_mask

    print(f"Low-confidence cells: {low_mask.sum():,}")
    print(f"High-confidence/reference cells: {high_mask.sum():,}")

    # ------------------------------------------------------------
    # Build full neighbour graph
    # ------------------------------------------------------------

    sc.pp.neighbors(
        adata,
        use_rep=embedding_key,
        n_neighbors=full_n_neighbors,
        key_added="neighbors_full",
    )

    full_graph = adata.obsp["neighbors_full_connectivities"].tocsr()

    # ------------------------------------------------------------
    # Cluster high-confidence cells using full-graph induced edges
    # ------------------------------------------------------------

    high = adata[high_mask].copy()
    high_graph = full_graph[high_mask][:, high_mask].copy()

    sc.tl.leiden(
        high,
        adjacency=high_graph,
        resolution=high_resolution,
        key_added="high_conf_leiden",
    )

    high.obs["high_conf_leiden"] = "HC_" + high.obs["high_conf_leiden"].astype(str)

    adata.obs["high_conf_leiden"] = "not_high_conf"
    adata.obs.loc[high.obs_names, "high_conf_leiden"] = high.obs[
        "high_conf_leiden"
    ].values

    # Label high-confidence clusters by majority cell type
    high_label_table = pd.crosstab(
        high.obs["high_conf_leiden"],
        high.obs[prediction_key],
    )

    high_majority_label = high_label_table.idxmax(axis=1)

    high_cluster_label = {
        cluster: f"{cluster} ({high_majority_label.loc[cluster]})"
        for cluster in high_label_table.index
    }

    # Build low-cell × high-cluster connectivity profiles

    low_high_graph = full_graph[low_mask][:, high_mask].copy()

    high_clusters = high.obs["high_conf_leiden"].astype(str)
    high_codes, high_names = pd.factorize(high_clusters, sort=True)

    high_design = sparse.csr_matrix(
        (
            np.ones(len(high_codes)),
            (np.arange(len(high_codes)), high_codes),
        ),
        shape=(len(high_codes), len(high_names)),
    )

    # rows = low-confidence cells
    # columns = high-confidence clusters
    low_profiles = low_high_graph @ high_design

    low_profiles = low_profiles.toarray()

    low_profiles = pd.DataFrame(
        low_profiles,
        index=adata.obs_names[low_mask],
        columns=[high_cluster_label[cluster] for cluster in high_names],
    )

    # Convert edge weights to fractions per low-confidence cell
    low_profiles = low_profiles.div(
        low_profiles.sum(axis=1).replace(0, np.nan),
        axis=0,
    ).fillna(0)

    low_profiles.to_csv(
        output_dir / "low_cell_to_high_cluster_connectivity_profiles.csv"
    )

    # Cluster low-confidence cells by their high-cluster profiles
    low_profile_adata = ad.AnnData(
        X=low_profiles.to_numpy(),
        obs=adata.obs.loc[low_profiles.index].copy(),
        var=pd.DataFrame(index=low_profiles.columns),
    )

    sc.pp.neighbors(
        low_profile_adata,
        n_neighbors=low_profile_n_neighbors,
        metric="cosine",
    )

    sc.tl.leiden(
        low_profile_adata,
        resolution=low_profile_resolution,
        key_added="low_conf_profile_leiden",
    )

    low_profile_adata.obs["low_conf_profile_leiden"] = "LC_" + low_profile_adata.obs[
        "low_conf_profile_leiden"
    ].astype(str)

    adata.obs["low_conf_profile_leiden"] = "not_low_conf"
    adata.obs.loc[
        low_profile_adata.obs_names,
        "low_conf_profile_leiden",
    ] = low_profile_adata.obs["low_conf_profile_leiden"].values

    # Summarise each LC cluster by high-confidence similarity

    lc_to_hc = (
        low_profiles.assign(
            low_conf_profile_leiden=low_profile_adata.obs[
                "low_conf_profile_leiden"
            ].values
        )
        .groupby("low_conf_profile_leiden")
        .mean()
    )

    lc_to_hc.to_csv(
        output_dir / "low_conf_cluster_to_high_conf_cluster_connectivity.csv"
    )

    lc_summary = pd.DataFrame(
        {
            "n_cells": low_profile_adata.obs["low_conf_profile_leiden"].value_counts(),
            "best_matching_high_conf_cluster": lc_to_hc.idxmax(axis=1),
            "best_matching_fraction": lc_to_hc.max(axis=1),
        }
    )

    lc_summary.to_csv(output_dir / "low_conf_cluster_summary.csv")

    adata.write_h5ad(output_dir / "adata_low_conf_profile_clusters.h5ad")

    return low_profile_adata, low_profiles, lc_to_hc, lc_summary


def save_clustered_heatmap(df, path, title, cmap="viridis", center=None):
    if df.shape[0] < 2 or df.shape[1] < 2:
        fig, ax = plt.subplots(figsize=(10, 4))
        sns.heatmap(df, cmap=cmap, center=center, ax=ax)
        ax.set_title(title)
        fig.savefig(path, bbox_inches="tight", dpi=300)
        plt.close(fig)
    else:
        grid = sns.clustermap(
            df,
            cmap=cmap,
            center=center,
            linewidths=0.2,
            figsize=(
                max(8, 0.35 * df.shape[1]),
                max(5, 0.35 * df.shape[0]),
            ),
            cbar_kws={"label": "Connectivity fraction"},
        )
        grid.fig.suptitle(title, y=1.02)
        grid.fig.savefig(path, bbox_inches="tight", dpi=300)
        plt.close(grid.fig)


def main() -> None:
    """Visualise and cluster low confidence cells."""
    adata = ad.read_h5ad(
        adata_path
        / "sysvi_tuning"
        / "c561826c_sysvi_label_spreading_alpha_0.2_n_5.h5ad",
    )

    # plot_low_conf_umap(adata, 0.8)

    (low_profile_adata, low_profiles, lc_to_hc, lc_summary) = (
        cluster_low_conf_by_high_conf_connectivity(
            adata,
            output_dir=adata_path / "low_confidence_profile_clustering",
            confidence_threshold=0.8,
            full_n_neighbors=75,
            high_resolution=0.4,
            low_profile_n_neighbors=15,
            low_profile_resolution=0.1,
        )
    )

    # Heatmap LC cluster × HC cluster connectivity
    save_clustered_heatmap(
        lc_to_hc,
        adata_path
        / "low_confidence_profile_clustering"
        / "heatmap_low_conf_cluster_to_high_conf_cluster_connectivity.pdf",
        title="Low-confidence clusters by high-confidence cluster connectivity",
        cmap="viridis",
    )

    # Heatmap column-scaled LC cluster × HC cluster connectivity
    lc_to_hc_z = (
        lc_to_hc.sub(lc_to_hc.mean(axis=0), axis=1)
        .div(
            lc_to_hc.std(axis=0).replace(0, np.nan),
            axis=1,
        )
        .fillna(0)
    )

    save_clustered_heatmap(
        lc_to_hc_z,
        adata_path
        / "low_confidence_profile_clustering"
        / "heatmap_low_conf_cluster_to_high_conf_cluster_connectivity_zscore.pdf",
        title="Column-scaled LC-to-HC connectivity",
        cmap="vlag",
        center=0,
    )

    # Heatmap low-confidence cells ordered by LC cluster
    low_cell_profiles = low_profiles.copy()
    low_cell_profiles["low_conf_profile_leiden"] = low_profile_adata.obs[
        "low_conf_profile_leiden",
    ].values

    low_cell_profiles = low_cell_profiles.sort_values("low_conf_profile_leiden")

    row_groups = low_cell_profiles.pop("low_conf_profile_leiden")

    plot_profiles = low_cell_profiles.assign(
        low_conf_profile_leiden=row_groups.values
    ).groupby("low_conf_profile_leiden", group_keys=False)

    plot_row_groups = plot_profiles.pop("low_conf_profile_leiden")

    grid = sns.clustermap(
        plot_profiles,
        row_cluster=False,
        col_cluster=True,
        cmap="viridis",
        linewidths=0,
        yticklabels=False,
        figsize=(max(8, 0.35 * plot_profiles.shape[1]), 10),
        cbar_kws={"label": "Connectivity fraction"},
    )

    grid.fig.suptitle(
        "Low-confidence cells by high-confidence cluster connectivity",
        y=1.02,
    )

    grid.fig.savefig(
        adata_path
        / "low_confidence_profile_clustering"
        / "heatmap_low_conf_cells_to_high_conf_clusters.png",
        bbox_inches="tight",
        dpi=300,
    )

    grid.fig.savefig(
        adata_path
        / "low_confidence_profile_clustering"
        / "heatmap_low_conf_cells_to_high_conf_clusters.pdf",
        bbox_inches="tight",
    )

    plt.close(grid.fig)


if __name__ == "__main__":
    main()
