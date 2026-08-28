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
    purity_threshold=0.8,
):
    """Compute LC-HC connectivities and cluster."""
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

    # ------------------------------------------------------------
    # Label high-confidence clusters by majority cell-type purity
    # ------------------------------------------------------------

    high_label_table = pd.crosstab(
        high.obs["high_conf_leiden"],
        high.obs[prediction_key],
    ).astype(int)

    high_cluster_size = high_label_table.sum(axis=1)
    high_majority_label = high_label_table.idxmax(axis=1)
    high_majority_count = high_label_table.max(axis=1)

    high_purity = high_majority_count / high_cluster_size
    high_is_pure = high_purity.ge(purity_threshold)

    high_cluster_label = {
        cluster: (
            f"{cluster} ({high_majority_label.loc[cluster]})"
            if high_is_pure.loc[cluster]
            else f"{cluster} (impure)"
        )
        for cluster in high_label_table.index
    }

    high_cluster_stats = pd.DataFrame(
        {
            "cluster_label": pd.Series(high_cluster_label),
            "n_cells": high_cluster_size,
            "majority_cell_type": high_majority_label,
            "majority_count": high_majority_count,
            "purity": high_purity,
            "is_pure": high_is_pure,
        }
    )

    high_cluster_stats.to_csv(
        output_dir / "high_conf_cluster_purity_stats.csv",
    )

    high.obs["high_conf_cluster_label"] = high.obs["high_conf_leiden"].map(
        high_cluster_label
    )

    adata.obs["high_conf_cluster_label"] = "not_high_conf"
    adata.obs.loc[
        high.obs_names,
        "high_conf_cluster_label",
    ] = high.obs["high_conf_cluster_label"].values

    # ------------------------------------------------------------
    # For each cell type, compute how many cells sit in pure HC clusters
    # ------------------------------------------------------------

    pure_clusters = high_is_pure[high_is_pure].index
    impure_clusters = high_is_pure[~high_is_pure].index

    hc_celltype_records = []

    for cell_type in high_label_table.columns:
        total_cells = high_label_table[cell_type].sum()

        n_in_pure_clusters = high_label_table.loc[
            pure_clusters,
            cell_type,
        ].sum()

        n_in_impure_clusters = high_label_table.loc[
            impure_clusters,
            cell_type,
        ].sum()

        own_pure_clusters = [
            cluster
            for cluster in pure_clusters
            if high_majority_label.loc[cluster] == cell_type
        ]

        n_in_own_pure_clusters = (
            high_label_table.loc[
                own_pure_clusters,
                cell_type,
            ].sum()
            if len(own_pure_clusters) > 0
            else 0
        )

        hc_celltype_records.append(
            {
                "cell_type": cell_type,
                "n_cells": int(total_cells),
                "n_cells_in_pure_clusters": int(n_in_pure_clusters),
                "pct_cells_in_pure_clusters": (
                    100 * n_in_pure_clusters / total_cells
                    if total_cells > 0
                    else np.nan
                ),
                "n_cells_in_impure_clusters": int(n_in_impure_clusters),
                "pct_cells_in_impure_clusters": (
                    100 * n_in_impure_clusters / total_cells
                    if total_cells > 0
                    else np.nan
                ),
                "n_cells_in_own_pure_majority_clusters": int(n_in_own_pure_clusters),
                "pct_cells_in_own_pure_majority_clusters": (
                    100 * n_in_own_pure_clusters / total_cells
                    if total_cells > 0
                    else np.nan
                ),
                "own_pure_majority_clusters": ";".join(own_pure_clusters),
            }
        )

    hc_celltype_stats = pd.DataFrame(hc_celltype_records).sort_values(
        "pct_cells_in_own_pure_majority_clusters",
        ascending=False,
    )

    hc_celltype_stats.to_csv(
        output_dir / "high_conf_celltype_purity_stats.csv",
        index=False,
    )

    # Build low-cell x high-cluster connectivity profiles

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
        output_dir / "low_cell_to_high_cluster_connectivity_profiles.csv",
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
    adata.obs.loc[low_profile_adata.obs_names, "low_conf_profile_leiden"] = (
        low_profile_adata.obs["low_conf_profile_leiden"].values
    )

    # Summarise each LC cluster by high-confidence similarity

    lc_to_hc = (
        low_profiles.assign(
            low_conf_profile_leiden=low_profile_adata.obs[
                "low_conf_profile_leiden"
            ].values,
        )
        .groupby("low_conf_profile_leiden")
        .mean()
    )

    lc_to_hc.to_csv(
        output_dir / "low_conf_cluster_to_high_conf_cluster_connectivity.csv",
    )

    adata.write_h5ad(output_dir / "adata_low_conf_profile_clusters.h5ad")

    return (
        low_profile_adata,
        low_profiles,
        lc_to_hc,
    )


def cluster_low_conf_by_hc_celltype_anchors(
    adata,
    output_dir,
    embedding_key="X_embeddings",
    confidence_key="label_spreading_confidence",
    prediction_key="label_spreading_prediction",
    platform_key="platform",
    query_platform="Parse",
    confidence_threshold=0.8,
    high_confidence_threshold=0.8,
    full_n_neighbors=75,
    low_profile_n_neighbors=15,
    low_profile_resolution=0.6,
    exclude_labels=("Unknown", "nan", "None"),
):
    """
    Assign low-confidence cells by connectivity to high-confidence cell-type anchors.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define low-confidence Parse cells
    is_parse = adata.obs[platform_key].astype(str).eq(query_platform)
    confidence = adata.obs[confidence_key].astype(float)
    low_mask = (is_parse & confidence.lt(confidence_threshold)).fillna(False).to_numpy()

    # Define high-confidence anchor cells
    labels = adata.obs[prediction_key].astype(str)
    valid_label = adata.obs[prediction_key].notna() & ~labels.isin(exclude_labels)
    high_mask = (
        (
            valid_label
            & (~is_parse | confidence.ge(high_confidence_threshold))
            & ~pd.Series(low_mask, index=adata.obs_names)
        )
        .fillna(False)
        .to_numpy()
    )

    if low_mask.sum() == 0:
        raise ValueError("No low-confidence cells found.")

    if high_mask.sum() == 0:
        raise ValueError("No high-confidence anchor cells found.")

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
    # Build low-cell x high-cell-type connectivity profiles
    # ------------------------------------------------------------

    low_high_graph = full_graph[low_mask][:, high_mask].copy()

    high_labels = adata.obs.loc[high_mask, prediction_key].astype(str)

    anchor_codes, anchor_names = pd.factorize(
        high_labels,
        sort=True,
    )

    anchor_design = sparse.csr_matrix(
        (
            np.ones(len(anchor_codes)),
            (np.arange(len(anchor_codes)), anchor_codes),
        ),
        shape=(len(anchor_codes), len(anchor_names)),
    )

    # rows = low-confidence cells
    # columns = high-confidence cell-type anchors
    low_profiles = low_high_graph @ anchor_design
    low_profiles = low_profiles.toarray()

    low_profiles = pd.DataFrame(
        low_profiles,
        index=adata.obs_names[low_mask],
        columns=anchor_names.astype(str),
    )

    # Convert raw graph connectivity to fractions per LC cell
    low_profiles = low_profiles.div(
        low_profiles.sum(axis=1).replace(0, np.nan),
        axis=0,
    ).fillna(0)

    low_profiles.to_csv(output_dir / "low_cell_to_hc_celltype_anchor_profiles.csv")

    # ------------------------------------------------------------
    # Cluster LC cells by high-confidence anchor profile
    # ------------------------------------------------------------

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

    # ------------------------------------------------------------
    # Summarise LC clusters by HC cell-type anchor connectivity
    # ------------------------------------------------------------

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
        output_dir / "low_conf_cluster_to_hc_celltype_anchor_connectivity.csv"
    )

    # ------------------------------------------------------------
    # Assignment summary
    # ------------------------------------------------------------

    lc_summary = pd.DataFrame(index=lc_to_hc.index)

    lc_summary["n_cells"] = (
        low_profile_adata.obs["low_conf_profile_leiden"]
        .value_counts()
        .reindex(lc_to_hc.index)
    )

    lc_summary["best_matching_hc_celltype"] = lc_to_hc.idxmax(axis=1)
    lc_summary["max_connectivity"] = lc_to_hc.max(axis=1)

    lc_summary["second_matching_hc_celltype"] = lc_to_hc.apply(
        lambda row: row.sort_values(ascending=False).index[1]
        if row.shape[0] > 1
        else np.nan,
        axis=1,
    )

    lc_summary["second_connectivity"] = lc_to_hc.apply(
        lambda row: row.sort_values(ascending=False).iloc[1]
        if row.shape[0] > 1
        else np.nan,
        axis=1,
    )

    lc_summary["connectivity_margin"] = (
        lc_summary["max_connectivity"] - lc_summary["second_connectivity"]
    )

    lc_summary.to_csv(output_dir / "low_conf_cluster_hc_celltype_anchor_summary.csv")

    # ------------------------------------------------------------
    # Per-cell assignment summary
    # ------------------------------------------------------------

    cell_summary = pd.DataFrame(index=low_profiles.index)

    cell_summary["best_matching_hc_celltype"] = low_profiles.idxmax(axis=1)
    cell_summary["max_connectivity"] = low_profiles.max(axis=1)

    cell_summary["low_conf_profile_leiden"] = low_profile_adata.obs[
        "low_conf_profile_leiden"
    ].values

    cell_summary.to_csv(output_dir / "low_conf_cell_hc_celltype_anchor_summary.csv")

    adata.write_h5ad(output_dir / "adata_low_conf_hc_celltype_anchor_profiles.h5ad")

    return adata, low_profile_adata, low_profiles, lc_to_hc, lc_summary, cell_summary


def summarise_lc_anchor_runs(
    output_base,
    summary_filename="low_conf_cluster_hc_celltype_anchor_summary.csv",
    recursive=True,
    max_threshold=0.8,
    relaxed_threshold=0.5,
    ratio_threshold=3.0,
    margin_threshold=0.5,
    output_csv=None,
):
    """
    Loop through output directories, load LC cluster summary files, and summarise
    how many LC cells fall in confidently assignable clusters.
    """

    output_base = Path(output_base)

    pattern = f"**/{summary_filename}" if recursive else summary_filename
    summary_files = sorted(output_base.glob(pattern))

    if len(summary_files) == 0:
        raise FileNotFoundError(
            f"No files named {summary_filename!r} found under {output_base}"
        )

    records = []

    for summary_file in summary_files:
        df = pd.read_csv(summary_file, index_col=0)

        required = {
            "n_cells",
            "max_connectivity",
            "second_connectivity",
            "connectivity_margin",
        }

        missing = required.difference(df.columns)

        if missing:
            raise KeyError(
                f"{summary_file} is missing required columns: {sorted(missing)}"
            )

        df["n_cells"] = pd.to_numeric(df["n_cells"])
        df["max_connectivity"] = pd.to_numeric(df["max_connectivity"])
        df["second_connectivity"] = pd.to_numeric(df["second_connectivity"])
        df["connectivity_margin"] = pd.to_numeric(df["connectivity_margin"])

        strict_mask = df["max_connectivity"].gt(max_threshold)

        relaxed_specific_mask = df["max_connectivity"].gt(relaxed_threshold) & (
            df["max_connectivity"].gt(ratio_threshold * df["second_connectivity"])
            | df["connectivity_margin"].gt(margin_threshold)
        )

        n_lc_cells_strict = df.loc[strict_mask, "n_cells"].sum()

        n_lc_cells_relaxed_specific = df.loc[
            relaxed_specific_mask,
            "n_cells",
        ].sum()

        union_mask = strict_mask | relaxed_specific_mask

        n_lc_cells_union = df.loc[
            union_mask,
            "n_cells",
        ].sum()

        total_lc_cells = df["n_cells"].sum()

        records.append(
            {
                "run_dir": summary_file.parent.name,
                "n_lc_clusters": df.shape[0],
                "total_lc_cells": int(total_lc_cells),
                "connectivity_exceeding_0.8": int(n_lc_cells_strict),
                "connectivity_exceeding_0.5_first-second_ratio_exceeding_3": int(
                    n_lc_cells_relaxed_specific
                ),
                "n_lc_cells_union": int(n_lc_cells_union),
                "pct_connectivity_exceeding_0.8": (
                    100 * n_lc_cells_strict / total_lc_cells
                    if total_lc_cells > 0
                    else np.nan
                ),
                "pct_exceeding_0.5_first-second_ratio_exceeding_3": (
                    100 * n_lc_cells_relaxed_specific / total_lc_cells
                    if total_lc_cells > 0
                    else np.nan
                ),
                "pct_lc_cells_union": (
                    100 * n_lc_cells_union / total_lc_cells
                    if total_lc_cells > 0
                    else np.nan
                ),
            }
        )

    summary = pd.DataFrame(records).sort_values(
        ["n_lc_cells_union"],
        ascending=False,
    )

    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(output_csv, index=False)

    return summary


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
    # adata = ad.read_h5ad(
    #     adata_path
    #     / "sysvi_tuning"
    #     / "c561826c_sysvi_label_spreading_alpha_0.2_n_5.h5ad",
    # )

    # adata_ref = ad.read_h5ad(
    #     adata_path / "query_concat_curated_clustered_reassigned.h5ad"
    # )

    # adata.obs["curated_labels"] = (
    #     adata_ref.obs["curated_labels"]
    #     .reindex(adata.obs_names)
    #     .astype("object")
    #     .fillna("low_conf")
    #     .astype("category")
    # )

    # for neighbors_grid in (10, 15, 30, 50, 75):
    #     for resolution_grid in (0.3, 0.4, 0.6, 0.8):
    # (low_profile_adata, low_profiles, lc_to_hc) = (
    #     cluster_low_conf_by_high_conf_connectivity(
    #         adata,
    #         output_dir=adata_path
    #         / "low_confidence_profile_clustering"
    #         / f"HC_n{neighbors_grid}_res{resolution_grid}",
    #         confidence_threshold=0.8,
    #         full_n_neighbors=neighbors_grid,
    #         high_resolution=resolution_grid,
    #         low_profile_n_neighbors=neighbors_grid,
    #         low_profile_resolution=resolution_grid,
    #     )
    # )

    # (
    #     adata,
    #     low_profile_adata,
    #     low_profiles,
    #     lc_to_hc,
    #     lc_summary,
    #     cell_summary,
    # ) = cluster_low_conf_by_hc_celltype_anchors(
    #     adata,
    #     output_dir=adata_path
    #     / "low_confidence_celltype_anchor_clustering"
    #     / f"HC_n{neighbors_grid}_res{resolution_grid}",
    #     confidence_threshold=0.8,
    #     high_confidence_threshold=0.9,
    #     full_n_neighbors=75,
    #     low_profile_n_neighbors=neighbors_grid,
    #     low_profile_resolution=resolution_grid,
    # )

    # # Heatmap LC cluster X HC cluster connectivity
    # save_clustered_heatmap(
    #     lc_to_hc,
    #     adata_path
    #     / "low_confidence_celltype_anchor_clustering"
    #     / f"HC_n{neighbors_grid}_res{resolution_grid}"
    #     / "heatmap_low_conf_cluster_to_high_conf_cluster_connectivity.pdf",
    #     title="Low-confidence clusters by high-confidence cluster connectivity",
    #     cmap="viridis",
    # )

    # # Heatmap column-scaled LC cluster x HC cluster connectivity
    # lc_to_hc_z = (
    #     lc_to_hc.sub(lc_to_hc.mean(axis=0), axis=1)
    #     .div(
    #         lc_to_hc.std(axis=0).replace(0, np.nan),
    #         axis=1,
    #     )
    #     .fillna(0)
    # )

    # save_clustered_heatmap(
    #     lc_to_hc_z,
    #     adata_path
    #     / "low_confidence_celltype_anchor_clustering"
    #     / f"HC_n{neighbors_grid}_res{resolution_grid}"
    #     / "heatmap_low_conf_cluster_to_high_conf_cluster_connectivity_zscore.pdf",
    #     title="Column-scaled LC-to-HC connectivity",
    #     cmap="vlag",
    #     center=0,
    # )

    # # Heatmap low-confidence cells ordered by LC cluster
    # low_cell_profiles = low_profiles.copy()

    # low_cell_profiles["low_conf_profile_leiden"] = low_profile_adata.obs[
    #     "low_conf_profile_leiden"
    # ].values

    # low_cell_profiles = low_cell_profiles.sort_values("low_conf_profile_leiden")

    # # Optional downsample per LC cluster for readability
    # max_cells_per_cluster = 300

    # plot_profiles = low_cell_profiles.groupby(
    #     "low_conf_profile_leiden", group_keys=False
    # ).apply(
    #     lambda x: x.sample(
    #         min(len(x), max_cells_per_cluster),
    #         random_state=0,
    #     )
    # )

    # plot_row_groups = plot_profiles["low_conf_profile_leiden"].copy()

    # plot_profiles = plot_profiles.drop(columns="low_conf_profile_leiden")

    # grid = sns.clustermap(
    #     plot_profiles,
    #     row_cluster=False,
    #     col_cluster=True,
    #     cmap="viridis",
    #     linewidths=0,
    #     yticklabels=False,
    #     figsize=(max(8, 0.35 * plot_profiles.shape[1]), 10),
    #     cbar_kws={"label": "Connectivity fraction"},
    # )

    # grid.fig.suptitle(
    #     "Low-confidence cells by high-confidence cluster connectivity",
    #     y=1.02,
    # )

    # grid.fig.savefig(
    #     adata_path
    #     / "low_confidence_celltype_anchor_clustering"
    #     / f"HC_n{neighbors_grid}_res{resolution_grid}"
    #     / "heatmap_low_conf_cells_to_high_conf_clusters.pdf",
    #     bbox_inches="tight",
    # )

    # plt.close(grid.fig)

    summary = summarise_lc_anchor_runs(
        output_base=adata_path / "low_confidence_celltype_anchor_clustering",
        output_csv=(
            adata_path
            / "low_confidence_celltype_anchor_clustering"
            / "lc_anchor_connectivity_summary.csv"
        ),
    )


if __name__ == "__main__":
    main()
