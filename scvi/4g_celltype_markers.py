#!/usr/bin/env python3
"""Generate heatmaps of top marker genes per cell type from curated cell identities."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse

mpl.use("Agg")

# Paths
scvi_dir = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/",
)

adata_in = scvi_dir / "query_concat_curated.h5ad"

output_dir = scvi_dir / "celltype_marker_gene_heatmaps"
output_dir.mkdir(parents=True, exist_ok=True)


def filter_cells(
    adata,
    confidence_key,
    label_key,
    high_conf_threshold,
    min_cells_per_type,
) -> sc.AnnData:
    """Filter cells based on label spreading confidence and cell type counts."""
    if label_key not in adata.obs:
        msg = f"{label_key!r} not found in adata.obs"
        raise KeyError(msg)

    if confidence_key not in adata.obs:
        msg = f"{confidence_key!r} not found in adata.obs"
        raise KeyError(msg)

    adata.obs[confidence_key] = pd.to_numeric(
        adata.obs[confidence_key],
        errors="coerce",
    )

    adata.obs[label_key] = adata.obs[label_key].astype(str)

    adata_filt = adata[adata.obs[label_key].notna()].copy()

    adata_filt = adata_filt[adata_filt.obs[label_key] != "Unknown"].copy()

    adata_filt = adata_filt[
        adata_filt.obs[confidence_key] >= high_conf_threshold
    ].copy()

    celltype_counts = adata_filt.obs[label_key].value_counts()

    keep_celltypes = celltype_counts[celltype_counts >= min_cells_per_type].index

    adata_filt = adata_filt[adata_filt.obs[label_key].isin(keep_celltypes)].copy()

    adata_filt.obs[label_key] = adata_filt.obs[label_key].astype("category")
    adata_filt.obs[label_key] = adata_filt.obs[label_key].cat.remove_unused_categories()

    celltype_order = list(adata_filt.obs[label_key].cat.categories)

    print(f"Cells retained: {adata_filt.n_obs:,}")
    print(f"Genes retained: {adata_filt.n_vars:,}")
    print(f"Cell types retained: {len(celltype_order)}")
    print(adata_filt.obs[label_key].value_counts())

    return adata_filt, celltype_order


def select_top_markers(
    deg_df,
    celltype_order,
    top_n_markers,
    min_log2fc,
    max_padj,
) -> pd.DataFrame:
    """Select top N marker genes per cell type by log2FC."""
    marker_tables = []

    for celltype in celltype_order:
        sub = deg_df[deg_df["group"] == celltype].copy()

        # Prefer significant positive markers
        filtered = sub[
            (sub["logfoldchanges"] >= min_log2fc) & (sub["pvals_adj"] <= max_padj)
        ].copy()

        # Fallback if too few genes pass filters
        if filtered.shape[0] < top_n_markers:
            filtered = sub[sub["logfoldchanges"] > 0].copy()

        filtered = filtered.sort_values(
            "logfoldchanges",
            ascending=False,
        )

        filtered = filtered.drop_duplicates("names")

        top = filtered.head(top_n_markers).copy()

        top["marker_for"] = celltype
        top["marker_rank"] = np.arange(1, top.shape[0] + 1)

        marker_tables.append(top)

    marker_df = pd.concat(
        marker_tables,
        ignore_index=True,
    )

    marker_df.to_csv(
        output_dir / "top5_marker_genes_per_celltype_by_log2fc.csv",
        index=False,
    )

    return marker_df


def compute_z_scored_expression(
    adata,
    label_key,
    marker_genes_present,
    celltype_order,
) -> pd.DataFrame:
    """Compute Z-scored expression per gene across cell types."""
    groups = adata.obs[label_key].astype(str).to_numpy()

    expr_matrix = pd.DataFrame(
        index=marker_genes_present,
        columns=celltype_order,
        dtype=float,
    )

    for celltype in celltype_order:
        mask = groups == celltype

        if sparse.issparse(adata.X[mask, :]):
            mean_tp10k = np.asarray(adata.X[mask, :]).mean(axis=0).ravel()
        else:
            mean_tp10k = np.asarray(adata.X[mask, :]).mean(axis=0)

        expr_matrix[celltype] = mean_tp10k + 1

    expr_matrix.to_csv(output_dir / "top_marker_genes_mean_TP10K_plus1_matrix.csv")

    # Z-score TP10K + 1 expression per gene across cell types
    row_means = expr_matrix.mean(axis=1)
    row_stds = expr_matrix.std(axis=1, ddof=0)

    expr_z = expr_matrix.sub(row_means, axis=0).div(row_stds, axis=0)

    expr_z = expr_z.replace([np.inf, -np.inf], np.nan).fillna(0)

    expr_z.to_csv(output_dir / "top_marker_genes_TP10K_plus1_zscore_matrix.csv")

    return expr_z


def plot_heatmap(
    expr_z,
    marker_df,
    celltype_order,
) -> None:
    """Plot heatmap of Z-scored expression of top marker genes per cell type."""
    plot_rows = []
    plot_values = []

    for _, row in marker_df.iterrows():
        gene = row["names"]
        marker_for = row["marker_for"]

        if gene not in expr_z.index:
            continue

        row_label = f"{marker_for} | {gene}"

        plot_rows.append(row_label)
        plot_values.append(expr_z.loc[gene, celltype_order].values)

    plot_matrix = pd.DataFrame(
        plot_values,
        index=plot_rows,
        columns=celltype_order,
    )

    # Plot heatmap
    fig_width = max(10, 0.35 * len(celltype_order))
    fig_height = max(8, 0.22 * plot_matrix.shape[0])

    plt.figure(figsize=(fig_width, fig_height))

    ax = sns.heatmap(
        plot_matrix,
        cmap="vlag",
        center=0,
        linewidths=0.2,
        linecolor="lightgrey",
        cbar_kws={"label": "Gene expression (TP10K + 1) Z-score"},
    )

    ax.set_title(
        "Top 5 marker genes per transferred cell type",
        pad=20,
    )

    ax.set_xlabel("Cell type")
    ax.set_ylabel("Marker gene")

    plt.xticks(
        rotation=45,
        ha="right",
    )

    plt.yticks(
        rotation=0,
        fontsize=6,
    )

    block_sizes = (
        marker_df.groupby("marker_for")
        .size()
        .reindex(celltype_order)
        .fillna(0)
        .astype(int)
    )

    block_ends = np.cumsum(block_sizes.values)

    for y in block_ends[:-1]:
        ax.axhline(
            y,
            color="black",
            linewidth=0.8,
        )

    plt.tight_layout()

    plt.savefig(
        output_dir / "top5_markers_per_celltype_TP10K_plus1_zscore_heatmap.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )

    plt.close()


def main() -> None:
    """Compute Z-scored gene expression and plot top N marker genes per cell type."""
    # User settings
    label_key = "label_spreading_prediction"
    confidence_key = "label_spreading_confidence"

    high_conf_threshold = 0.8
    remove_unknown = True

    top_n_markers = 5
    min_cells_per_type = 20

    target_sum = 1e4

    # Marker filtering
    min_log2fc = 0.25
    max_padj = 0.05

    # Load AnnData
    adata = sc.read_h5ad(adata_in)

    adata.var_names_make_unique()

    # Filter cells
    adata_filt, celltype_order = filter_cells(
        adata,
        confidence_key,
        label_key,
        high_conf_threshold,
        remove_unknown,
        min_cells_per_type,
    )

    # Prepare expression matrix for marker testing.
    adata_de = adata_filt.copy()

    sc.pp.normalize_total(
        adata_de,
        target_sum=target_sum,
    )

    sc.pp.log1p(adata_de)

    # Differential expression: one-vs-rest markers per cell type
    sc.tl.rank_genes_groups(
        adata_de,
        groupby=label_key,
        method="wilcoxon",
        n_genes=adata_de.n_vars,
        use_raw=False,
    )

    deg_df = sc.get.rank_genes_groups_df(
        adata_de,
        group=None,
    )

    deg_df = deg_df.replace([np.inf, -np.inf], np.nan)

    deg_df = deg_df.dropna(subset=["group", "names", "logfoldchanges", "pvals_adj"])

    deg_df["group"] = deg_df["group"].astype(str)

    deg_df.to_csv(
        output_dir / "all_celltype_marker_gene_tests.csv",
        index=False,
    )

    # Select top N markers per cell type by log2FC
    marker_df = select_top_markers(
        deg_df,
        celltype_order,
        top_n_markers,
        min_log2fc,
        max_padj,
    )

    # Prepare expression matrix for plotting
    adata_expr = adata_filt.copy()

    marker_genes_present = [
        gene
        for gene in marker_df[["names"].drop_duplicates().tolist()]
        if gene in adata_expr.var_names
    ]

    adata_expr = adata_expr[:, marker_genes_present].copy()

    sc.pp.normalize_total(
        adata_expr,
        target_sum=target_sum,
    )

    # Compute Z-scored expression per gene across cell types
    expr_z = compute_z_scored_expression(
        adata_expr,
        label_key,
        marker_genes_present,
        celltype_order,
    )

    # Build plotting matrix
    plot_heatmap(expr_z, marker_df, celltype_order)


if __name__ == "__main__":
    main()
