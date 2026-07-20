#!/usr/bin/env python3
"""Compute cluster assignment probability for non-negligble predictions."""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

scvi_dir = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/",
)

adata_in = scvi_dir / "query_concat_curated.h5ad"

output_dir = scvi_dir / "celltype_cluster_membership"
output_dir.mkdir(parents=True, exist_ok=True)


embedding_key = "X_embeddings"
label_key = "label_spreading_prediction"
confidence_key = "label_spreading_confidence"

high_conf_threshold = 0.8

n_neighbors_list = [15, 30, 50]
resolution_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]

adata = sc.read_h5ad(adata_in)


# Filter to high-confidence transferred labels

adata = adata[adata.obs[confidence_key] >= high_conf_threshold].copy()

label_key = "label_spreading_prediction"
confidence_key = "label_spreading_confidence"

# Make sure labels are categorical
adata.obs[label_key] = adata.obs[label_key].astype("category")

# Make sure confidence is numeric
adata.obs[confidence_key] = pd.to_numeric(
    adata.obs[confidence_key],
    errors="coerce",
)


# Helper functions
def safe_name(value):
    """Create filename-safe parameter strings."""
    return str(value).replace(".", "p")


def plot_heatmap(
    matrix,
    title,
    output_file,
    cbar_label,
    figsize=(12, 10),
):
    """Plot clustered heatmap."""
    if matrix.empty:
        print(f"Skipping empty matrix: {output_file}")
        return

    g = sns.clustermap(
        matrix,
        cmap="viridis",
        linewidths=0.2,
        figsize=figsize,
        cbar_kws={"label": cbar_label},
        row_cluster=True,
        col_cluster=True,
    )

    g.fig.suptitle(title, y=1.02)

    g.fig.savefig(
        output_file,
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )

    plt.close(g.fig)


def clusters_covering_fraction(row, target_fraction=0.8):
    """
    Compute number of clusters needed to account for target_fraction of a cell type.

    A value of 1 means most cells of that label sit in one cluster.
    A larger value suggests that the label is spread across multiple clusters.
    """
    vals = row.sort_values(ascending=False)
    cumulative = vals.cumsum()
    return int((cumulative < target_fraction).sum() + 1)


def summarise_membership(counts, label_norm, cluster_norm, cluster_key):
    """Produce summary tables of cluster purity and label spread."""
    # Cell type summar (does one label spread across clusters)

    label_summary = pd.DataFrame(index=label_norm.index)

    label_summary["n_cells"] = counts.sum(axis=1)
    label_summary["dominant_cluster"] = label_norm.idxmax(axis=1)
    label_summary["max_cluster_fraction"] = label_norm.max(axis=1)
    label_summary["n_clusters_nonzero"] = (counts > 0).sum(axis=1)
    label_summary["n_clusters_covering_80pct"] = label_norm.apply(
        clusters_covering_fraction,
        axis=1,
        target_fraction=0.8,
    )

    label_summary = label_summary.sort_values(
        ["n_clusters_covering_80pct", "max_cluster_fraction"],
        ascending=[False, True],
    )

    label_summary.to_csv(output_dir / f"{cluster_key}_celltype_split_review.csv")

    # Cluster summary
    non_negligible_threshold = 0.1

    cluster_summary = pd.DataFrame(index=cluster_norm.columns)

    cluster_summary["n_cells"] = counts.sum(axis=0)
    cluster_summary["dominant_celltype"] = cluster_norm.idxmax(axis=0)
    cluster_summary["cluster_purity"] = cluster_norm.max(axis=0)
    cluster_summary["n_celltypes_nonzero"] = (counts > 0).sum(axis=0)
    cluster_summary[f"n_celltypes_above_{int(non_negligible_threshold * 100)}pct"] = (
        cluster_norm >= non_negligible_threshold
    ).sum(axis=0)

    cluster_summary = cluster_summary.sort_values(
        [
            "cluster_purity",
            f"n_celltypes_above_{int(non_negligible_threshold * 100)}pct",
        ],
        ascending=[True, False],
    )

    cluster_summary.to_csv(output_dir / f"{cluster_key}_cluster_mixing_review.csv")

    return label_summary, cluster_summary


def compute_cluster_membership(
    adata,
    n_neighbors_list,
    resolution_list,
    all_label_summaries,
    all_cluster_summaries,
    output_dir=output_dir,
) -> tuple:
    """Compute cluster membership summaries for a range of parameters."""
    for n_neighbors in n_neighbors_list:
        print(f"\nComputing neighbours: n_neighbors={n_neighbors}")

        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            use_rep=embedding_key,
            key_added=f"neighbors_n{n_neighbors}",
        )

        for resolution in resolution_list:
            res_name = safe_name(resolution)

            cluster_key = f"leiden_n{n_neighbors}_r{res_name}"

            print(f"Running Leiden: {cluster_key}")

            sc.tl.leiden(
                adata,
                resolution=resolution,
                neighbors_key=f"neighbors_n{n_neighbors}",
                key_added=cluster_key,
            )

            adata.write_h5ad(output_dir / "query_concat_curated_clustered.h5ad")

            # Raw cell type * cluster counts
            counts = pd.crosstab(
                adata.obs[label_key],
                adata.obs[cluster_key],
            )

            counts.to_csv(output_dir / f"{cluster_key}_celltype_cluster_counts.csv")

            # Probability of cluster given cell type
            label_norm = counts.div(
                counts.sum(axis=1),
                axis=0,
            ).fillna(0)

            label_norm.to_csv(
                output_dir / f"{cluster_key}_p_cluster_given_celltype.csv",
            )

            plot_heatmap(
                matrix=label_norm,
                title=(f"P(cluster | cell type)\n{cluster_key}"),
                output_file=output_dir
                / f"{cluster_key}_p_cluster_given_celltype_heatmap.pdf",
                cbar_label="Proportion of cell type in cluster",
                figsize=(12, 10),
            )

            # Probability of cell type given cluster
            cluster_norm = counts.div(
                counts.sum(axis=0),
                axis=1,
            ).fillna(0)

            cluster_norm.to_csv(
                output_dir / f"{cluster_key}_p_celltype_given_cluster.csv",
            )

            plot_heatmap(
                matrix=cluster_norm,
                title=(f"P(cell type | cluster)\n{cluster_key}"),
                output_file=output_dir
                / f"{cluster_key}_p_celltype_given_cluster_heatmap.pdf",
                cbar_label="Proportion of cluster assigned to cell type",
                figsize=(12, 10),
            )

            # Summary tables
            label_summary, cluster_summary = summarise_membership(
                counts=counts,
                label_norm=label_norm,
                cluster_norm=cluster_norm,
                cluster_key=cluster_key,
            )

            label_summary = label_summary.reset_index().rename(
                columns={label_key: "celltype"},
            )
            label_summary["n_neighbors"] = n_neighbors
            label_summary["resolution"] = resolution
            label_summary["cluster_key"] = cluster_key

            cluster_summary = cluster_summary.reset_index().rename(
                columns={cluster_key: "cluster"},
            )
            cluster_summary["n_neighbors"] = n_neighbors
            cluster_summary["resolution"] = resolution
            cluster_summary["cluster_key"] = cluster_key

            all_label_summaries.append(label_summary)
            all_cluster_summaries.append(cluster_summary)

    return all_label_summaries, all_cluster_summaries


# Save combined summaries across all parameter settings
def main():
    """Run cluster membership analysis."""
    all_label_summaries = []
    all_cluster_summaries = []

    all_label_summaries, all_cluster_summaries = compute_cluster_membership(
        adata=adata,
        n_neighbors_list=n_neighbors_list,
        resolution_list=resolution_list,
        all_label_summaries=all_label_summaries,
        all_cluster_summaries=all_cluster_summaries,
    )

    all_label_summaries = pd.concat(
        all_label_summaries,
        ignore_index=True,
    )

    all_cluster_summaries = pd.concat(
        all_cluster_summaries,
        ignore_index=True,
    )

    all_label_summaries.to_csv(
        output_dir / "all_parameter_celltype_split_review.csv",
        index=False,
    )

    all_cluster_summaries.to_csv(
        output_dir / "all_parameter_cluster_mixing_review.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
