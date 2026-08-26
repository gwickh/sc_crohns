#!/usr/bin/env python3
"""Identify marker genes and generate volcano plots for cluster-cell types."""

from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from adjustText import adjust_text
from matplotlib.lines import Line2D

SCVI_PATH = Path("project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output")
INPUT_FILE = SCVI_PATH / "query_concat_curated_clustered.h5ad"

OUTPUT_DIR = SCVI_PATH / "candidate_subcluster_markers"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


candidate_clusters = {
    "B cells": ["8", "20", "21"],
    "ADAMDEC1+ stromal": ["17", "19"],
    "Macrophages": ["34", "14", "2", "3"],
    "Plasma cells": ["4", "16"],
}

GROUPBY = "temp_subtypes"

# Significance thresholds
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0
MIN_IN_GROUP_FRACTION = 0.2
MIN_FRACTION_DIFFERENCE = 0.2

# Groups smaller than this are excluded
MIN_CELLS_PER_GROUP = 20

# Maximum number of significant genes annotated
MAX_LABELS = 50


def safe_filename(value: object) -> str:
    """Convert a group label into a safe filename."""
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    return value.strip("_") or "unnamed_group"


def clean_name(value: str) -> str:
    """Convert a label into a safe subtype name."""
    return str(value).strip().replace(" ", "_").replace("+", "pos").replace("/", "_")


def create_temp_subtypes(
    adata,
    celltype_obs,
    cluster_obs,
    sample_id,
    candidate_clusters,
) -> ad.AnnData:
    """Create temp_subtypes obs for candidate celltype-cluster comparisons."""
    for col in [celltype_obs, cluster_obs, sample_id]:
        if col not in adata.obs.columns:
            msg = f"Column '{col}' not found in AnnData.obs."
            raise ValueError(msg)

        adata.obs[col] = adata.obs[col].astype("string")

    temp_subtypes = pd.Series(
        "Not_candidate",
        index=adata.obs_names,
        dtype="string",
    )

    assignment_summary = []

    for celltype, candidate_cluster in candidate_clusters.items():
        candidate_cluster_str = [str(cluster) for cluster in candidate_cluster]

        for cluster in candidate_cluster_str:
            mask = adata.obs[celltype_obs].eq(celltype) & adata.obs[cluster_obs].eq(
                cluster,
            )

            subtype = f"{clean_name(celltype)}_c{cluster}"

            temp_subtypes.loc[mask] = subtype

            assignment_summary.append(
                {
                    "cell_type": celltype,
                    "cluster": cluster,
                    "temp_subtype": subtype,
                    "n_cells": int(mask.sum()),
                    "n_samples": int(adata.obs.loc[mask, sample_id].nunique()),
                },
            )

    adata.obs["temp_subtypes"] = pd.Categorical(temp_subtypes)

    pd.DataFrame(assignment_summary).to_csv(
        OUTPUT_DIR / "candidate_subcluster_cell_counts.csv",
        index=False,
    )

    return adata


def validate_anndata(
    adata: ad.AnnData,
    groupby: str,
    min_cells_per_group: int,
) -> ad.AnnData:
    """Validate the AnnData object and remove groups with too few cells."""
    if groupby not in adata.obs.columns:
        raise KeyError(
            f"{groupby!r} is not present in adata.obs. "
            f"Available columns: {adata.obs.columns.tolist()}",
        )

    adata = adata[adata.obs[groupby].notna()].copy()

    adata.obs[groupby] = adata.obs[groupby].astype(str)

    group_counts = adata.obs[groupby].value_counts()

    excluded_groups = group_counts.loc[group_counts < min_cells_per_group]

    if not excluded_groups.empty:
        print(f"\nExcluding groups with fewer than {min_cells_per_group} cells:")
        print(excluded_groups.to_string())

    valid_groups = group_counts.loc[group_counts >= min_cells_per_group].index

    adata = adata[adata.obs[groupby].isin(valid_groups)].copy()

    adata.obs[groupby] = pd.Categorical(adata.obs[groupby])

    print("\nCells per retained cluster-cell type:")
    print(adata.obs[groupby].value_counts().sort_index().to_string())

    return adata


def ensembl_id_to_gene_name(adata, dict_map):
    """Check var_names are gene names."""
    ensembl_fraction = np.mean(adata.var_names.astype(str).str.startswith("ENSG"))

    if ensembl_fraction > 0.5:
        print(f"{ensembl_fraction * 100}% of adata.var_names appear to be Ensembl IDs.")

        mapping_df = pd.read_csv(SCVI_PATH / dict_map, sep="\t", header=None)

        id_to_name = dict(zip(mapping_df[0], mapping_df[1], strict=True))

        mapped_names = pd.Series(
            adata.var_names.map(id_to_name),
            index=adata.var_names,
        )

        original_names = pd.Series(
            adata.var_names,
            index=adata.var_names,
        )

        adata.var["gene_name"] = mapped_names.fillna(original_names).to_numpy()

        adata.var_names = adata.var["gene_name"].astype(str)

    return adata


def extract_marker_results(
    adata: ad.AnnData,
    groupby: str,
    padj_threshold: float,
    log2fc_threshold: float,
    min_expression_fraction: float,
    min_fraction_difference: float,
) -> pd.DataFrame:
    """Extract and classify marker results for every group."""
    frames = []

    for group in adata.obs[groupby].cat.categories:
        result = sc.get.rank_genes_groups_df(
            adata,
            group=group,
            key="cluster_celltype_markers",
        )

        result = result.rename(
            columns={
                "names": "gene",
                "scores": "test_statistic",
                "logfoldchanges": "log2fc",
                "pvals": "pvalue",
                "pvals_adj": "padj",
                "pct_nz_group": "fraction_in_group",
                "pct_nz_reference": "fraction_in_reference",
            },
        )

        result.insert(
            0,
            groupby,
            str(group),
        )

        numeric_columns = [
            "test_statistic",
            "log2fc",
            "pvalue",
            "padj",
            "fraction_in_group",
            "fraction_in_reference",
        ]

        for column in numeric_columns:
            if column in result.columns:
                result[column] = pd.to_numeric(
                    result[column],
                    errors="coerce",
                )

        # Positive values mean more cells express the gene in the target group.
        # Negative values mean more cells express it in the reference.
        result["expression_fraction_difference"] = (
            result["fraction_in_group"] - result["fraction_in_reference"]
        )

        result["direction"] = "nonsignificant"

        upregulated = (
            (result["padj"] < padj_threshold)
            & (result["log2fc"] >= log2fc_threshold)
            & (result["fraction_in_group"] >= min_expression_fraction)
            & (result["expression_fraction_difference"] >= min_fraction_difference)
        )

        downregulated = (
            (result["padj"] < padj_threshold)
            & (result["log2fc"] <= -log2fc_threshold)
            & (result["fraction_in_reference"] >= min_expression_fraction)
            & (result["expression_fraction_difference"] <= -min_fraction_difference)
        )

        result.loc[
            upregulated,
            "direction",
        ] = "upregulated"

        result.loc[
            downregulated,
            "direction",
        ] = "downregulated"

        result["significant"] = result["direction"] != "nonsignificant"

        frames.append(result)

    return pd.concat(
        frames,
        ignore_index=True,
    )


def select_labels(
    results: pd.DataFrame,
    max_labels: int,
) -> pd.DataFrame:
    """Select significant genes for annotation."""
    significant = results.loc[results["significant"]].copy()

    if significant.empty or max_labels == 0:
        return significant

    significant["label_priority"] = significant["log2fc"].abs() * -np.log10(
        significant["padj"].clip(lower=np.finfo(float).tiny),
    )

    n_up = max_labels // 2
    n_down = max_labels - n_up

    up = significant.loc[significant["direction"] == "upregulated"].nlargest(
        n_up,
        "label_priority",
    )

    down = significant.loc[significant["direction"] == "downregulated"].nlargest(
        n_down,
        "label_priority",
    )

    selected = pd.concat(
        [up, down],
        ignore_index=True,
    )

    if len(selected) < max_labels:
        remaining = significant.loc[
            ~significant["gene"].isin(selected["gene"])
        ].nlargest(
            max_labels - len(selected),
            "label_priority",
        )

        selected = pd.concat(
            [selected, remaining],
            ignore_index=True,
        )

    return selected


def plot_volcano(
    results: pd.DataFrame,
    group_name: str,
    output_path: Path,
    padj_threshold: float,
    log2fc_threshold: float,
    max_labels: int,
    figsize: tuple[float, float],
    dpi: int,
) -> None:
    """Plot a volcano plot for one cluster-cell type."""
    plot_df = results.copy()

    plot_df["minus_log10_padj"] = -np.log10(
        plot_df["padj"].clip(lower=np.finfo(float).tiny),
    )

    colours = {
        "upregulated": "red",
        "downregulated": "blue",
        "nonsignificant": "grey",
    }

    fig, ax = plt.subplots(figsize=figsize)

    min_point_size = 12
    max_point_size = 120

    plot_df["relevant_fraction"] = np.where(
        plot_df["log2fc"] >= 0,
        plot_df["fraction_in_group"],
        plot_df["fraction_in_reference"],
    )

    plot_df["point_size"] = min_point_size + plot_df["relevant_fraction"].fillna(
        0,
    ).clip(0, 1) * (max_point_size - min_point_size)

    for category in [
        "nonsignificant",
        "downregulated",
        "upregulated",
    ]:
        subset = plot_df.loc[plot_df["direction"] == category]

        ax.scatter(
            subset["log2fc"],
            subset["minus_log10_padj"],
            c=colours[category],
            s=subset["point_size"],
            alpha=0.35,
            linewidths=0,
            rasterized=True,
        )

    ax.axvline(
        -log2fc_threshold,
        linestyle="--",
        linewidth=1,
        color="black",
    )

    ax.axvline(
        log2fc_threshold,
        linestyle="--",
        linewidth=1,
        color="black",
    )

    ax.axhline(
        -np.log10(padj_threshold),
        linestyle="--",
        linewidth=1,
        color="black",
    )

    labels = select_labels(
        plot_df,
        max_labels=max_labels,
    )

    # Establish final axis limits before running adjust_text.
    finite_y = plot_df["minus_log10_padj"].replace(
        [np.inf, -np.inf],
        np.nan,
    )

    maximum_y = finite_y.max()

    if pd.notna(maximum_y):
        headroom = max(
            maximum_y * 0.18,
            5.0,
        )

        ax.set_ylim(
            bottom=0,
            top=maximum_y + headroom,
        )

    # Complete the normal Matplotlib layout before adjusting labels.
    fig.tight_layout()
    fig.canvas.draw()

    text_objects = []
    target_x = []
    target_y = []

    # Use a small offset measured in data coordinates.
    initial_vertical_offset = max(
        maximum_y * 0.012,
        1.0,
    )

    for label_number, (_, row) in enumerate(labels.iterrows()):
        x = float(row["log2fc"])
        y = float(row["minus_log10_padj"])

        # Slightly stagger initial label positions.
        stagger = (label_number % 4) * initial_vertical_offset * 0.6

        text_objects.append(
            ax.text(
                x,
                y + initial_vertical_offset + stagger,
                str(row["gene"]),
                fontsize=5,
                ha="center",
                va="bottom",
                clip_on=True,
            ),
        )

        target_x.append(x)
        target_y.append(y)

    if text_objects and adjust_text is not None:
        adjust_text(
            text_objects,
            target_x=np.asarray(target_x),
            target_y=np.asarray(target_y),
            ax=ax,
            only_move={
                "text": "xy",
                "static": "xy",
                "explode": "xy",
                "pull": "xy",
            },
            # Increase vertical separation between labels.
            expand=(1.05, 2.1),
            force_text=(0.0, 0.8),
            force_explode=(0.0, 0.4),
            # Keep a weak attraction towards the corresponding point.
            force_pull=(0.0, 0.02),
            pull_threshold=10,
            # Do not repel labels from every point in the volcano.
            # The labelled points themselves are already accounted for.
            avoid_self=True,
            # Prevent the label adjustment from changing the plot limits.
            ensure_inside_axes=True,
            expand_axes=False,
            max_move=(8, 15),
            iter_lim=1000,
            prevent_crossings=True,
            min_arrow_len=4,
            arrowprops={
                "arrowstyle": "-",
                "linewidth": 0.4,
                "color": "black",
                "alpha": 0.65,
            },
        )

    elif text_objects:
        print(
            "Warning: adjustText is not installed; gene labels may overlap.",
            flush=True,
        )

    n_up = (plot_df["direction"] == "upregulated").sum()

    n_down = (plot_df["direction"] == "downregulated").sum()

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="red",
            markeredgecolor="none",
            markersize=7,
            label=f"Upregulated ({n_up:,})",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="blue",
            markeredgecolor="none",
            markersize=7,
            label=f"Downregulated ({n_down:,})",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="grey",
            markeredgecolor="none",
            markersize=7,
            label=f"Non-significant \n\
                (padj > {padj_threshold} | \n\
                |log2FC| < {log2fc_threshold} | \n\
                in-group expressed fraction < {MIN_IN_GROUP_FRACTION} | \n\
                Δ group-vs-rest expressed fraction < {MIN_FRACTION_DIFFERENCE})",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper left",
    )

    ax.set_xlabel(f"log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(
        f"{group_name}-vs-rest gene expression (in-group fraction ≥ {MIN_IN_GROUP_FRACTION}, Δ fraction ≥ {MIN_FRACTION_DIFFERENCE})",
    )

    ax.grid(False)

    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(fig)


def save_marker_tables(
    markers: pd.DataFrame,
    groupby: str,
    output_dir: Path,
) -> None:
    """Save combined and group-specific marker tables."""
    tables_dir = output_dir / "marker_tables"

    tables_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    markers.to_csv(
        output_dir / "all_cluster_markers.csv",
        index=False,
    )

    significant = markers.loc[markers["significant"]].sort_values(
        [groupby, "direction", "padj", "log2fc"],
        ascending=[True, True, True, False],
    )

    significant.to_csv(
        output_dir
        / f"significant_cluster_markers_group_frac_{MIN_IN_GROUP_FRACTION}.csv",
        index=False,
    )

    for group_name, group_results in markers.groupby(
        groupby,
        observed=True,
    ):
        filename = safe_filename(group_name)

        group_results.sort_values(
            ["padj", "log2fc"],
            ascending=[True, False],
        ).to_csv(
            tables_dir / f"{filename}_all_markers.csv",
            index=False,
        )

        (
            group_results.loc[group_results["significant"]]
            .sort_values(
                ["direction", "padj", "log2fc"],
                ascending=[True, True, False],
            )
            .to_csv(
                tables_dir
                / f"{filename}_significant_markers_group_frac_{MIN_IN_GROUP_FRACTION}.csv",
                index=False,
            )
        )


def make_volcano_plots(
    markers: pd.DataFrame,
    groupby: str,
    output_dir: Path,
    padj_threshold: float,
    log2fc_threshold: float,
    max_labels: int,
    figsize: tuple[float, float],
    dpi: int,
) -> None:
    """Generate one volcano plot per cluster-cell type."""
    volcano_dir = output_dir / "volcano_plots"

    volcano_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for group_name, group_results in markers.groupby(
        groupby,
        observed=True,
    ):
        output_path = (
            volcano_dir
            / f"{safe_filename(group_name)}_group_frac_{MIN_IN_GROUP_FRACTION}volcano.png"
        )

        print(f"Plotting volcano for {group_name}: {output_path}")

        plot_volcano(
            results=group_results,
            group_name=str(group_name),
            output_path=output_path,
            padj_threshold=padj_threshold,
            log2fc_threshold=log2fc_threshold,
            max_labels=max_labels,
            figsize=figsize,
            dpi=dpi,
        )


def main() -> None:
    """Run marker analysis and volcano plotting."""
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(f"Reading {INPUT_FILE}")
    adata = sc.read_h5ad(INPUT_FILE)

    adata = ensembl_id_to_gene_name(adata, "ensembl_dict.tsv")

    # Works for both sparse and dense matrices.
    values = adata.X.data if hasattr(adata.X, "data") else np.asarray(adata.X).ravel()

    # Check whether all values are effectively integers.
    is_integer_matrix = np.allclose(
        values,
        np.round(values),
        atol=1e-8,
    )

    if is_integer_matrix:
        print(
            "adata.X contains integer counts. Applying log1p normalization.",
        )
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print(
            "adata.X does not contain integer counts. Assuming it is already log-normalized.",
        )

    adata = create_temp_subtypes(
        adata,
        "label_spreading_prediction_filtered",
        "leiden_n30_r1p0",
        "sample_id",
        candidate_clusters,
    )

    adata = validate_anndata(
        adata=adata,
        groupby=GROUPBY,
        min_cells_per_group=MIN_CELLS_PER_GROUP,
    )

    sc.tl.rank_genes_groups(
        adata,
        groupby=GROUPBY,
        groups="all",
        reference="rest",
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        pts=True,
        use_raw=False,
        key_added="cluster_celltype_markers",
    )

    markers = extract_marker_results(
        adata=adata,
        groupby=GROUPBY,
        padj_threshold=PADJ_THRESHOLD,
        log2fc_threshold=LOG2FC_THRESHOLD,
        min_expression_fraction=MIN_IN_GROUP_FRACTION,
        min_fraction_difference=MIN_FRACTION_DIFFERENCE,
    )

    save_marker_tables(
        markers=markers,
        groupby=GROUPBY,
        output_dir=OUTPUT_DIR,
    )

    make_volcano_plots(
        markers=markers,
        groupby=GROUPBY,
        output_dir=OUTPUT_DIR,
        padj_threshold=PADJ_THRESHOLD,
        log2fc_threshold=LOG2FC_THRESHOLD,
        max_labels=MAX_LABELS,
        figsize=(7, 8.5),
        dpi=300,
    )


if __name__ == "__main__":
    main()
