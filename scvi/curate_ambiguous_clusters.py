#!/usr/bin/env python3
"""Subcluster, predict labels from CellTypist and pathway analysis by decoupler."""

from pathlib import Path

import anndata as ad
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import TwoSlopeNorm
from scipy import sparse
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

SCVI_PATH = Path("project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output")
OUTPUT_DIR = SCVI_PATH / "candidate_subcluster_hallmark"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

celltype_col = "label_spreading_prediction_filtered"

cluster_col = "leiden_n30_r1p0"
sample_col = "sample_id"
counts_layer = "counts"

candidate_clusters = {
    "B cells": ["8", "20", "21"],
    "ADAMDEC1+ stromal": ["17", "19"],
    "Macrophages": ["34", "14", "2", "3"],
    "Plasma cells": ["4", "16"],
    "Monocytes": ["26", "33"],
}


# Helper functions
def clean_name(value: str) -> str:
    """Convert a label into a safe subtype name."""
    return str(value).strip().replace(" ", "_").replace("+", "pos").replace("/", "_")


def matrix_to_dataframe(adata_obj: ad.AnnData) -> pd.DataFrame:
    """Convert an AnnData expression matrix into a sample-by-gene dataframe."""
    matrix = adata_obj.X

    if sparse.issparse(matrix):
        matrix = matrix.toarray()

    return pd.DataFrame(
        matrix,
        index=adata_obj.obs_names,
        columns=adata_obj.var_names,
    )


def load_hallmark_network(hallmark_path):
    """Load a local OmniPath/MSigDB Hallmark table for decoupler."""
    hallmark_path = Path(hallmark_path)

    if not hallmark_path.exists():
        msg = f"Hallmark resource not found: {hallmark_path}"
        raise FileNotFoundError(msg)

    hallmark = pd.read_csv(
        hallmark_path,
        sep="\t",
        compression="infer",
        dtype=str,
    )

    print("\nHallmark resource columns:")
    print(hallmark.columns.tolist())

    if {"geneset", "genesymbol"}.issubset(hallmark.columns):
        hallmark = hallmark.rename(
            columns={
                "geneset": "source",
                "genesymbol": "target",
            },
        )

    elif not {"source", "target"}.issubset(hallmark.columns):
        msg = "The Hallmark table must contain either ['geneset', 'genesymbol'] or \
            ['source', 'target'].Found: {hallmark.columns.tolist()}"
        raise ValueError(msg)

    # Retain gene-level entries only, where this metadata is available.
    if "entity_type" in hallmark.columns:
        hallmark = hallmark.loc[hallmark["entity_type"].eq("protein")].copy()

    # Retain Hallmark collection entries, where this metadata is available.
    if "collection" in hallmark.columns:
        print(
            "\nCollections found:",
            hallmark["collection"].dropna().unique().tolist(),
        )

        hallmark = hallmark.loc[
            hallmark["collection"].str.lower().str.contains("hallmark", na=False)
        ].copy()

    hallmark = (
        hallmark[["source", "target"]]
        .dropna()
        .assign(
            source=lambda frame: frame["source"].astype(str).str.strip(),
            target=lambda frame: frame["target"].astype(str).str.strip(),
        )
        .loc[lambda frame: frame["source"].ne("") & frame["target"].ne("")]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    if hallmark.empty:
        msg = "No Hallmark pathway-gene relationships after loading and filtering the \
            local table."
        raise ValueError(msg)

    print(
        f"\nLoaded {hallmark['source'].nunique():,} Hallmark pathways "
        f"and {hallmark.shape[0]:,} pathway-gene relationships.",
    )

    print("\nExample pathway-gene relationships:")
    print(hallmark.head().to_string(index=False))

    print("\nHallmark resource columns:")

    print(hallmark.columns.tolist())

    # Decoupler enrichment methods require source-target column names.
    required_columns = {"source", "target"}

    if not required_columns.issubset(hallmark.columns):
        msg = f"The local Hallmark table does not contain the expected columns \
            {sorted(required_columns)}. Found: {hallmark.columns.tolist()}"
        raise ValueError(msg)

    hallmark = hallmark[["source", "target"]].dropna().drop_duplicates()

    print(
        f"Loaded {hallmark['source'].nunique()} Hallmark pathways and "
        f"{hallmark.shape[0]:,} pathway-gene relationships.",
    )

    return hallmark


def signed_stouffer(group, score_col="score", pvalue_col="pvalue"):
    """Compute signed Stouffer statistic."""
    pvalues = group[pvalue_col].to_numpy(dtype=float)
    scores = group[score_col].to_numpy(dtype=float)
    # Convert two-sided p-values into signed z-statistics.

    valid = np.isfinite(scores) & np.isfinite(pvalues) & (pvalues > 0) & (pvalues <= 1)
    scores = scores[valid]
    pvalues = pvalues[valid]

    if len(pvalues) == 0:
        return pd.Series(
            {
                "stouffer_z": np.nan,
                "stouffer_pvalue": np.nan,
                "n_stouffer": 0,
            },
        )

    # Avoid infinite z-scores when p-values underflow to exactly zero.

    pvalues = np.clip(
        pvalues,
        np.finfo(float).tiny,
        1.0,
    )

    # Assumes the individual p-values are two-sided.
    signed_z = np.sign(scores) * norm.isf(pvalues / 2)
    combined_z = signed_z.sum() / np.sqrt(len(signed_z))
    combined_pvalue = 2 * norm.sf(abs(combined_z))

    return pd.Series(
        {
            "stouffer_z": combined_z,
            "stouffer_pvalue": combined_pvalue,
            "n_stouffer": len(signed_z),
        },
    )


def ensembl_id_to_gene_name(adata, dict_map):
    """Check var_names are gene names."""
    ensembl_fraction = np.mean(adata.var_names.astype(str).str.startswith("ENSG"))

    if ensembl_fraction > 0.5:
        print(
            f"{ensembl_fraction * 100}% of adata.var_names appear to be Ensembl IDs. \
            Hallmark gene sets use gene symbols",
        )

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


def mask_non_candidate_cells(
    adata,
    celltype_obs,
    cluster_obs,
):
    """Subset to only cells in celltype-clusters of interest."""
    subset_mask = adata.obs["temp_subtypes"].astype(str).ne("Not_candidate")

    adata_subset = adata[subset_mask].copy()

    if adata_subset.n_obs == 0:
        msg = f"No candidate cells were found. Check the exact values \
            in {celltype_obs!r} and {cluster_obs!r}."
        raise ValueError(msg)

    print(
        f"\nRetained {adata_subset.n_obs:,} candidate cells "
        f"from {adata_subset.obs['temp_subtypes'].nunique()} subclusters.",
    )

    # Pseudobulking must use raw counts.
    adata_subset.X = adata_subset.layers["counts"].copy()

    # Remove genes with no counts anywhere in the candidate dataset.
    sc.pp.filter_genes(adata_subset, min_counts=1)

    return adata_subset


def pseudobulk_samples(
    adata_subset,
    sample_id,
    min_cell_per_pseudobulk=10,
    min_counts_per_pseudobulk=1000,
):
    """Pseudobulk cells by sample and temporary subtype, then normalise."""
    adata_subset = adata_subset.copy()

    required_obs = [sample_id, "temp_subtypes"]

    missing_obs = [
        column for column in required_obs if column not in adata_subset.obs.columns
    ]
    if missing_obs:
        msg = f"Missing required adata.obs columns: {missing_obs}"
        raise KeyError(msg)

    # decoupler constructs boolean masks internally. Pandas nullable string
    # columns produce BooleanArray masks, which cannot index SciPy sparse
    # matrices. Convert the grouping columns to ordinary object dtype.

    for column in required_obs:
        n_missing = int(adata_subset.obs[column].isna().sum())

        if n_missing > 0:
            msg = f"adata.obs[{column!r}] contains {n_missing} missing values."
            raise ValueError(msg)

        adata_subset.obs[column] = adata_subset.obs[column].astype(str).astype(object)

    print("\nGrouping-column dtypes:")

    print(adata_subset.obs[required_obs].dtypes)
    pdata = dc.pp.pseudobulk(
        adata=adata_subset,
        sample_col=sample_id,
        groups_col="temp_subtypes",
        layer=None,
        mode="sum",
    )
    print("\nPseudobulk object before filtering:")
    print(pdata)

    # Record all profiles before filtering.
    pseudobulk_qc = pdata.obs.copy()
    pseudobulk_qc.to_csv(OUTPUT_DIR / "pseudobulk_qc_before_filtering.csv")

    dc.pp.filter_samples(
        pdata,
        min_cells=min_cell_per_pseudobulk,
        min_counts=min_counts_per_pseudobulk,
    )

    if pdata.n_obs == 0:
        msg = "No pseudobulk profiles passed filtering. Reduce \
            min_cell_per_pseudobulk or min_counts_per_pseudobulk \
            after inspecting the QC table."
        raise ValueError(msg)

    pdata.obs.to_csv(OUTPUT_DIR / "pseudobulk_qc_after_filtering.csv")

    # Retain raw summed counts.
    pdata.layers["counts"] = pdata.X.copy()

    sc.pp.normalize_total(
        pdata,
        target_sum=1e6,
    )

    sc.pp.log1p(pdata)

    pdata.layers["log1p_cpm"] = pdata.X.copy()

    return pdata


def run_hallmark_enrichment(
    pdata,
    hallmark_db,
    sample_id,
    n_pathways=1000,
):
    """Run Hallmark enrichment analysis on a pseudobulk object."""
    hallmark = load_hallmark_network(OUTPUT_DIR / hallmark_db)
    hallmark = hallmark.drop_duplicates(subset=["source", "target"])

    genes_in_data = set(pdata.var_names)
    genes_in_hallmark = set(hallmark["target"])
    overlap = genes_in_data.intersection(genes_in_hallmark)

    print(
        f"\nGene overlap with Hallmark: "
        f"{len(overlap):,}/{pdata.n_vars:,} genes in the pseudobulk object."
    )

    if len(overlap) < 100:
        raise ValueError("Very few genes overlap the Hallmark resource.")

    expression = matrix_to_dataframe(pdata)

    # The second returned matrix contains raw ULM p-values.
    hallmark_scores, hallmark_pvalues = dc.mt.ulm(
        data=expression,
        net=hallmark,
    )

    hallmark_scores.index = pdata.obs_names
    hallmark_pvalues.index = pdata.obs_names

    hallmark_scores.to_csv(OUTPUT_DIR / "hallmark_scores_per_pseudobulk.csv")
    hallmark_pvalues.to_csv(OUTPUT_DIR / "hallmark_pvalues_per_pseudobulk.csv")

    score_long = (
        hallmark_scores.rename_axis(
            index="pseudobulk_id",
            columns="hallmark",
        )
        .reset_index()
        .melt(
            id_vars="pseudobulk_id",
            var_name="hallmark",
            value_name="score",
        )
    )

    pvalue_long = (
        hallmark_pvalues.rename_axis(
            index="pseudobulk_id",
            columns="hallmark",
        )
        .reset_index()
        .melt(
            id_vars="pseudobulk_id",
            var_name="hallmark",
            value_name="pvalue",
        )
    )

    results_long = score_long.merge(
        pvalue_long,
        on=["pseudobulk_id", "hallmark"],
        how="left",
        validate="one_to_one",
    )

    metadata = (
        pdata.obs[
            [
                sample_id,
                "temp_subtypes",
                "psbulk_cells",
                "psbulk_counts",
            ]
        ]
        .copy()
        .rename_axis("pseudobulk_id")
        .reset_index()
    )

    results_long = results_long.merge(
        metadata,
        on="pseudobulk_id",
        how="left",
        validate="many_to_one",
    )

    # BH correction across Hallmark pathways separately within each
    # pseudobulk expression profile.
    results_long["padj"] = results_long.groupby(
        "pseudobulk_id",
        observed=True,
    )["pvalue"].transform(
        lambda p: multipletests(
            p.to_numpy(),
            method="fdr_bh",
        )[1]
    )

    results_long.to_csv(
        OUTPUT_DIR / "hallmark_results_long.csv",
        index=False,
    )

    pathway_summary = results_long.groupby(
        ["temp_subtypes", "hallmark"],
        observed=True,
        as_index=False,
    ).agg(
        mean_score=("score", "mean"),
        median_score=("score", "median"),
        score_sd=("score", "std"),
        fraction_positive=(
            "score",
            lambda x: np.mean(x > 0),
        ),
        fraction_sig=(
            "padj",
            lambda x: np.mean(x < 0.05),
        ),
        n_samples=(sample_id, "nunique"),
    )

    stouffer_summary = (
        results_long.groupby(
            ["temp_subtypes", "hallmark"],
            observed=True,
        )
        .apply(
            signed_stouffer,
            score_col="score",
            pvalue_col="pvalue",
            include_groups=False,
        )
        .reset_index()
    )

    pathway_summary = pathway_summary.merge(
        stouffer_summary,
        on=["temp_subtypes", "hallmark"],
        how="left",
        validate="one_to_one",
    )

    pathway_summary["stouffer_padj"] = np.nan

    # Correct the combined pathway p-values across Hallmark pathways,
    # independently within each subtype.
    for _, index in pathway_summary.groupby(
        "temp_subtypes",
        observed=True,
    ).groups.items():
        index = np.asarray(list(index))

        valid_index = index[
            pathway_summary.loc[
                index,
                "stouffer_pvalue",
            ]
            .notna()
            .to_numpy()
        ]

        if valid_index.size == 0:
            continue

        pathway_summary.loc[
            valid_index,
            "stouffer_padj",
        ] = multipletests(
            pathway_summary.loc[
                valid_index,
                "stouffer_pvalue",
            ],
            method="fdr_bh",
        )[1]

    pathway_summary = pathway_summary.loc[
        pathway_summary["stouffer_padj"] < 0.05
    ].copy()

    # Rank positive pathways by median rather than mean activity.
    pathway_summary["rank"] = (
        pathway_summary.groupby(
            "temp_subtypes",
            observed=True,
        )["median_score"]
        .rank(
            method="first",
            ascending=False,
        )
        .astype(int)
    )

    pathway_summary = pathway_summary.sort_values(["temp_subtypes", "rank"])

    pathway_summary.to_csv(
        OUTPUT_DIR / "hallmark_pathway_summary_all.csv",
        index=False,
    )

    top_pathways = pathway_summary.loc[pathway_summary["rank"] <= n_pathways].copy()

    top_pathways.to_csv(
        OUTPUT_DIR / f"top_{n_pathways}_hallmark_pathways_per_subtype.csv",
        index=False,
    )

    return top_pathways


def plot_score_heatmap(
    top_pathways,
    subtype_col="temp_subtypes",
    pathway_col="hallmark",
    score_col="median_score",
    padj_col="stouffer_padj",
    top_n=None,
    figsize=None,
    min_dot_size=20,
    max_dot_size=300,
    significance_cap=10,
    cmap="coolwarm",
    outfile="top_hallmark_pathways_per_subtype.pdf",
):
    """Plot heatmap with significant enriched pathways."""
    plot_df = top_pathways[
        [
            subtype_col,
            pathway_col,
            score_col,
            padj_col,
        ]
    ].copy()

    plot_df = plot_df.dropna(
        subset=[
            subtype_col,
            pathway_col,
            score_col,
            padj_col,
        ]
    )

    # select the strongest significant pathways per subtype.
    if top_n is not None:
        plot_df["absolute_score"] = plot_df[score_col].abs()

        plot_df = (
            plot_df.sort_values(
                [subtype_col, "absolute_score", padj_col],
                ascending=[True, False, True],
            )
            .groupby(subtype_col, observed=True)
            .head(top_n)
            .drop(columns="absolute_score")
        )

    # Keep every selected pathway across all displayed subtypes.

    pathways = (
        plot_df.groupby(pathway_col, observed=True)[score_col]
        .apply(lambda x: x.abs().max())
        .sort_values(ascending=True)
        .index.tolist()
    )

    subtypes = plot_df[subtype_col].drop_duplicates().astype(str).tolist()

    pathway_positions = {pathway: position for position, pathway in enumerate(pathways)}

    subtype_positions = {subtype: position for position, subtype in enumerate(subtypes)}

    plot_df["_x"] = plot_df[subtype_col].astype(str).map(subtype_positions)

    plot_df["_y"] = plot_df[pathway_col].map(pathway_positions)

    # Avoid log10(0), then cap extreme values for plotting.

    minimum_positive = np.nextafter(0, 1)

    plot_df["_neg_log10_padj"] = -np.log10(
        plot_df[padj_col].clip(lower=minimum_positive)
    )

    plot_df["_size_value"] = plot_df["_neg_log10_padj"].clip(upper=significance_cap)

    size_min = 0

    size_max = 10

    if np.isclose(size_min, size_max):
        plot_df["_dot_size"] = (min_dot_size + max_dot_size) / 2

    else:
        plot_df["_dot_size"] = min_dot_size + (plot_df["_size_value"] - size_min) / (
            size_max - size_min
        ) * (max_dot_size - min_dot_size)

    max_absolute_score = plot_df[score_col].abs().max()

    if max_absolute_score == 0:
        max_absolute_score = 1

    colour_norm = TwoSlopeNorm(
        vmin=0,
        vcenter=7.5,
        vmax=15,
    )

    if figsize is None:
        figsize = (
            max(7, 0.8 * len(subtypes) + 3),
            max(5, 0.35 * len(pathways) + 2),
        )

    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        plot_df["_x"],
        plot_df["_y"],
        c=plot_df[score_col],
        s=plot_df["_dot_size"],
        cmap=cmap,
        norm=colour_norm,
        edgecolors="black",
        linewidths=0.4,
    )

    ax.set_xticks(range(len(subtypes)))

    ax.set_xticklabels(
        subtypes,
        rotation=45,
        ha="right",
    )

    ax.set_yticks(range(len(pathways)))

    ax.set_yticklabels(pathways)

    ax.set_xlim(-0.5, len(subtypes) - 0.5)

    ax.set_ylim(-0.5, len(pathways) - 0.5)

    ax.set_xlabel("Cluster")

    ax.set_ylabel("Hallmark pathway")

    ax.set_title("Significantly enriched Hallmark pathways")

    ax.grid(
        visible=True,
        axis="both",
        linewidth=0.5,
        alpha=0.3,
    )

    ax.set_axisbelow(True)

    colour_bar = fig.colorbar(
        scatter,
        ax=ax,
        pad=0.02,
    )

    colour_bar.set_label("Median ULM score")

    legend_values = np.linspace(
        size_min,
        size_max,
        num=min(4, len(plot_df)),
    )

    legend_values = np.unique(np.round(legend_values, decimals=1))

    legend_handles = []

    for value in legend_values:
        if np.isclose(size_min, size_max):
            marker_size = (min_dot_size + max_dot_size) / 2

        else:
            marker_size = min_dot_size + (value - size_min) / (size_max - size_min) * (
                max_dot_size - min_dot_size
            )

        legend_handles.append(
            ax.scatter(
                [],
                [],
                s=marker_size,
                facecolors="none",
                edgecolors="black",
                label=f"{value:g}",
            )
        )

    ax.legend(
        handles=legend_handles,
        title="−log10 adjusted p",
        bbox_to_anchor=(1.18, 1),
        loc="upper left",
        frameon=False,
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / outfile,
        dpi=300,
        bbox_inches="tight",
    )


def main():
    """Run subclustering, predict labels and pathway analysis by decoupler."""
    adata = sc.read_h5ad(
        SCVI_PATH / "query_concat_curated_clustered.h5ad",
    )

    adata = ensembl_id_to_gene_name(adata, "ensembl_dict.tsv")
    adata = create_temp_subtypes(
        adata,
        celltype_col,
        cluster_col,
        sample_col,
        candidate_clusters,
    )
    adata_subset = mask_non_candidate_cells(adata, celltype_col, cluster_col)

    pdata = pseudobulk_samples(adata_subset, sample_col)

    top_pathways = run_hallmark_enrichment(pdata, "msigdb-hallmark.tsv.gz", sample_col)

    plot_score_heatmap(
        top_pathways,
        subtype_col="temp_subtypes",
        pathway_col="hallmark",
        score_col="median_score",
        padj_col="stouffer_padj",
        top_n=20,
        figsize=(10, 16),
        min_dot_size=20,
        max_dot_size=300,
        significance_cap=10,
        cmap="coolwarm",
        outfile="top_hallmark_pathways_per_subtype.pdf",
    )

    pdata.write_h5ad(OUTPUT_DIR / "candidate_subclusters_pseudobulk.h5ad")


if __name__ == "__main__":
    main()
