#!/usr/bin/env python3
"""Get markers for cell type identity and differential expression for Xenium panel."""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

#  config
PROJECT_AREA = Path("project-area/data/crohns_scrnaseq/3c_4n_analysis")
h5ad_path = (
    PROJECT_AREA / "scvi_tools_output/Integrated_05_label/query_concat_curated.h5ad"
)

OUTPATH = PROJECT_AREA / "marker_selection"
Path.mkdir(OUTPATH, exist_ok=True)

# drop tiny groups params
min_cells = 25
use_raw = False
layer = None


# compute identity markers (one-vs-rest per group)
def identity_markers(
    adata,
    to_drop,
    groupby="curated",
    method="wilcoxon",
    n_markers=50,
) -> None:
    """Compute identity markers for each cell type."""
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=n_markers,
        use_raw=False,
    )

    markers = sc.get.rank_genes_groups_df(adata, group=None)

    # compute detection fractions
    count = adata.raw.X if use_raw else adata.X
    var_names = adata.raw.var_names if use_raw else adata.var_names
    g2i = {g: i for i, g in enumerate(var_names)}
    labels = adata.obs[groupby].astype(str).to_numpy()

    def pct(mask, gene) -> float:
        j = g2i.get(gene)
        if j is None:
            return np.nan
        x = count[mask, j]
        if hasattr(x, "toarray"):
            x = x.toarray()
        return float((x > 0).mean() * 100)

    # compute pct in group, rest and difference
    markers["pct_in_group"] = [
        pct(labels == g, gene)
        for g, gene in zip(markers["group"], markers["names"], strict=False)
    ]
    markers["pct_in_rest"] = [
        pct(labels != g, gene)
        for g, gene in zip(markers["group"], markers["names"], strict=False)
    ]
    markers["pct_diff"] = markers["pct_in_group"] - markers["pct_in_rest"]
    markers["detection_ratio"] = markers["pct_in_group"] / (markers["pct_in_rest"] + 1)

    # filter low discriminatory/ abundance markers
    within_group_threshold = 60
    between_group_threshold = 20

    markers = markers[markers["pct_in_group"] >= within_group_threshold]
    markers = markers[markers["pct_in_rest"] <= between_group_threshold]

    # drop MT-/ribosomal genes
    markers = markers[~markers["names"].str.match(r"^(MT-|RPL|RPS)")]

    markers_filtered = markers[~markers["names"].astype(str).isin(set(to_drop))].copy()

    markers_filtered.to_csv(OUTPATH / "markers_identity.csv", index=False)


# compute diagnosis markers within each cell type
def deg_markers(
    adata,
    to_drop,
    groupby="diagnosis",
    method="wilcoxon",
    n_markers=200,
) -> None:
    """Compute DEG markers within each cell type."""
    sid = adata.obs["sample_id"].astype(str)
    status = np.where(
        sid.str.contains("C", case=False, na=False),
        "crohns",
        np.where(sid.str.contains("N", case=False, na=False), "normal", "other"),
    )
    adata.obs["diagnosis"] = pd.Categorical(
        status,
        categories=["normal", "crohns", "other"],
    )

    rows = []
    min_cells = 25
    min_log2fc = 2
    min_qval = 0.05

    for ct in adata.obs["curated"].astype(str).unique():
        sub = adata[adata.obs["curated"].astype(str) == ct].copy()

        # require 25 cells in each condition for comparison
        counts = sub.obs["diagnosis"].value_counts()
        if counts.get("crohns", 0) < min_cells or counts.get("normal", 0) < min_cells:
            continue

        sc.tl.rank_genes_groups(
            sub,
            groupby=groupby,
            groups=["crohns"],
            reference="normal",
            method=method,
            n_genes=n_markers,
            use_raw=use_raw,
        )

        df = sc.get.rank_genes_groups_df(sub, group="crohns")
        df.insert(0, "celltype", ct)
        rows.append(df)

    de = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    de = de[(de["logfoldchanges"] >= min_log2fc) & (de["pvals_adj"] <= min_qval)]
    de_filtered = de[~de["names"].astype(str).isin(set(to_drop))].copy()

    de_filtered.to_csv(OUTPATH / "markers_diagnosis.csv", index=False)


def main() -> None:
    """Run marker selection."""
    adata = sc.read_h5ad(h5ad_path)

    # load list of genes in Xenium core panel
    to_drop = pd.read_csv(
        OUTPATH / "XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv",
    ).iloc[:, 0]

    vc = adata.obs["curated"].value_counts()
    keep = vc[vc >= min_cells].index
    adata = adata[adata.obs["curated"].isin(keep)].copy()

    # choose matrix
    if layer is not None:
        adata.X = adata.layers[layer]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    identity_markers(adata, to_drop)
    deg_markers(adata, to_drop)


if __name__ == "__main__":
    main()
