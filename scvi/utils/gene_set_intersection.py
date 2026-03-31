#!/usr/bin/env python3
"""Compute intersection of gene sets between 10X Chromium and Parse samples."""

import numpy as np
import scanpy as sc

ADATA_OBJ_PATH = (
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/,"
    "scanpy/adata_merged_clustered.h5ad"
)


def main() -> None:
    """Compute gene set intersection."""
    adata = sc.read_h5ad(ADATA_OBJ_PATH)

    adata_chromium = adata[adata.obs["platform"] == "10X_Chromium"]
    adata_chromium = adata_chromium[
        :,
        np.asarray(adata_chromium.X.sum(axis=0)).ravel() > 0,
    ].copy()
    chromium_gene_set = set(adata_chromium.var["gene_id"])

    adata_parse = adata[adata.obs["platform"] == "Parse"]

    parse_gene_set = set(adata_parse.var["gene_id"])

    print(f"Parse gene set size: {len(parse_gene_set)}")
    print(f"Chromium gene set size: {len(chromium_gene_set)}")
    print(
        (
            "Intersection gene set size:",
            f"{len(parse_gene_set.intersection(chromium_gene_set))}",
        ),
    )


if __name__ == "__main__":
    main()
