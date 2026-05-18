#!/usr/bin/env python3
"""Curate cell type annotations and generate UMAP visualisations."""

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from utils.curated_UMAPs_utils import (
    compute_celltype_props,
    compute_joint_umap,
    compute_marginal_umap,
    compute_umap,
)

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

PATH = Path("project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/")

# mapping dict
map_category = {
    "CD4 T cell": "T cells",
    "CD8 T cell": "T cells",
    "Neuronal": "Neuronal",
    "Endothelial": "Endothelial",
    "Goblet": "Epithelial",
    "BEST4+ epithelial": "Epithelial",
    "CD8 T mem": "T cells",
    "cDC": "Myeloid",
    "Enterocyte": "Epithelial",
    "Pericyte": "Mesenchymal",
    "B cells": "B cells",
    "Plasma cells": "Plasma cells",
    "Enteroendocrine cell": "Epithelial",
    "ADAMDEC1+ stromal": "Mesenchymal",
    "NK T cell": "T cells",
    "NK cell": "T cells",
    "Macrophages": "Myeloid",
    "Mast cell": "Myeloid",
    "Mesenchymal": "Mesenchymal",
    "Monocytes": "Myeloid",
    "Myofibroblast": "Mesenchymal",
    "Paneth": "Epithelial",
    "Red blood cells": "Red blood cells",
    "Naïve CD8 T cell": "T cells",
    "Stem cell": "Epithelial",
    "NPY+ Stromal": "Mesenchymal",
    "Tfh": "T cells",
    "Transitional Stromal C3+": "Mesenchymal",
    "Tuft": "Epithelial",
}

reduct_name = "X_embeddings"
min_cells = 10
label = "label_spreading_prediction_filtered"


def main():
    """Plot UMAPs and cell type proporations from curated cell identities."""
    adata = sc.read_h5ad(
        PATH / "c561826c_sysvi_label_spreading_UMAP_X_embeddings_umaps.h5ad",
    )

    counts = adata.obs[label].value_counts()
    keep_labels = counts[counts >= min_cells].index
    drop_labels = counts[counts < min_cells].index

    print(
        f"labels with < {min_cells} cells: {drop_labels.tolist()}",
    )

    mask = adata.obs[label].isin(keep_labels)
    adata = adata[mask].copy()

    # count number of cells in crohn's disease and normal
    sample_id = adata.obs["sample_id"].astype(str)

    crohns_samples = sorted(
        sample_id[sample_id.str.contains("crohns", case=False, na=False)]
        .unique()
        .tolist(),
    )

    normal_samples = sorted(
        sample_id[sample_id.str.contains("normal", case=False, na=False)]
        .unique()
        .tolist(),
    )

    num_diseased_cells = adata[adata.obs["sample_id"].isin(crohns_samples)].shape[0]
    num_normal_cells = adata[~adata.obs["sample_id"].isin(crohns_samples)].shape[0]

    print(f"Number of cells in Crohn's Disease: {num_diseased_cells}")
    print(f"Number of cells in Normal: {num_normal_cells}")

    # reassign labels and save
    adata.obs["curated"] = adata.obs["label_spreading_prediction_filtered"]

    adata.obs["category"] = adata.obs["label_spreading_prediction_filtered"].map(
        map_category,
    )

    adata.obs["Diagnosis"] = adata.obs["sample_id"].isin(crohns_samples)
    adata.obs["Diagnosis"] = (
        adata.obs["Diagnosis"]
        .map({True: "Crohn's Disease", False: "Normal"})
        .astype("category")
    )

    # adata.write_h5ad(PATH / "query_concat_curated.h5ad")

    # Joint UMAP
    sc.pp.neighbors(adata, use_rep=reduct_name)
    sc.tl.umap(adata, min_dist=0.3)

    # coloured by diagnosis
    compute_umap(adata, color="Diagnosis", outdir=PATH, reduct_name=reduct_name)

    plot_configs = [
        {"color": "sample_id", "legends": True, "annotations": False},
        {"color": "category", "legends": True, "annotations": False},
        {"color": "curated", "legends": True, "annotations": False},
        {"color": "curated", "legends": False, "annotations": True},
    ]

    for config in plot_configs:
        compute_joint_umap(
            adata,
            reduct_name=reduct_name,
            outdir=PATH,
            **config,
        )

        compute_marginal_umap(
            adata,
            reduct_name=reduct_name,
            outdir=PATH,
            **config,
        )

    compute_celltype_props(adata, "category", crohns_samples, normal_samples, PATH)
    compute_celltype_props(adata, "curated", crohns_samples, normal_samples, PATH)


if __name__ == "__main__":
    main()
