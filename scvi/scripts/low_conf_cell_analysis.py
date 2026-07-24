#!/usr/bin/env python3
"""Investigate distribution of low confidence cells."""

from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

adata_path = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output",
)


def main() -> None:
    """Plot UMAPs with low confidence cells highlighted."""
    adata = ad.read_h5ad(
        adata_path / "c561826c_sysvi_label_spreading_UMAP_X_embeddings_umaps.h5ad",
    )

    adata_parse = adata[adata.obs["platform"] == "Parse"]
    print(adata_parse.obs["label_spreading_confidence"].isna().sum())

    sc.pl.umap(
        adata,
        color=[
            "assignment_status",
            "label_spreading_confidence",
            "label_spreading_prediction",
            "sample_id",
        ],
        wspace=0.4,
        save="umap_all_cells.png",
    )

    sc.pl.umap(
        adata[adata.obs["assignment_status"] == "low_confidence"],
        color=[
            "assignment_status",
            "label_spreading_confidence",
            "label_spreading_prediction",
            "sample_id",
        ],
        wspace=0.4,
        save="umap_low_conf_cells.png",
    )


if __name__ == "__main__":
    main()
