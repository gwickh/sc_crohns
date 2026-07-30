#!/usr/bin/env python3
"""Investigate distribution of low confidence cells."""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

adata_path = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/",
)


def main() -> None:
    """Plot full-object UMAP and low-confidence-only UMAP."""
    adata = ad.read_h5ad(
        adata_path
        / "sysvi_tuning"
        / "c561826c_sysvi_label_spreading_alpha_0.2_n_5.h5ad",
    )

    adata.obs["assignment_status"] = np.where(
        adata.obs["label_spreading_confidence"] < 0.8,
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

    # UMAP using only low-confidence cells
    low_conf = adata[adata.obs["assignment_status"].eq("low_conf")].copy()

    sc.pp.neighbors(
        low_conf,
        use_rep="X_embeddings",
        n_neighbors=15,
    )

    sc.tl.umap(
        low_conf,
        min_dist=0.3,
    )

    fig = sc.pl.embedding(
        low_conf,
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
        adata_path / "umap_low_conf_subset_0p8.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
