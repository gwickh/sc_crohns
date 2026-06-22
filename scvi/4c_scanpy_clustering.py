#!/usr/bin/env python3

"""Scanpy clustering and UMAP visualization for sysVI latent space."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

tuning_dir = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/sysvi_tuning/",
)


def compute_umap(
    adata: ad.AnnData,
    use_rep: str,
    neighbors_key: str,
    umap_key: str,
    n_neighbors: int = 15,
    min_dist: float = 0.3,
) -> ad.AnnData:
    """Compute UMAP embedding."""
    sc.pp.neighbors(
        adata,
        use_rep=use_rep,
        n_neighbors=n_neighbors,
        key_added=neighbors_key,
    )
    sc.tl.umap(adata, neighbors_key=neighbors_key, min_dist=min_dist)
    adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
    return adata


def save_umap_pdf(
    adata: ad.AnnData,
    umap_key: str,
    color: str,
    out_pdf: Path,
    ncols: int,
) -> None:
    """Save UMAP plot as PDF."""
    fig = sc.pl.embedding(
        adata,
        basis=umap_key,
        color=color,
        ncols=ncols,
        frameon=True,
        return_fig=True,
        show=False,
    )
    fig.set_size_inches(6, 6)
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def run_umap_and_save_pdfs(
    adata: ad.AnnData,
    prefix: str,
    sample_set: str,
) -> ad.AnnData:
    """Run UMAP coloured by diagnosis, sample_id, platform, and predictions."""
    adata = compute_umap(
        adata,
        use_rep="X_embeddings",
        neighbors_key=f"neighbors_{sample_set}",
        umap_key=f"X_umap_{sample_set}",
    )

    for color, out_pdf in [
        (
            "diagnosis",
            tuning_dir / f"{prefix}_UMAP_{sample_set}_X_embeddings_diagnosis.pdf",
        ),
        (
            "sample_id",
            tuning_dir / f"{prefix}_UMAP_{sample_set}_X_embeddings_metadata.pdf",
        ),
        (
            "platform",
            tuning_dir / f"{prefix}_UMAP_{sample_set}_X_embeddings_platform.pdf",
        ),
        (
            "label_spreading_prediction_filtered",
            tuning_dir / f"{prefix}_UMAP_{sample_set}_X_embeddings_celltypes.pdf",
        ),
    ]:
        save_umap_pdf(
            adata,
            umap_key=f"X_umap_{sample_set}",
            color=color,
            out_pdf=out_pdf,
            ncols=1,
        )

    return adata


def main() -> None:
    """Run UMAP and save PDFs."""
    # h5ad_files = list(tuning_dir.glob("*.h5ad"))

    # prefixes = [p.name[:8] for p in h5ad_files]

    for prefix in ["c561826c_sysvi_label_spreading"]:
        adata_file = tuning_dir / "c561826c_sysvi_label_spreading_alpha_0.2_n_5.h5ad"

        adata = sc.read_h5ad(adata_file)

        # add diagnosis column based on sample_id
        adata.obs["diagnosis"] = np.where(
            adata.obs["sample_id"].str.contains("crohns", case=False, na=False),
            "Crohn's disease",
            "Normal",
        )
        # remove low confidence label transfer predictions
        adata = adata[adata.obs["label_spreading_prediction_filtered"] != "Unknown"]

        # subset to Crohn's disease and normal samples for marginal UMAPs
        sample_id = adata.obs["sample_id"].astype(str)

        crohns_samples = sorted(
            sample_id[sample_id.str.contains("crohns", case=False, na=False)]
            .unique()
            .tolist(),
        )
        adata_c = adata[adata.obs["sample_id"].isin(crohns_samples)].copy()

        # Subset to normal and compute marginal UMAPs
        normal_samples = sorted(
            sample_id[sample_id.str.contains("normal", case=False, na=False)]
            .unique()
            .tolist(),
        )
        adata_n = adata[adata.obs["sample_id"].isin(normal_samples)].copy()

        # compute UMAPs and save PDFs
        adata = run_umap_and_save_pdfs(adata, prefix, sample_set="full")
        adata_c = run_umap_and_save_pdfs(adata_c, prefix, sample_set="crohns")
        adata_n = run_umap_and_save_pdfs(adata_n, prefix, sample_set="normal")

        adata.write_h5ad(tuning_dir / f"{prefix}_UMAP_X_embeddings_umaps.h5ad")


if __name__ == "__main__":
    main()
