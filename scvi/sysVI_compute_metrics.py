#!/usr/bin/env python3
"""Compute integration and biological preservation metrics for sysVI."""

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

# set pandas string handling to use builtin str type, not pyarrow to avoid IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

ref_dir = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/"
tuning_dir = (
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/sysvi_tuning"
)


def add_cell_type_annotation(
    cell_type_ref_adata: sc.AnnData,
    query_adata_file: Path,
    cell_type_key: str = "Integrated_05",
) -> sc.AnnData:
    """Add cell type annotation from reference data to query data."""
    query_adata = sc.read_h5ad(query_adata_file)

    if cell_type_key in query_adata.obs.columns:
        print(f"{cell_type_key} already exists in query data, skipping annotation.")
        return

    if cell_type_key not in cell_type_ref_adata.obs.columns:
        msg = f"Reference AnnData object must have '{cell_type_key}' column in .obs"
        raise ValueError(msg)

    # ensure gene names are unique
    if (
        not cell_type_ref_adata.obs_names.is_unique
        or not query_adata.obs_names.is_unique
    ):
        msg = "Observation names in either reference or query data are not unique."
        raise ValueError(msg)

    missing_ref_obs = query_adata[
        query_adata.obs["platform"] == "10X"
    ].obs_names.difference(
        cell_type_ref_adata[cell_type_ref_adata.obs["platform"] == "10X"].obs_names,
    )

    if len(missing_ref_obs) > 0:
        msg = f"""{len(missing_ref_obs)} query cells are missing from reference.\n
        query: {query_adata.obs_names[:10].tolist()}\n
        ref: {cell_type_ref_adata.obs_names[:10].tolist()}
        """
        pd.DataFrame(missing_ref_obs).to_csv(
            Path(tuning_dir) / "missing_from_reference.csv",
            index=False,
        )
        raise ValueError(msg)

    query_adata.obs[cell_type_key] = (
        cell_type_ref_adata.obs[cell_type_key].reindex(query_adata.obs_names).to_numpy()
    )

    query_adata.write_h5ad(f"{query_adata_file}_with_10X_labels.h5ad")


def main():
    """Compute metrics for sysVI."""
    cell_type_ref_adata = sc.read_h5ad(Path(ref_dir) / "query_concat.h5ad")

    adatas = list(Path(tuning_dir).glob("*_gca_sysvi.h5ad"))
    print(f"Found {len(adatas)} sysVI files.")

    for file in adatas:
        try:
            add_cell_type_annotation(
                cell_type_ref_adata=cell_type_ref_adata,
                query_adata_file=file,
            )
        except OSError as err:
            print(f"[skip] Could not read {file}: {err}")
        continue


if __name__ == "__main__":
    main()
