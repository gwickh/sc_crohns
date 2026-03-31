#!/usr/bin/env python3
"""Convert .h5ad files to 10x .mtx format."""

import gzip
import re
import tempfile
import zipfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy import sparse

PROJECT_AREA = Path("project-area/data/crohns_scrnaseq/3c_4n_analysis")
h5ad_path = PROJECT_AREA / "scvi_tools_output/Integrated_05_label"


# load count matrices from given path into AnnData objects
def load_count_matrices(path: Path) -> list:
    """Load count matrices from given path into AnnData objects."""
    matrix_paths = sorted(
        [
            p
            for p in path.glob("** / outs / filtered_feature_bc_matrix.h5")
            if p.is_file() and p.suffix == ".h5"
        ],
    )
    sample_names = [Path(p).parent.parent.name for p in matrix_paths]

    adata_list = []
    for matrix_path, sample in zip(matrix_paths, sample_names, strict=False):
        if ".h5" in matrix_path:
            adata = sc.read_10x_h5(matrix_path)
            adata.obs["platform"] = "10X_Chromium"
            adata.var["gene_name"] = adata.var.index
        else:
            error_msg = f"Expected .h5 file for sample '{sample}', but got '{path}'."
            raise ValueError(error_msg)
        adata.obs["sample_id"] = sample  # add sample_id column
        adata.var_names_make_unique()  # make var names unique
        adata.obs_names = [
            f"{sample}_{bc}" for bc in adata.obs_names
        ]  # make barcodes unique
        adata.raw = adata.copy()  # store raw counts
        adata_list.append(adata)
    return adata_list


def match_barcodes(raw_ad, full_ad) -> ad.AnnData:
    """Check matching barcodes between raw and full data."""
    raw_core = pd.Index([core10x(x) for x in raw_ad.obs_names])
    full_core = pd.Index([core10x(x) for x in full_ad.obs_names])

    # map core barcode to first matching full barcode
    full_lookup = {}
    for c, orig in zip(full_core, full_ad.obs_names, strict=False):
        full_lookup.setdefault(c, orig)

    mapped = [full_lookup.get(c) for c in raw_core]
    keep = [m is not None for m in mapped]

    print("matched:", sum(keep), "/", raw_ad.n_obs)

    raw_sub = raw_ad[keep].copy()
    full_sub = full_ad[[m for m in mapped if m is not None], :].copy()

    full_sub.obs = raw_sub.obs.copy()
    full_sub.obs_names = raw_sub.obs_names

    return full_sub


# strip sample id from barcodes to match 10x format
def core10x(x: str) -> str:
    """Extract the core 10x barcode from a string, handling various formats."""
    m = re.search(r"([ACGT]{16})(?:-(\d+))?", str(x))
    if not m:
        return str(x)
    suf = m.group(2) if m.group(2) else "1"
    return f"{m.group(1)}-{suf}"


# convert h5ad to 10x format and save as mtx in zip file with features and barcodes
def h5ad_to_10x(
    ad,
    gene_id_key="gene_ids",
    cell_type_key="curated",
    barcode_key=None,
    subsample_rate=None,
):
    """Save as mtx in zip file with features and barcodes files."""
    if subsample_rate:
        sc.pp.subsample(ad, subsample_rate)

    count = ad.layers.get("counts", ad.count)
    if not sparse.issparse(count):
        count = sparse.csr_matrix(count)
    count = count.astype(np.int32)

    # features.tsv.gz
    genes = pd.DataFrame(
        {
            0: ad.var[gene_id_key].astype(str).to_numpy(),  # feature_id
            1: ad.var_names.astype(str).to_numpy(),  # feature_name (symbols)
            2: ad.var["feature_types"].astype(str).to_numpy()
            if "feature_types" in ad.var.columns
            else ["Gene Expression"] * ad.n_vars,  # feature_type
        },
    )

    # barcodes.tsv.gz
    barcodes = (
        pd.DataFrame(ad.obs[barcode_key].astype(str).to_numpy())
        if barcode_key
        else pd.DataFrame(ad.obs_names.astype(str).to_numpy())
    )

    # annotation CSV
    ann = pd.DataFrame(
        {
            "Barcode": ad.obs_names.astype(str).to_numpy(),
            "Annotation": ad.obs[cell_type_key].astype(str).to_numpy(),
        },
    )

    with tempfile.TemporaryDirectory() as tmp_dir:
        with gzip.open(Path.join(tmp_dir, "matrix.mtx.gz"), "wb") as handle:
            scipy.io.mmwrite(handle, sparse.csc_matrix(count.T))

        genes.to_csv(
            Path.join(tmp_dir, "features.tsv.gz"),
            sep="\t",
            index=False,
            header=False,
            compression="gzip",
        )
        barcodes.to_csv(
            Path.join(tmp_dir, "barcodes.tsv.gz"),
            sep="\t",
            index=False,
            header=False,
            compression="gzip",
        )
        ann.to_csv(Path.join(tmp_dir, "celltype_annotations.csv"), index=False)

        output_path = Path.join(PROJECT_AREA, "marker_selection/matrix.zip")
        with zipfile.ZipFile(output_path, "w") as zip_handle:
            for fn in [
                "matrix.mtx.gz",
                "features.tsv.gz",
                "barcodes.tsv.gz",
                "celltype_annotations.csv",
            ]:
                zip_handle.write(Path.join(tmp_dir, fn), arcname=fn)


def main() -> None:
    """Perform conversion."""
    adata_list = load_count_matrices(
        "project-area/data/crohns_scrnaseq/3c_4n_analysis/crohns_samples",
    )

    full_ad = ad.concat(
        adata_list,
        keys=[a.obs["sample_id"].iloc[0] for a in adata_list],
        join="outer",
        merge="unique",
    )

    # load annotated query data
    raw_ad = sc.read_h5ad(Path.join(h5ad_path, "query_concat_curated.h5ad"))

    # match barcodes between raw and full data, keeping only those that match
    raw_ad = match_barcodes(raw_ad, full_ad)

    # convert to 10x format and save as mtx in zip file with features and barcodes files
    h5ad_to_10x(
        raw_ad,
        gene_id_key="gene_ids",
        cell_type_key="curated",
    )


if __name__ == "__main__":
    main()
