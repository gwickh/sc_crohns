import gzip
import os
import re
import tempfile
import zipfile
from glob import glob

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import scipy.sparse as sparse

PROJECT_AREA = "project-area/data/crohns_scrnaseq/3c_4n_analysis"
h5ad_path = os.path.join(PROJECT_AREA, "scvi_tools_output/Integrated_05_label")


def load_count_matrices(path: str) -> list:
    """
    Load count matrices from given path into AnnData objects.
    """
    matrix_paths = sorted(
        [
            p
            for p in sorted(
                glob(
                    os.path.join(path, "**", "outs", "filtered_feature_bc_matrix.h5"),
                    recursive=True,
                )
            )
            if not os.path.isdir(p) and os.path.basename(p).endswith((".h5", ".h5ad"))
        ]
    )
    sample_names = [
        os.path.basename(os.path.dirname(os.path.dirname(p))) for p in matrix_paths
    ]

    adata_list = []
    for path, sample in zip(matrix_paths, sample_names):
        if ".h5ad" in path:
            adata = sc.read_h5ad(path)
            adata.obs["platform"] = "Parse"
        elif ".h5" in path:
            adata = sc.read_10x_h5(path)
            adata.obs["platform"] = "10X_Chromium"
            adata.var["gene_name"] = adata.var.index
        else:
            raise ValueError(f"Unknown file type: {sample}")
        adata.obs["sample_id"] = sample  # add sample_id column
        adata.var_names_make_unique()  # make var names unique
        adata.obs_names = [
            f"{sample}_{bc}" for bc in adata.obs_names
        ]  # make barcodes unique
        adata.raw = adata.copy()  # store raw counts
        adata_list.append(adata)
    return adata_list


adata_list = load_count_matrices(
    "project-area/data/crohns_scrnaseq/3c_4n_analysis/crohns_samples"
)

full_ad = ad.concat(
    adata_list,
    label="sample_id",
    keys=[a.obs["sample_id"].iloc[0] for a in adata_list],
    join="outer",
    merge="unique",
)

# load annotated query data
raw_ad = sc.read_h5ad(os.path.join(h5ad_path, "query_concat_curated.h5ad"))


# strip sample id from barcodes to match 10x format
def core10x(x: str) -> str:
    """
    Extract the core 10x barcode from a string, handling various formats.
    """
    m = re.search(r"([ACGT]{16})(?:-(\d+))?", str(x))
    if not m:
        return str(x)
    suf = m.group(2) if m.group(2) else "1"
    return f"{m.group(1)}-{suf}"


raw_core = pd.Index([core10x(x) for x in raw_ad.obs_names])
full_core = pd.Index([core10x(x) for x in full_ad.obs_names])

# map core barcode to first matching full barcode
full_lookup = {}
for c, orig in zip(full_core, full_ad.obs_names):
    full_lookup.setdefault(c, orig)

mapped = [full_lookup.get(c, None) for c in raw_core]
keep = [m is not None for m in mapped]

print("matched:", sum(keep), "/", raw_ad.n_obs)

raw_sub = raw_ad[keep].copy()
full_sub = full_ad[[m for m in mapped if m is not None], :].copy()

full_sub.obs = raw_sub.obs.copy()
full_sub.obs_names = raw_sub.obs_names

raw_ad = full_sub


def h5ad_to_10x(
    ad,
    gene_id_key="gene_ids",  # Ensembl IDs
    cell_type_key="curated",
    output_path=os.path.join(PROJECT_AREA, "marker_selection/matrix.zip"),
    barcode_key=None,
    subsample_rate=None,
):
    if subsample_rate:
        sc.pp.subsample(ad, subsample_rate)

    X = ad.layers["counts"] if "counts" in ad.layers else ad.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.astype(np.int32)

    # features.tsv.gz
    genes = pd.DataFrame(
        {
            0: ad.var[gene_id_key].astype(str).values,  # feature_id
            1: ad.var_names.astype(str).values,  # feature_name (symbols)
            2: ad.var["feature_types"].astype(str).values
            if "feature_types" in ad.var.columns
            else ["Gene Expression"] * ad.n_vars,  # feature_type
        }
    )

    # barcodes.tsv.gz
    barcodes = (
        pd.DataFrame(ad.obs[barcode_key].astype(str).values)
        if barcode_key
        else pd.DataFrame(ad.obs_names.astype(str))
    )

    # annotation CSV
    ann = pd.DataFrame(
        {
            "Barcode": ad.obs_names.astype(str),
            "Annotation": ad.obs[cell_type_key].astype(str).values,
        }
    )

    with tempfile.TemporaryDirectory() as tmp_dir:
        with gzip.open(os.path.join(tmp_dir, "matrix.mtx.gz"), "wb") as handle:
            scipy.io.mmwrite(handle, sparse.csc_matrix(X.T))

        genes.to_csv(
            os.path.join(tmp_dir, "features.tsv.gz"),
            sep="\t",
            index=False,
            header=False,
            compression="gzip",
        )
        barcodes.to_csv(
            os.path.join(tmp_dir, "barcodes.tsv.gz"),
            sep="\t",
            index=False,
            header=False,
            compression="gzip",
        )
        ann.to_csv(os.path.join(tmp_dir, "celltype_annotations.csv"), index=False)

        with zipfile.ZipFile(output_path, "w") as zip_handle:
            for fn in [
                "matrix.mtx.gz",
                "features.tsv.gz",
                "barcodes.tsv.gz",
                "celltype_annotations.csv",
            ]:
                zip_handle.write(os.path.join(tmp_dir, fn), arcname=fn)


# call
h5ad_to_10x(raw_ad, gene_id_key="gene_ids", cell_type_key="curated")

h5ad_to_10x(
    raw_ad,
    gene_id_key="gene_ids",
    gene_name_key="feature_types",
    cell_type_key="curated",
)
