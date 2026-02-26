import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# ---- config ----
PROJECT_AREA = "project-area/data/crohns_scrnaseq/3c_4n_analysis"
h5ad_path = os.path.join(
    PROJECT_AREA, "scvi_tools_output/Integrated_05_label/query_concat_curated.h5ad"
)
groupby = "curated"
method = "wilcoxon"
n_markers = 50
min_cells = 25
use_raw = False
layer = None
OUTPATH = os.path.join(PROJECT_AREA, "marker_selection")
os.makedirs(OUTPATH, exist_ok=True)

adata = sc.read_h5ad(h5ad_path)

# drop tiny groups
vc = adata.obs[groupby].value_counts()
keep = vc[vc >= min_cells].index
adata = adata[adata.obs[groupby].isin(keep)].copy()

# choose matrix
if layer is not None:
    adata.X = adata.layers[layer]
    use_raw = False

# compute identity markers (one-vs-rest per group)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.tl.rank_genes_groups(
    adata,
    groupby=groupby,
    method=method,
    n_genes=n_markers,
    use_raw=use_raw,
)

markers = sc.get.rank_genes_groups_df(adata, group=None)

# compute detection fractions
X = adata.raw.X if use_raw else adata.X
var_names = adata.raw.var_names if use_raw else adata.var_names
g2i = {g: i for i, g in enumerate(var_names)}
labels = adata.obs[groupby].astype(str).values


def pct(mask, gene):
    j = g2i.get(gene)
    if j is None:
        return np.nan
    x = X[mask, j]
    if hasattr(x, "toarray"):
        x = x.toarray()
    return float((x > 0).mean() * 100)


# compute pct in group, rest and difference
markers["pct_in_group"] = [
    pct(labels == g, gene) for g, gene in zip(markers["group"], markers["names"])
]

markers["pct_in_rest"] = [
    pct(labels != g, gene) for g, gene in zip(markers["group"], markers["names"])
]

markers["pct_diff"] = markers["pct_in_group"] - markers["pct_in_rest"]

markers["detection_ratio"] = markers["pct_in_group"] / (markers["pct_in_rest"] + 1)

# filter low discriminatory/ abundance markers
markers = markers[markers["pct_in_group"] >= 60]
markers = markers[markers["pct_in_rest"] <= 20]

# drop MT-/ribosomal genes
markers = markers[~markers["names"].str.match(r"^(MT-|RPL|RPS)")]

to_drop = pd.read_csv(
    os.path.join(OUTPATH, "XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")
).columns[0]
markers_filtered = markers[~markers["names"].astype(str).isin(set(to_drop))].copy()

markers_filtered.to_csv(os.path.join(OUTPATH, "markers_identity.csv"), index=False)


# compute diagnosis markers within each cell type
sid = adata.obs["sample_id"].astype(str)
status = np.where(
    sid.str.contains("C", case=False, na=False),
    "crohns",
    np.where(sid.str.contains("N", case=False, na=False), "normal", "other"),
)
adata.obs["diagnosis"] = pd.Categorical(
    status, categories=["normal", "crohns", "other"]
)

rows = []
for ct in adata.obs["curated"].astype(str).unique():
    sub = adata[adata.obs["curated"].astype(str) == ct].copy()

    # require 25 cells in each condition for comparison
    counts = sub.obs["diagnosis"].value_counts()
    if counts.get("crohns", 0) < 25 or counts.get("normal", 0) < 25:
        continue

    sc.tl.rank_genes_groups(
        sub,
        groupby="diagnosis",
        groups=["crohns"],
        reference="normal",
        method="wilcoxon",
        n_genes=200,
        use_raw=use_raw,
    )

    df = sc.get.rank_genes_groups_df(sub, group="crohns")
    df.insert(0, "celltype", ct)
    rows.append(df)

de = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
de = de[(de["logfoldchanges"] >= 2) & (de["pvals_adj"] <= 0.05)]
de_filtered = de[~de["names"].astype(str).isin(set(to_drop))].copy()

de_filtered.to_csv(os.path.join(OUTPATH, "markers_diagnosis.csv"), index=False)
