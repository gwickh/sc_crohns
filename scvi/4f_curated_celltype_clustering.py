#!/usr/bin/env python3
"""Curate cell type annotations and generate UMAP visualisations."""

from pathlib import Path

import anndata as ad
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

PATH = Path("project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/")
adata = sc.read_h5ad(
    PATH / "c561826c_sysvi_label_spreading_UMAP_X_embeddings_umaps.h5ad",
)

# filter rare labels
REDUCT_NAME = "X_embeddings"

min_cells = 10
label = "label_spreading_prediction_filtered"
counts = adata.obs[label].value_counts()
keep_labels = counts[counts >= min_cells].index
mask = adata.obs[label].isin(keep_labels)

adata = adata[mask].copy()

# count number of cells in crohn's disease and normal
sample_id = adata.obs["sample_id"].astype(str)

crohns_samples = sorted(
    sample_id[sample_id.str.contains("crohns", case=False, na=False)].unique().tolist(),
)

normal_samples = sorted(
    sample_id[sample_id.str.contains("normal", case=False, na=False)].unique().tolist(),
)

num_diseased_cells = adata[adata.obs["sample_id"].isin(crohns_samples)].shape[0]
num_normal_cells = adata[~adata.obs["sample_id"].isin(crohns_samples)].shape[0]

print(f"Number of cells in Crohn's Disease: {num_diseased_cells}")
print(f"Number of cells in Normal: {num_normal_cells}")

# mapping dicts
print("mapping dicts")
map_category = {
    "Activated CD4 T": "T cells",
    "Activated CD8 T": "T cells",
    "Adult Glia": "Neuronal",
    "arterial capillary": "Endothelial",
    "BEST2+ Goblet cell": "Epithelial",
    "BEST4+ epithelial": "Epithelial",
    "CD8 Tmem": "T cells",
    "cDC1": "Myeloid",
    "cDC2": "Myeloid",
    "Colonocyte": "Epithelial",
    "Contractile pericyte (PLN+)": "Mesenchymal",
    "CX3CR1+ CD8 Tmem": "T cells",
    "Cycling B cell": "B cells",
    "Cycling plasma cell": "Plasma cells",
    "D cells (SST+)": "Epithelial",
    "DZ GC cell": "B cells",
    "EC cells (TAC1+)": "Epithelial",
    "EECs": "Epithelial",
    "Enterocyte": "Epithelial",
    "FCRL4+ Memory B": "B cells",
    "FDC": "Mesenchymal",
    "GC B cell": "B cells",
    "gdT": "T cells",
    "Goblet cell": "Epithelial",
    "IgA plasma cell": "Plasma cells",
    "IgG plasma cell": "Plasma cells",
    "IgM plasma cell": "Plasma cells",
    "ILC3": "T cells",
    "Immature B": "B cells",
    "L cells (PYY+)": "Epithelial",
    "LEC1 (ACKR4+)": "Endothelial",
    "LEC3 (ADGRG3+)": "Endothelial",
    "LEC5 (CLDN11+)": "Endothelial",
    "LEC6 (ADAMTS4+)": "Endothelial",
    "Lymphoid DC": "Myeloid",
    "LYVE1+ Macrophage": "Myeloid",
    "LZ GC cell": "B cells",
    "Macrophages": "Myeloid",
    "MAIT cell": "T cells",
    "Mast cell": "Myeloid",
    "Mature arterial EC": "Endothelial",
    "Mature venous EC": "Endothelial",
    "Memory B": "B cells",
    "Mesothelium (PRG4+)": "Mesenchymal",
    "Microfold cell": "Epithelial",
    "mLN Stroma (FMO2+)": "Mesenchymal",
    "mLTo": "Mesenchymal",
    "MMP9+ Inflammatory macrophage": "Myeloid",
    "Monocytes": "Myeloid",
    "myofibroblast": "Mesenchymal",
    "N cells (NTS+)": "Epithelial",
    "Naive B": "B cells",
    "NK cell": "T cells",
    "NK T cell": "T cells",
    "Paneth": "Epithelial",
    "pDC": "Myeloid",
    "Pericyte": "Mesenchymal",
    "RBC": "Red blood cells",
    "SELL+ CD4 T": "T cells",
    "SELL+ CD8 T": "T cells",
    "STAT1+ Naive B": "B cells",
    "Stem cells": "Epithelial",
    "Stromal 1 (ADAMDEC1+)": "Mesenchymal",
    "Stromal 1 (CCL11+)": "Mesenchymal",
    "Stromal 2 (NPY+)": "Mesenchymal",
    "Stromal 3 (C7+)": "Mesenchymal",
    "Stromal 4 (MMP1+)": "Mesenchymal",
    "T reticular": "Mesenchymal",
    "TA": "Epithelial",
    "Tfh": "T cells",
    "Th1": "T cells",
    "Th17": "T cells",
    "Transitional Stromal 3 (C3+)": "Mesenchymal",
    "Treg": "T cells",
    "TRGV2 gdT": "T cells",
    "Tuft": "Epithelial",
    "Fetal venous EC": "Endothelial",
    "TRGV4 gdT": "T cells",
    "TRGV5/7 gdT": "T cells",
    "myofibroblast (RSPO2+)": "Mesenchymal",
    "LEC4 (STAB2+)": "Endothelial",
    "LEC2 (MADCAM1+)": "Endothelial",
    "Progenitor (NEUROG3+)": "Epithelial",
}
map_curated = {
    "Activated CD4 T": "CD4 T cell",
    "Activated CD8 T": "CD8 T cell",
    "Adult Glia": "Neuronal",
    "arterial capillary": "Endothelial",
    "BEST2+ Goblet cell": "Goblet",
    "BEST4+ epithelial": "BEST4+ epithelial",
    "CD8 Tmem": "CD8 T mem",
    "cDC1": "cDC",
    "cDC2": "cDC",
    "Colonocyte": "Enterocyte",
    "Contractile pericyte (PLN+)": "Pericyte",
    "CX3CR1+ CD8 Tmem": "CD8 T mem",
    "Cycling B cell": "B cells",
    "Cycling plasma cell": "Plasma cells",
    "D cells (SST+)": "Enteroendocrine cell",
    "DZ GC cell": "B cells",
    "EC cells (TAC1+)": "Enteroendocrine cell",
    "EECs": "Enteroendocrine cell",
    "Enterocyte": "Enterocyte",
    "FCRL4+ Memory B": "B cells",
    "FDC": "ADAMDEC1+ stromal",
    "GC B cell": "B cells",
    "gdT": "NK T cell",
    "Goblet cell": "Goblet",
    "IgA plasma cell": "Plasma cells",
    "IgG plasma cell": "Plasma cells",
    "IgM plasma cell": "Plasma cells",
    "ILC3": "NK cell",
    "Immature B": "B cells",
    "L cells (PYY+)": "Enteroendocrine cell",
    "LEC1 (ACKR4+)": "Endothelial",
    "LEC3 (ADGRG3+)": "Endothelial",
    "LEC5 (CLDN11+)": "Endothelial",
    "LEC6 (ADAMTS4+)": "Endothelial",
    "Lymphoid DC": "cDC",
    "LYVE1+ Macrophage": "Macrophages",
    "LZ GC cell": "B cells",
    "Macrophages": "Macrophages",
    "MAIT cell": "CD8 T cell",
    "Mast cell": "Mast cell",
    "Mature arterial EC": "Endothelial",
    "Mature venous EC": "Endothelial",
    "Memory B": "B cells",
    "Mesothelium (PRG4+)": "Mesenchymal",
    "Microfold cell": "Enterocyte",
    "mLN Stroma (FMO2+)": "Mesenchymal",
    "mLTo": "Mesenchymal",
    "MMP9+ Inflammatory macrophage": "Macrophages",
    "Monocytes": "Monocytes",
    "myofibroblast": "Myofibroblast",
    "N cells (NTS+)": "Enteroendocrine cell",
    "Naive B": "B cells",
    "NK cell": "NK cell",
    "NK T cell": "NK T cell",
    "Paneth": "Paneth",
    "pDC": "cDC",
    "Pericyte": "Pericyte",
    "RBC": "Red blood cells",
    "SELL+ CD4 T": "CD4 T cell",
    "SELL+ CD8 T": "Naïve CD8 T cell",
    "STAT1+ Naive B": "B cells",
    "Stem cells": "Stem cell",
    "Stromal 1 (ADAMDEC1+)": "ADAMDEC1+ stromal",
    "Stromal 1 (CCL11+)": "ADAMDEC1+ stromal",
    "Stromal 2 (NPY+)": "NPY+ Stromal",
    "Stromal 3 (C7+)": "ADAMDEC1+ stromal",
    "Stromal 4 (MMP1+)": "ADAMDEC1+ stromal",
    "T reticular": "ADAMDEC1+ stromal",
    "TA": "Stem cell",
    "Tfh": "Tfh",
    "Th1": "CD4 T cell",
    "Th17": "CD4 T cell",
    "Transitional Stromal 3 (C3+)": "Transitional Stromal C3+",
    "Treg": "CD4 T cell",
    "TRGV2 gdT": "NK T cell",
    "Tuft": "Tuft",
    "Fetal venous EC": "Endothelial",
    "TRGV4 gdT": "NK T cell",
    "TRGV5/7 gdT": "NK T cell",
    "myofibroblast (RSPO2+)": "Myofibroblast",
    "LEC4 (STAB2+)": "Endothelial",
    "LEC2 (MADCAM1+)": "Endothelial",
    "Progenitor (NEUROG3+)": "Enteroendocrine cell",
}

# reassign labels and save
adata.obs["category"] = adata.obs["label_spreading_prediction_filtered"].map(
    map_category,
)
adata.obs["curated"] = adata.obs["label_spreading_prediction_filtered"].map(map_curated)

adata.obs["Diagnosis"] = adata.obs["sample_id"].isin(crohns_samples)
adata.obs["Diagnosis"] = (
    adata.obs["Diagnosis"]
    .map({True: "Crohn's Disease", False: "Normal"})
    .astype("category")
)

adata.write_h5ad(PATH / "query_concat_curated.h5ad")

# Joint UMAP on Integrated_05
sc.pp.neighbors(adata, use_rep=REDUCT_NAME)
sc.tl.umap(adata, min_dist=0.3)


def umap_function(color, save_name) -> None:
    """Generate and save UMAP plot colored by specified variable."""
    fig = sc.pl.umap(
        adata,
        color=color,
        frameon=True,
        legend_loc="right margin",
        return_fig=True,
        show=False,
    )

    fig.set_size_inches(6, 6)
    fig.savefig(
        PATH / (save_name + ".pdf"),
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )


# coloured by diagnosis
umap_function(color="Diagnosis", save_name="joint_UMAP_" + REDUCT_NAME + "_diagnosis")

# coloured by sample_id
adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
sample_cats = list(adata.obs["sample_id"].cat.categories)

cmap = plt.get_cmap("tab10", len(sample_cats))
palette = {sid: cmap(i) for i, sid in enumerate(sample_cats)}
adata.uns["sample_id_colors"] = [palette[sid] for sid in sample_cats]

mask = adata.obs["Diagnosis"] == "Crohn's Disease"
adata_crohns = adata[mask].copy()
adata_normal = adata[~mask].copy()

xy = adata.obsm["X_umap"]
pad = 2
xlim = (xy[:, 0].min() - pad, xy[:, 0].max() + pad)
ylim = (xy[:, 1].min() - pad, xy[:, 1].max() + pad)

fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

sc.pl.umap(
    adata_crohns,
    color="sample_id",
    frameon=True,
    legend_loc="right margin",
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

sc.pl.umap(
    adata_normal,
    color="sample_id",
    frameon=True,
    legend_loc="right margin",
    title="Normal",
    ax=axes[1],
    show=False,
)

for ax in axes:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

fig.savefig(
    PATH / f"joint_UMAP_{REDUCT_NAME}_sample_id.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)

# coloured by sample_id
adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
sample_cats = list(adata.obs["sample_id"].cat.categories)

cmap = plt.get_cmap("tab10", len(sample_cats))
palette = {sid: cmap(i) for i, sid in enumerate(sample_cats)}
adata.uns["sample_id_colors"] = [palette[sid] for sid in sample_cats]

mask = adata.obs["Diagnosis"] == "Crohn's Disease"
adata_crohns = adata[mask].copy()
adata_normal = adata[~mask].copy()

xy = adata.obsm["X_umap"]
pad = 2
xlim = (xy[:, 0].min() - pad, xy[:, 0].max() + pad)
ylim = (xy[:, 1].min() - pad, xy[:, 1].max() + pad)

fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

sc.pl.umap(
    adata_crohns,
    color="sample_id",
    frameon=True,
    legend_loc="right margin",
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

sc.pl.umap(
    adata_normal,
    color="sample_id",
    frameon=True,
    legend_loc="right margin",
    title="Normal",
    ax=axes[1],
    show=False,
)

for ax in axes:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

fig.savefig(
    PATH / f"joint_UMAP_{REDUCT_NAME}_sample_id.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)

# coloured by category
adata.obs["category"] = adata.obs["category"].astype("category")
cats = list(adata.obs["category"].cat.categories)

cmap = plt.get_cmap("tab20", len(cats))
cat_palette = {c: cmap(i) for i, c in enumerate(cats)}
adata.uns["category_colors"] = [cat_palette[c] for c in cats]

mask = adata.obs["Diagnosis"] == "Crohn's Disease"
adata_crohns = adata[mask].copy()
adata_normal = adata[~mask].copy()

xy = adata.obsm["X_umap"]
pad = 2
xlim = (xy[:, 0].min() - pad, xy[:, 0].max() + pad)
ylim = (xy[:, 1].min() - pad, xy[:, 1].max() + pad)

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 6),
    gridspec_kw={"width_ratios": [1, 1]},
    constrained_layout=True,
)

# left panel (no legend)
sc.pl.umap(
    adata_crohns,
    color="category",
    frameon=True,
    legend_loc=None,
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

# right panel (with legend)
sc.pl.umap(
    adata_normal,
    color="category",
    frameon=True,
    legend_loc="right margin",
    title="Normal",
    ax=axes[1],
    show=False,
)

for ax in axes:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")

fig.savefig(
    PATH / f"joint_UMAP_{REDUCT_NAME}_category.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)


# coloured by curated
adata.obs["curated"] = adata.obs["curated"].astype("category")
cats = list(adata.obs["curated"].cat.categories)

cmap = sns.color_palette("hls", len(cats))
cat_palette = {c: cmap[i] for i, c in enumerate(cats)}
adata.uns["curated_colors"] = [cat_palette[c] for c in cats]

mask = adata.obs["Diagnosis"] == "Crohn's Disease"
adata_crohns = adata[mask].copy()
adata_normal = adata[~mask].copy()

xy = adata.obsm["X_umap"]
pad = 2
xlim = (xy[:, 0].min() - pad, xy[:, 0].max() + pad)
ylim = (xy[:, 1].min() - pad, xy[:, 1].max() + pad)

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 6),
    gridspec_kw={"width_ratios": [1, 1]},
    constrained_layout=True,
)

# left panel (no legend)
sc.pl.umap(
    adata_crohns,
    color="curated",
    frameon=True,
    legend_loc=None,
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

# right panel (with legend)
sc.pl.umap(
    adata_normal,
    color="curated",
    frameon=True,
    legend_loc="right margin",
    title="Normal",
    ax=axes[1],
    show=False,
)

for ax in axes:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")

fig.savefig(
    PATH / f"joint_UMAP_{REDUCT_NAME}_curated.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)

# coloured by curated and annotated
adata.obs["curated"] = adata.obs["curated"].astype("category")
cats = list(adata.obs["curated"].cat.categories)

cmap = sns.color_palette("hls", len(cats))
cat_palette = {c: cmap[i] for i, c in enumerate(cats)}
adata.uns["curated_colors"] = [cat_palette[c] for c in cats]

mask = adata.obs["Diagnosis"] == "Crohn's Disease"
adata_crohns = adata[mask].copy()
adata_normal = adata[~mask].copy()

xy = adata.obsm["X_umap"]
pad = 2
xlim = (xy[:, 0].min() - pad, xy[:, 0].max() + pad)
ylim = (xy[:, 1].min() - pad, xy[:, 1].max() + pad)

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 6),
    gridspec_kw={"width_ratios": [1, 1]},
    constrained_layout=True,
)

# left panel (no legend)
sc.pl.umap(
    adata_crohns,
    color="curated",
    frameon=True,
    legend_loc="on data",
    legend_fontsize="small",
    legend_fontweight="normal",
    legend_fontoutline=2,
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

# right panel (with legend)
sc.pl.umap(
    adata_normal,
    color="curated",
    frameon=True,
    legend_loc="on data",
    legend_fontsize="small",
    legend_fontweight="normal",
    legend_fontoutline=2,
    title="Normal",
    ax=axes[1],
    show=False,
)

for ax in axes:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")

fig.savefig(
    PATH / f"joint_UMAP_{REDUCT_NAME}_curated_annotated.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)


# Marginal UMAP on Integrated_05
def compute_umap(adata_input, use_rep=REDUCT_NAME, min_dist=0.3, random_state=0):
    """Compute UMAP embedding for given AnnData object."""
    sc.pp.neighbors(adata_input, use_rep=use_rep)
    sc.tl.umap(adata_input, min_dist=min_dist, random_state=random_state)


compute_umap(adata_crohns)
compute_umap(adata_normal)


def set_centered_limits(ax, xy, pad=0.5):
    """Center axis limits around data with equal aspect ratio."""
    x_min, x_max = xy[:, 0].min(), xy[:, 0].max()
    y_min, y_max = xy[:, 1].min(), xy[:, 1].max()
    x_mid, y_mid = (x_min + x_max) / 2, (y_min + y_max) / 2
    span = max(x_max - x_min, y_max - y_min) / 2 + pad
    ax.set_xlim(x_mid - span, x_mid + span)
    ax.set_ylim(y_mid - span, y_mid + span)
    ax.set_aspect("equal")
    return ax


# coloured by sample_id
adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
sample_cats = list(adata.obs["sample_id"].cat.categories)

cmap = plt.get_cmap("tab10", len(sample_cats))
palette = {sid: cmap(i) for i, sid in enumerate(sample_cats)}
adata.uns["sample_id_colors"] = [palette[sid] for sid in sample_cats]


fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

sc.pl.umap(
    adata_crohns,
    color="sample_id",
    frameon=True,
    legend_loc="right margin",
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

sc.pl.umap(
    adata_normal,
    color="sample_id",
    frameon=True,
    legend_loc="right margin",
    title="Normal",
    ax=axes[1],
    show=False,
)

set_centered_limits(axes[0], adata_crohns.obsm["X_umap"])
set_centered_limits(axes[1], adata_normal.obsm["X_umap"])

fig.savefig(
    PATH / f"marginal_UMAP_{REDUCT_NAME}_sample_id.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)

# coloured by category
adata.obs["category"] = adata.obs["category"].astype("category")
cats = list(adata.obs["category"].cat.categories)

cmap = plt.get_cmap("tab20", len(cats))
cat_palette = {c: cmap(i) for i, c in enumerate(cats)}
adata.uns["category_colors"] = [cat_palette[c] for c in cats]

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 6),
    gridspec_kw={"width_ratios": [1, 1]},
    constrained_layout=True,
)

# left panel (no legend)
sc.pl.umap(
    adata_crohns,
    color="category",
    frameon=True,
    legend_loc=None,
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

# right panel (with legend)
sc.pl.umap(
    adata_normal,
    color="category",
    frameon=True,
    legend_loc="right margin",
    title="Normal",
    ax=axes[1],
    show=False,
)

set_centered_limits(axes[0], adata_crohns.obsm["X_umap"])
set_centered_limits(axes[1], adata_normal.obsm["X_umap"])

fig.savefig(
    PATH / f"marginal_UMAP_{REDUCT_NAME}_category.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)


# coloured by curated
adata.obs["curated"] = adata.obs["curated"].astype("category")
cats = list(adata.obs["curated"].cat.categories)

cmap = plt.get_cmap("tab20", len(cats))
cat_palette = {c: cmap(i) for i, c in enumerate(cats)}
adata.uns["curated_colors"] = [cat_palette[c] for c in cats]

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 6),
    gridspec_kw={"width_ratios": [1, 1]},
    constrained_layout=True,
)

# left panel (no legend)
sc.pl.umap(
    adata_crohns,
    color="curated",
    frameon=True,
    legend_loc=None,
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

# right panel (with legend)
sc.pl.umap(
    adata_normal,
    color="curated",
    frameon=True,
    legend_loc="right margin",
    title="Normal",
    ax=axes[1],
    show=False,
)

set_centered_limits(axes[0], adata_crohns.obsm["X_umap"])
set_centered_limits(axes[1], adata_normal.obsm["X_umap"])

fig.savefig(
    PATH / f"marginal_UMAP_{REDUCT_NAME}_curated.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)

# coloured by curated and annotated
adata.obs["curated"] = adata.obs["curated"].astype("category")
cats = list(adata.obs["curated"].cat.categories)

cmap = sns.color_palette("hls", len(cats))
cat_palette = {c: cmap[i] for i, c in enumerate(cats)}
adata.uns["curated_colors"] = [cat_palette[c] for c in cats]

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 6),
    gridspec_kw={"width_ratios": [1, 1]},
    constrained_layout=True,
)

# left panel (no legend)
sc.pl.umap(
    adata_crohns,
    color="curated",
    frameon=True,
    legend_loc="on data",
    legend_fontsize="small",
    legend_fontweight="normal",
    legend_fontoutline=2,
    title="Crohn's Disease",
    ax=axes[0],
    show=False,
)

# right panel (with legend)
sc.pl.umap(
    adata_normal,
    color="curated",
    frameon=True,
    legend_loc="on data",
    legend_fontsize="small",
    legend_fontweight="normal",
    legend_fontoutline=2,
    title="Normal",
    ax=axes[1],
    show=False,
)

set_centered_limits(axes[0], adata_crohns.obsm["X_umap"])
set_centered_limits(axes[1], adata_normal.obsm["X_umap"])

fig.savefig(
    PATH / f"marginal_UMAP_{REDUCT_NAME}_curated_annotated.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2,
)

# celltype proportions per sample
df = adata.obs[["sample_id", "category"]].copy()

# proportions per sample
prop = (
    df.groupby(["sample_id", "category"], observed=True)
    .size()
    .groupby(level=0, observed=True)
    .apply(lambda s: s / s.sum())
    .pivot_table(fill_value=0)  # samples x celltypes
    .sort_index()
)

if isinstance(prop.index, pd.MultiIndex):
    prop = prop.copy()
    prop.index = prop.index.get_level_values(0)
prop.index = prop.index.astype(str)

# color mapping: alphabetical by cell type
cts_alpha = sorted(prop.columns)
cmap = plt.get_cmap("tab10", len(cts_alpha))
color_by_ct = {ct: cmap(i) for i, ct in enumerate(cts_alpha)}

# split panels
crohns_subset = prop.reindex([s for s in crohns_samples if s in prop.index]).dropna(
    how="all",
)
normal_subset = prop.reindex([s for s in normal_samples if s in prop.index]).dropna(
    how="all",
)

width_ratios = [max(1, crohns_subset.shape[0]), max(1, normal_subset.shape[0])]
fig_w = max(10, 0.6 * (crohns_subset.shape[0] + normal_subset.shape[0]) + 2)
fig, axes = plt.subplots(
    1,
    2,
    figsize=(fig_w, 5),
    sharey=True,
    gridspec_kw={"width_ratios": width_ratios},
)

for ax, sub, title in [
    (axes[0], crohns_subset, "Crohn's Disease"),
    (axes[1], normal_subset, "Normal"),
]:
    x = np.arange(sub.shape[0])
    bottom = np.zeros(sub.shape[0])

    # stacked bars in abundance order
    for ct in prop.columns[::-1]:
        vals = sub[ct].to_numpy() if ct in sub.columns else np.zeros(sub.shape[0])
        ax.bar(
            x,
            vals,
            bottom=bottom,
            alpha=0.75,
            label=ct,
            color=color_by_ct[ct],
            edgecolor="dimgray",
            linewidth=1,
            zorder=2,
        )
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(sub.index.tolist(), rotation=0, ha="center")
    ax.set_xlabel("Sample")
    ax.set_title(title)

    ax.set_ylim(0, 1)
    ax.yaxis.set_major_locator(mtick.MultipleLocator(0.25))
    ax.grid(
        axis="y",
        which="major",
        color="black",
        linestyle="-",
        linewidth=0.75,
        zorder=3,
    )
    ax.set_axisbelow(False)  # grid above bars

    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# one legend for both panels, alphabetical order
handles, labels_ = axes[0].get_legend_handles_labels()
order = np.argsort(labels_)
handles = [handles[i] for i in order]
labels_ = [labels_[i] for i in order]
ncol = min(4, max(1, len(labels_) // 10 + 1))
fig.legend(
    handles,
    labels_,
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False,
    ncol=ncol,
)

fig.tight_layout()
out_pdf = PATH / "stacked_bar_category.pdf"
fig.savefig(out_pdf, bbox_inches="tight")


df = adata.obs[["sample_id", "curated"]].copy()

# proportions per sample
prop = (
    df.groupby(["sample_id", "curated"], observed=True)
    .size()
    .groupby(level=0, observed=True)
    .apply(lambda s: s / s.sum())
    .pivot_table(fill_value=0)  # samples x celltypes
    .sort_index()
)

if isinstance(prop.index, pd.MultiIndex):
    prop = prop.copy()
    prop.index = prop.index.get_level_values(0)
prop.index = prop.index.astype(str)

# color mapping: alphabetical by cell type
cts_alpha = sorted(prop.columns)
cmap = sns.color_palette("hls", len(cts_alpha))
color_by_ct = {ct: cmap[i] for i, ct in enumerate(cts_alpha)}

# split panels
crohns_subset = prop.reindex([s for s in crohns_samples if s in prop.index]).dropna(
    how="all",
)
normal_subset = prop.reindex([s for s in normal_samples if s in prop.index]).dropna(
    how="all",
)

width_ratios = [max(1, crohns_subset.shape[0]), max(1, normal_subset.shape[0])]
fig_w = max(10, 0.6 * (crohns_subset.shape[0] + normal_subset.shape[0]) + 2)
fig, axes = plt.subplots(
    1,
    2,
    figsize=(fig_w, 5),
    sharey=True,
    gridspec_kw={"width_ratios": width_ratios},
)

for ax, sub, title in [
    (axes[0], crohns_subset, "Crohn's Disease"),
    (axes[1], normal_subset, "Normal"),
]:
    x = np.arange(sub.shape[0])
    bottom = np.zeros(sub.shape[0])

    # stacked bars in abundance order
    for ct in prop.columns[::-1]:
        vals = sub[ct].to_numpy() if ct in sub.columns else np.zeros(sub.shape[0])
        ax.bar(
            x,
            vals,
            bottom=bottom,
            label=ct,
            color=color_by_ct[ct],
            edgecolor="dimgray",
            linewidth=0.75,
            zorder=2,
        )
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(sub.index.tolist(), rotation=0, ha="center")
    ax.set_xlabel("Sample")
    ax.set_title(title)

    ax.set_ylim(0, 1)
    ax.yaxis.set_major_locator(mtick.MultipleLocator(0.25))
    ax.grid(
        axis="y",
        which="major",
        color="black",
        linestyle="-",
        linewidth=0.75,
        zorder=3,
    )
    ax.set_axisbelow(False)  # grid above bars

    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# one legend for both panels, alphabetical order
handles, labels_ = axes[0].get_legend_handles_labels()
order = np.argsort(labels_)
handles = [handles[i] for i in order]
labels_ = [labels_[i] for i in order]
ncol = min(4, max(1, len(labels_) // 10 + 1))
fig.legend(
    handles,
    labels_,
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False,
    ncol=ncol,
)

fig.tight_layout()
out_pdf = PATH / "stacked_bar_curated.pdf"
fig.savefig(out_pdf, bbox_inches="tight")
