import scanpy as sc
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
import os
import numpy as np
import pandas as pd

PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/Integrated_05_label" 
label =  "Integrated_05"
REDUCT_NAME = "X_embeddings"+label
adata = sc.read_h5ad(os.path.join(PATH, "query_concat.h5ad"))

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
    "Cycling plasma cell": "B cells",
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
}
map_curated = {
    "Activated CD4 T": "CD4 T cell",
    "Activated CD8 T": "CD8 T cell",
    "Adult Glia": "Neuronal",
    "arterial capillary": "Endothelial",
    "BEST2+ Goblet cell": "Goblet",
    "BEST4+ epithelial": "BEST4+ epitheial",
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
    "Stem cells": "Progenitor epithelial",
    "Stromal 1 (ADAMDEC1+)": "ADAMDEC1+ stromal",
    "Stromal 1 (CCL11+)": "ADAMDEC1+ stromal",
    "Stromal 2 (NPY+)": "NPY+ Stromal",
    "Stromal 3 (C7+)": "ADAMDEC1+ stromal",
    "Stromal 4 (MMP1+)": "ADAMDEC1+ stromal",
    "T reticular": "ADAMDEC1+ stromal",
    "TA": "Progenitor epithelial",
    "Tfh": "CD4 T cell",
    "Th1": "CD4 T cell",
    "Th17": "CD4 T cell",
    "Transitional Stromal 3 (C3+)": "Transitional Stromal C3+",
    "Treg": "CD4 T cell",
    "TRGV2 gdT": "NK T cell",
    "Tuft": "Progenitor epithelial",
}

# reassign labels
adata.obs["category"] = adata.obs["Integrated_05"].map(map_category)
adata.obs["curated"]  = adata.obs["Integrated_05"].map(map_curated)


# Joint UMAP on Integrated_05 coloured by diagnosis
# Joint UMAP on Integrated_05 coloured by sample_id
# Joint UMAP on Integrated_05 coloured by category
# Joint UMAP on Integrated_05 coloured by curated

# Marginal UMAP on Integrated_05 coloured by diagnosis
# Marginal UMAP on Integrated_05 coloured by sample_id
# Marginal UMAP on Integrated_05 coloured by category

# Marginal UMAP on Integrated_05 coloured by sample_id
# Marginal UMAP on Integrated_05 coloured by curated



print("performing UMAP")
sc.pp.neighbors(adata, use_rep=REDUCT_NAME)
sc.tl.umap(adata, min_dist=0.3)

plot = sc.pl.umap(
    adata,
    color="category",
    frameon=True,
    legend_loc="on data",   # put labels on clusters if possible
    legend_fontsize=6,
    title="category",
    show=True,
    save=None,
    return_fig=True
)
plot.set_size_inches(6, 6)
plot.savefig(
    os.path.join(PATH, "UMAP_" + REDUCT_NAME + "_celltype_category.pdf"),
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

plot = sc.pl.umap(
    adata,
    color="curated",
    frameon=True,
    legend_loc="on data",   # put labels on clusters if possible
    legend_fontsize=8,
    title="curate celltypes",
    show=True,
    save=None,
    return_fig=True
)
plot.set_size_inches(6, 6)
plot.savefig(
    os.path.join(PATH, "UMAP_" + REDUCT_NAME + "_celltype_curated.pdf"),
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)


# per sample
selected = ["C_03", "C_08", "C_13"]
mask = adata.obs["sample_id"].isin(selected)

# Subset into selected and complement
adata_sel = adata[mask].copy()
adata_rest = adata[~mask].copy()

# UMAP needs neighbors for each subset
for ad in (adata_sel, adata_rest):
    sc.pp.neighbors(ad, use_rep=REDUCT_NAME)
    sc.tl.umap(ad, min_dist=0.3)

# Create side-by-side plots
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

sc.pl.umap(
    adata_sel,
    color="category",
    frameon=True,
    legend_loc="on data",
    legend_fontsize=6,
    title="Selected samples (C_03, C_08, C_13)",
    show=False,
    ax=axes[0]
)

sc.pl.umap(
    adata_rest,
    color="category",
    frameon=True,
    legend_loc="on data",
    legend_fontsize=6,
    title="Complement samples",
    show=False,
    ax=axes[1]
)

fig.tight_layout()
fig.savefig(
    os.path.join(PATH, f"UMAP_{REDUCT_NAME}_category_selected_vs_rest.pdf"),
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2
)

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

sc.pl.umap(
    adata_sel,
    color="curated",
    frameon=True,
    legend_loc="on data",
    legend_fontsize=6,
    title="Selected samples (C_03, C_08, C_13)",
    show=False,
    ax=axes[0]
)

sc.pl.umap(
    adata_rest,
    color="curated",
    frameon=True,
    legend_loc="on data",
    legend_fontsize=6,
    title="Complement samples",
    show=False,
    ax=axes[1]
)

fig.tight_layout()
fig.savefig(
    os.path.join(PATH, f"UMAP_{REDUCT_NAME}_curated_selected_vs_rest.pdf"),
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.2
)