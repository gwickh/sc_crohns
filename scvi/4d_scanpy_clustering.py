import scanpy as sc
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
import os
import numpy as np
import pandas as pd

PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/Integrated_05_label"
label =  "Integrated_05" # "category" 
adata = sc.read_h5ad(os.path.join(PATH, "query_concat.h5ad"))
REDUCT_NAME = "X_embeddings"+label

# use scVI latent space for UMAP generation of full dataset coloured by high granularity labels
sc.pp.neighbors(adata, use_rep=REDUCT_NAME)
sc.tl.umap(adata, min_dist=0.3)

plot = sc.pl.umap(
    adata,
    color=label,
    ncols=2,
    frameon=True, 
    return_fig=True,
    show=False
)
plot.set_size_inches(6, 6)
plot.savefig(
    os.path.join(PATH, "UMAP_" + REDUCT_NAME + "_celltypes.pdf"),
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

plot = sc.pl.umap(
    adata,
    color="sample_id",
    ncols=4,
    frameon=True, 
    return_fig=True,
    show=False
)
plot.set_size_inches(6, 6)
plot.savefig(
    os.path.join(PATH, "UMAP_" + REDUCT_NAME + "_metadata.pdf"), 
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

# use scVI latent space for UMAP generation of Crohn's samples
selected = ["C_03", "C_08", "C_13"]
mask = adata.obs["sample_id"].isin(selected)

sc.pp.neighbors(adata[mask], use_rep=REDUCT_NAME)
sc.tl.umap(adata[mask], min_dist=0.3)

plot = sc.pl.umap(
    adata[mask],
    color=label,
    ncols=2,
    frameon=True, 
    return_fig=True,
    show=False
)
plot.set_size_inches(6, 6)
plot.savefig(
    os.path.join(PATH, "UMAP_Crohns" + REDUCT_NAME + "_celltypes.pdf"), 
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

plot = sc.pl.umap(
    adata[mask],
    color="sample_id",
    ncols=4,
    frameon=True, 
    return_fig=True,
    show=False
)
plot.set_size_inches(6, 6)
plot.savefig(
    os.path.join(PATH, "UMAP_Crohns_" + REDUCT_NAME + "_metadata.pdf"), 
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

# use scVI latent space for UMAP generation of normal samples
selected = ["N_04", "N_06", "N_07", "N_12"]
mask = adata.obs["sample_id"].isin(selected)

sc.pp.neighbors(adata[mask], use_rep=REDUCT_NAME)
sc.tl.umap(adata[mask], min_dist=0.3)

plot = sc.pl.umap(
    adata[mask],
    color=label,
    ncols=2,
    frameon=True, 
    return_fig=True,
    show=False
)
plot.set_size_inches(6, 6)
plot.savefig(
    os.path.join(PATH, "UMAP_Normal" + REDUCT_NAME + "_celltypes.pdf"),
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

plot = sc.pl.umap(
    adata[mask],
    color="sample_id",
    ncols=4,
    frameon=True, 
    return_fig=True,
    show=False
)
plot.set_size_inches(6, 6)
plot.savefig(
    os.path.join(PATH, "UMAP_Normal" + REDUCT_NAME + "_metadata.pdf"), 
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)


# celltype proportions per sample
assert "sample_id" in adata.obs, "sample_id not in adata.obs"
assert label in adata.obs, label

df = adata.obs[["sample_id", label]].copy()

# proportions per sample
prop = (
    df.groupby(["sample_id", label], observed=True).size()
      .groupby(level=0, observed=True).apply(lambda s: s / s.sum())
      .unstack(fill_value=0)    # samples x celltypes
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
group_A = ["C_03", "C_08", "C_13"]
group_B = ["N_04", "N_06", "N_07", "N_12"]
subA = prop.reindex([s for s in group_A if s in prop.index]).dropna(how="all")
subB = prop.reindex([s for s in group_B if s in prop.index]).dropna(how="all")

width_ratios = [max(1, subA.shape[0]), max(1, subB.shape[0])]
fig_w = max(10, 0.6 * (subA.shape[0] + subB.shape[0]) + 2)
fig, axes = plt.subplots(1, 2, figsize=(fig_w, 5), sharey=True,
                         gridspec_kw={"width_ratios": width_ratios})

for ax, sub, title in [(axes[0], subA, "Crohn's Disease Terminal Ileum"), (axes[1], subB, "Normal Terminal Ileum")]:
    x = np.arange(sub.shape[0])
    bottom = np.zeros(sub.shape[0])

    # stacked bars in abundance order 
    for ct in prop.columns:
        vals = sub[ct].values if ct in sub.columns else np.zeros(sub.shape[0])
        ax.bar(x, vals, bottom=bottom, label=ct, color=color_by_ct[ct], linewidth=0)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(sub.index.tolist(), rotation=0, ha="center")
    ax.set_xlabel("Sample")
    ax.set_title(title)

    ax.set_ylim(0, 1)
    ax.yaxis.set_major_locator(mtick.MultipleLocator(0.25))
    ax.grid(axis="y", which="major", color="black", linestyle="-", linewidth=1.2, zorder=3)
    ax.set_axisbelow(False)  # grid above bars

    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# one legend for both panels, alphabetical order
handles, labels_ = axes[0].get_legend_handles_labels()
order = np.argsort(labels_)
handles = [handles[i] for i in order]
labels_  = [labels_[i] for i in order]
ncol = min(4, max(1, len(labels_) // 10 + 1))
fig.legend(handles, labels_, bbox_to_anchor=(1.02, 0.5), loc="center left",
           frameon=False, ncol=ncol, title=label)

fig.tight_layout()
out_pdf = os.path.join(PATH, "stacked_bar_celltypes_by_sample_panels.pdf")
fig.savefig(out_pdf, bbox_inches="tight")