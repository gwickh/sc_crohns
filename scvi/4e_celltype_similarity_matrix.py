import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
import seaborn as sns


PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/Integrated_05_label"
MIN_CELLS = 10
SOFT_Q_LAB = os.path.join(PATH, "query_soft_labelsIntegrated_05.csv")
GCA_Integrated_05 = os.path.join(PATH, "gca_ref_scvi_labels.csv")
GCA_category = os.path.join(PATH, "../category_label/gca_ref_scvi_labels.csv")
OUT_PREFIX = os.path.join(PATH, "query_soft_labelsIntegrated_05")

### generate similarity matrix ###
# load soft matrix
df_raw = pd.read_csv(SOFT_Q_LAB, sep="\t", header=None)
soft = df_raw.iloc[1:, 1:].copy()
soft.index = df_raw.iloc[1:, 0].astype(str)
soft.columns = df_raw.iloc[0, 1:].astype(str)
soft = soft.apply(pd.to_numeric, errors="coerce").fillna(0.0)

# filter rare labels
soft = soft.loc[:, soft.sum(axis=0) >= MIN_CELLS]
labels = soft.columns.tolist()

# fuzzy overlap M[i,j]
A = soft.to_numpy(dtype=float)              # (cells x labels)
n_cells, n_labels = A.shape
M = np.zeros((n_labels, n_labels), float)
for i in range(n_labels):
    ai = A[:, i]
    for j in range(i, n_labels):
        v = float(np.minimum(ai, A[:, j]).sum())
        M[i, j] = M[j, i] = v

# normalize: M_norm = M / sqrt(v_i * v_j)
v = A.sum(axis=0)
den = np.sqrt(np.outer(v, v))
den[den == 0] = np.nan
M_norm = M / den

# map labels to categories from GCA_STATS (cols: Integrated_05, category)
gca_Integrated_05 = pd.read_csv(GCA_Integrated_05, sep="\t", header=0)
gca_category = pd.read_csv(GCA_category, sep="\t", header=0)

map_cat = dict(zip(gca_Integrated_05.iloc[:, 1].astype(str), gca_category.iloc[:, 1].astype(str)))
categories = [map_cat.get(l, "") for l in labels]
cats = pd.Series(categories, index=labels, name="category")

# write outputs
pd.DataFrame(M, index=labels, columns=labels).to_csv(f"{OUT_PREFIX}.tsv", sep="\t")
pd.DataFrame(M_norm, index=labels, columns=labels).to_csv(f"{OUT_PREFIX}_norm.tsv", sep="\t")
cats.to_csv(f"{OUT_PREFIX}_label_map.tsv", sep="\t")   

### plot normalised heatmap ###
# Convert numpy arrays to DataFrames with labels
M = pd.DataFrame(M, index=labels, columns=labels)
M_norm = pd.DataFrame(M_norm, index=labels, columns=labels)

# align order
M      = M.loc[M.index, M.index]
M_norm = M_norm.loc[M.index, M.index]
labels = M.index.astype(str).tolist()
n = len(labels)

# mask
mask_upper   = np.triu(np.ones((n, n), dtype=bool), k=1)

# scales
vmin = max(M.values[M.values > 0].min(), 1e-6) if (M.values > 0).any() else 1e-6
vmax = M.values.max()
norm_power = mcolors.PowerNorm(gamma=0.1, vmin=0, vmax=np.nanmax(M_norm.values))

fig, axes = plt.subplots(1, 2, figsize=(24, 9))

# M: log-scaled
sns.heatmap(
    M, ax=axes[0], cmap="magma",
    norm=mcolors.LogNorm(vmin=vmin, vmax=vmax),
    mask=mask_upper, square=True,
    xticklabels=labels, yticklabels=labels, 
    cbar_kws={"label": "Fuzzy overlap (M)"}
)
axes[0].set_title("Fuzzy overlap of cell types")
axes[0].tick_params(axis="x", labelrotation=90, labelsize=8)
axes[0].tick_params(axis="y", labelsize=8)

# M_norm: power-stretched
sns.heatmap(
    M_norm, ax=axes[1], cmap="viridis",
    norm=norm_power, mask=mask_upper, square=True,
    xticklabels=labels, yticklabels=labels,
    cbar_kws={"label": "Normalized similarity (Power γ=0.1)"}
)
axes[1].set_title("Normalized Ochiai similarity")
axes[1].tick_params(axis="x", labelrotation=90, labelsize=8)
axes[1].tick_params(axis="y", labelsize=8)

plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}_heatmaps.pdf", bbox_inches="tight")
plt.close()