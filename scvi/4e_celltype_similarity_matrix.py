import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

# vars
PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/Integrated_05_label"
MIN_CELLS = 10
SOFT_Q_LAB = os.path.join(PATH, "query_soft_labelsIntegrated_05.csv")
GCA_Integrated_05 = os.path.join(PATH, "gca_ref_scvi_labels.csv")
GCA_category = os.path.join(PATH, "../category_label/gca_ref_scvi_labels.csv")
OUT_PREFIX = os.path.join(PATH, "query_soft_labelsIntegrated_05")

os.makedirs(PATH, exist_ok=True)

# load soft label matrix
df_raw = pd.read_csv(SOFT_Q_LAB, sep="\t", header=None)
soft = df_raw.iloc[1:, 1:].copy()
soft.index = df_raw.iloc[1:, 0].astype(str)
soft.columns = df_raw.iloc[0, 1:].astype(str)
soft = soft.apply(pd.to_numeric, errors="coerce").fillna(0.0)

# Filter rare labels (by total soft count across cells)
soft = soft.loc[:, soft.sum(axis=0) >= MIN_CELLS]
labels = soft.columns.tolist()

if len(labels) < 2:
    raise ValueError(
        f"After filtering MIN_CELLS={MIN_CELLS}, only {len(labels)} label(s) remain. "
        f"Lower MIN_CELLS or check input."
    )

# fuzzy overlap
A = soft.to_numpy(dtype=float)  # (cells x labels)
n_cells, n_labels = A.shape

# Fuzzy overlap M[i,j] = sum_c min(A[c,i], A[c,j])
M = np.zeros((n_labels, n_labels), dtype=float)
for i in range(n_labels):
    ai = A[:, i]
    # loop j >= i and mirror (symmetric)
    for j in range(i, n_labels):
        v = float(np.minimum(ai, A[:, j]).sum())
        M[i, j] = M[j, i] = v

# Normalize with Ochiai / cosine on counts: M_norm = M / sqrt(v_i * v_j)
v = A.sum(axis=0)  # per label total
den = np.sqrt(np.outer(v, v))
den[den == 0] = np.nan
M_norm = M / den

# save raw outputs
pd.DataFrame(M, index=labels, columns=labels).to_csv(f"{OUT_PREFIX}.tsv", sep="\t")
pd.DataFrame(M_norm, index=labels, columns=labels).to_csv(f"{OUT_PREFIX}_norm.tsv", sep="\t")

# hierarchical clustering (on M_norm)
S = np.array(M_norm, dtype=float)   # Convert similarity (in [0,1]) to distance
np.fill_diagonal(S, 1.0)            # Ensure diagonals are 1 for a proper similarity
S = np.clip(S, 0.0, 1.0)

D = 1.0 - S  # distance in [0,1], symmetric, 0 on diagonal

# Condensed distance vector for linkage
# (squareform reads upper triangle by default; checks=False for speed since we know it's valid)
Z = linkage(squareform(D, checks=False), method="average")  # UPGMA; try 'complete' if you prefer tighter clusters
order = leaves_list(Z)

# Reorder labels and matrices
labels_ord = [labels[i] for i in order]
M_df = pd.DataFrame(M, index=labels, columns=labels).loc[labels_ord, labels_ord]
Mnorm_df = pd.DataFrame(M_norm, index=labels, columns=labels).loc[labels_ord, labels_ord]
labels = labels_ord  # reuse downstream
n = len(labels)

# 2-panel masked heatmaps (clustered)
mask_upper = np.triu(np.ones((n, n), dtype=bool), k=1)

# Log norm for M
vmin_M = max(M_df.values[M_df.values > 0].min(), 1e-6) if (M_df.values > 0).any() else 1e-6
vmax_M = M_df.values.max()
norm_M = mcolors.LogNorm(vmin=vmin_M, vmax=vmax_M)

# Power-stretch for M_norm
vmax_Mn = np.nanmax(Mnorm_df.values)
norm_Mn = mcolors.PowerNorm(gamma=0.1, vmin=0, vmax=vmax_Mn)

fig, axes = plt.subplots(1, 2, figsize=(24, 9))

sns.heatmap(
    M_df, ax=axes[0], cmap="magma",
    norm=norm_M, mask=mask_upper, square=True,
    xticklabels=labels, yticklabels=labels,
    cbar_kws={"label": "Fuzzy overlap (M)"}
)
axes[0].set_title("Fuzzy overlap of cell types (clustered)")
axes[0].tick_params(axis="x", labelrotation=90, labelsize=8)
axes[0].tick_params(axis="y", labelsize=8)

sns.heatmap(
    Mnorm_df, ax=axes[1], cmap="viridis",
    norm=norm_Mn, mask=mask_upper, square=True,
    xticklabels=labels, yticklabels=labels,
    cbar_kws={"label": "Normalized similarity (Ochiai, γ=0.1)"}
)
axes[1].set_title("Normalized Ochiai similarity (clustered)")
axes[1].tick_params(axis="x", labelrotation=90, labelsize=8)
axes[1].tick_params(axis="y", labelsize=8)

plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}_heatmaps_clustered.pdf", bbox_inches="tight")
plt.close()