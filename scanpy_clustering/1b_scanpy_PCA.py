import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from utils.scanpy_pca_utils import pca_dispersion, pca_varfeatures

SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy"

# Load existing AnnData object
adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))

# Create output directory
PCA_OUTPUT_PATH = os.path.join(SCANPY_OBJECT_PATH, "PCA_stats")
os.makedirs(PCA_OUTPUT_PATH) if not os.path.exists(PCA_OUTPUT_PATH) else None

# Parameters for variable gene selection
disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
xmin = 0.1
xmax = 10

dr_list = []
vfeature_objects = {}

# Loop over dispersion cutoffs using mean.var.plot equivalent
pca_dispersion(adata, PCA_OUTPUT_PATH, disps, xmin, xmax, dr_list, vfeature_objects)

# Loop over nfeatures with vst method
pca_varfeatures(adata, PCA_OUTPUT_PATH, n_features, dr_list, vfeature_objects)

# Plot elbow plots
for dr in dr_list:
    var_ratio = adata.uns[f"{dr}_pca"]["variance_ratio"]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(18, 5))
    axs.plot(np.arange(1, len(var_ratio) + 1), var_ratio, marker="o")
    axs.set_xlabel("PC")
    axs.set_ylabel("Variance ratio")
    axs.set_title(dr, fontsize=10)
    axs.set_xticks(np.arange(1, len(var_ratio) + 1))
    plt.tight_layout()
    plt.savefig(os.path.join(PCA_OUTPUT_PATH, f"elbow_plots_{dr}.png"), dpi=300)
    plt.close()

# Plot loadings (PC1-3)
pcs = adata.varm[f"{dr}_PCs"]  # shape: (n_genes, n_pcs)
genes = adata.var_names

# Number of PCs and top genes to plot
n_pcs_to_plot = 5
top_n = 10

# Create subplots
for dr in dr_list:
    pcs = adata.varm[f"{dr}_PCs"]  # shape: (n_genes, n_pcs)
    genes = adata.var_names

    n_pcs_to_plot = 3
    top_n = 10

    fig, axs = plt.subplots(nrows=n_pcs_to_plot, figsize=(10, 3 * n_pcs_to_plot))

    for i in range(n_pcs_to_plot):
        pc_loadings = pcs[:, i]
        loading_df = pd.DataFrame({"gene": genes, "loading": pc_loadings}).set_index(
            "gene"
        )

    # Get top positive and negative contributing genes
    top_genes = pd.concat(
        [
            loading_df.sort_values("loading", ascending=False).head(top_n),
            loading_df.sort_values("loading", ascending=True).head(top_n),
        ]
    )

    sns.barplot(
        x="loading",
        y=top_genes.index,
        data=top_genes.reset_index(),
        ax=axs[i],
        palette="vlag",
    )
    axs[i].axvline(0, color="black", linestyle="--")
    axs[i].set_title(f"{dr} | PC{i + 1}", fontsize=12)
    axs[i].set_xlabel("Loading weight")

plt.tight_layout()
plt.savefig(os.path.join(PCA_OUTPUT_PATH, f"{dr}_pca_loadings.png"), dpi=300)
plt.close()

# Variable feature plots
fig, axs = plt.subplots(nrows=1, ncols=min(len(vfeature_objects), 5), figsize=(18, 5))
axs = axs.flatten()

for i, (name, features) in enumerate(vfeature_objects.items()):
    sc.pl.highly_variable_genes(adata, ax=axs[i], show=False)
    axs[i].set_title(name, fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(PCA_OUTPUT_PATH, "v_features.png"), dpi=300)
plt.close()

# Save the AnnData object
adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_reduced.h5ad"))
