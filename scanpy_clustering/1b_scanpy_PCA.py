import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/scanpy_clustering_output"

# Load existing AnnData object
adata_path = os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad")
adata = sc.read_h5ad(adata_path)

# Create output directory
PCA_OUTPUT_PATH = os.path.join(SCANPY_OBJECT_PATH, "PCA_stats")
os.makedirs(PCA_OUTPUT_PATH, exist_ok=True)

# Parameters for variable gene selection
disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
res = [0.2, 0.5, 0.8, 1.0, 1.2]
neighbors = [10, 20, 30, 50]
xmin = 0.1
xmax = 10

dr_li = []
vfeature_objects = {}

# Loop over dispersion cutoffs using mean.var.plot equivalent
for disp in disps:
    method_name = f"mean.var.plot_disp_{disp}"
    dr_name = f"pca_{method_name}"
    dr_li.append(dr_name)

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        min_mean=xmin,
        max_mean=xmax,
        min_disp=disp,
    )
    hv_adata = adata[:, adata.var.highly_variable].copy()
    features = hv_adata.var_names.tolist()
    vfeature_objects[method_name + "_vfeatures"] = features

    # Save features
    with open(os.path.join(PCA_OUTPUT_PATH, f"{method_name}_features.txt"), "w") as f:
        f.write("\n".join(features))

    # PCA
    sc.tl.pca(
        hv_adata, 
        svd_solver='auto', 
        n_comps=50,
        random_state=42
    )
    hv_adata.obsm[f"X_{dr_name}"] = hv_adata.obsm["X_pca"]
    adata.obsm[f"X_{dr_name}"] = hv_adata.obsm["X_pca"]

# Loop over nfeatures with vst method
for n in n_features:
    method_name = f"vst_top_{n}"
    dr_name = f"pca_{method_name}"
    dr_li.append(dr_name)

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n)
    hv_adata = adata[:, adata.var.highly_variable].copy()
    features = hv_adata.var_names.tolist()
    vfeature_objects[method_name + "_vfeatures"] = features

    with open(os.path.join(PCA_OUTPUT_PATH, f"{method_name}_features.txt"), "w") as f:
        f.write("\n".join(features))

    sc.tl.pca(hv_adata, svd_solver='auto', n_comps=50, random_state=42)
    hv_adata.obsm[f"X_{dr_name}"] = hv_adata.obsm["X_pca"]
    adata.obsm[f"X_{dr_name}"] = hv_adata.obsm["X_pca"]

# Plot elbow plots
fig, axs = plt.subplots(nrows=int(np.ceil(len(dr_li)/2)), ncols=2, figsize=(18, 4 * int(np.ceil(len(dr_li)/2))))
axs = axs.flatten()

for i, dr in enumerate(dr_li):
    sc.pl.pca_variance_ratio(adata, n_pcs=50, ax=axs[i], obsm=f"X_{dr}", show=False)
    axs[i].set_title(dr, fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(PCA_OUTPUT_PATH, "elbow_plots.png"), dpi=300)
plt.close()

# Plot loadings (PC1-3)
loading_plots = []
for dr in dr_li:
    sc.pl.pca_loadings(adata, components=[1, 2, 3], obsm=f"X_{dr}", show=False, title=dr)

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
adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))
