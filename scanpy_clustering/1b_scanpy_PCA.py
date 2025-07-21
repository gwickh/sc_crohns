import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/scanpy_clustering_output"

# Load existing AnnData object
adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))

# Create output directory
PCA_OUTPUT_PATH = os.path.join(SCANPY_OBJECT_PATH, "PCA_stats")
os.makedirs(PCA_OUTPUT_PATH, exist_ok=True)

# Parameters for variable gene selection
disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
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
        min_disp=disp
    )
    
    hv_adata = adata[:, adata.var.highly_variable].copy()

    # Select top 2000 by dispersion
    top_2000 = (
        adata.var[adata.var['highly_variable']]
        .sort_values(by='dispersions_norm', ascending=False)
        .head(2000).index
    )

    hv_adata = adata[:, top_2000].copy()

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
    adata.uns[f'{dr_name}_pca'] = hv_adata.uns['pca']
    
    n_genes_total = adata.n_vars
    n_pcs = hv_adata.varm['PCs'].shape[1]
    full_pcs = np.full((n_genes_total, n_pcs), np.nan)

    # Get indices of HVGs in adata.var_names
    hvg_mask = adata.var_names.isin(hv_adata.var_names)
    full_pcs[hvg_mask, :] = hv_adata.varm['PCs']

    # Store in adata.varm
    adata.varm[f'{dr_name}_PCs'] = full_pcs

# Loop over nfeatures with vst method
for n in n_features:
    method_name = f"vst_top_{n}"
    dr_name = f"pca_{method_name}"
    dr_li.append(dr_name)

    sc.pp.highly_variable_genes(adata, n_top_genes=n, batch_key="sample_id")
    hv_adata = adata[:, adata.var.highly_variable].copy()
    features = hv_adata.var_names.tolist()
    vfeature_objects[method_name + "_vfeatures"] = features

    with open(os.path.join(PCA_OUTPUT_PATH, f"{method_name}_features.txt"), "w") as f:
        f.write("\n".join(features))

    sc.tl.pca(hv_adata, svd_solver='auto', n_comps=50, random_state=42)
    hv_adata.obsm[f"X_{dr_name}"] = hv_adata.obsm["X_pca"]
    adata.obsm[f"X_{dr_name}"] = hv_adata.obsm["X_pca"]
    adata.uns[f'{dr_name}_pca'] = hv_adata.uns['pca']
    
    n_genes_total = adata.n_vars
    n_pcs = hv_adata.varm['PCs'].shape[1]
    full_pcs = np.full((n_genes_total, n_pcs), np.nan)

    # Get indices of HVGs in adata.var_names
    hvg_mask = adata.var_names.isin(hv_adata.var_names)
    full_pcs[hvg_mask, :] = hv_adata.varm['PCs']

    # Store in adata.varm
    adata.varm[f'{dr_name}_PCs'] = full_pcs

# Plot elbow plots
for dr in dr_li:
    var_ratio = adata.uns[f'{dr}_pca']['variance_ratio']
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(18, 5))
    axs.plot(np.arange(1, len(var_ratio)+1), var_ratio, marker='o')
    axs.set_xlabel('PC')
    axs.set_ylabel('Variance ratio')
    axs.set_title(dr, fontsize=10)
    axs.set_xticks(np.arange(1, len(var_ratio)+1))
    plt.tight_layout()
    plt.savefig(os.path.join(PCA_OUTPUT_PATH, f"elbow_plots_{dr}.png"), dpi=300)
    plt.close()

# Plot loadings (PC1-3)
pcs = adata.varm[f'{dr}_PCs']   # shape: (n_genes, n_pcs)
genes = adata.var_names

# Number of PCs and top genes to plot
n_pcs_to_plot = 5
top_n = 10

# Create subplots
for dr in dr_li:
    pcs = adata.varm[f'{dr}_PCs']   # shape: (n_genes, n_pcs)
    genes = adata.var_names

    n_pcs_to_plot = 3
    top_n = 10

    fig, axs = plt.subplots(nrows=n_pcs_to_plot, figsize=(10, 3 * n_pcs_to_plot))

    for i in range(n_pcs_to_plot):
        pc_loadings = pcs[:, i]
        loading_df = pd.DataFrame({
            'gene': genes,
            'loading': pc_loadings
        }).set_index('gene')

    # Get top positive and negative contributing genes
    top_genes = pd.concat([
        loading_df.sort_values('loading', ascending=False).head(top_n),
        loading_df.sort_values('loading', ascending=True).head(top_n)
    ])

    sns.barplot(
        x='loading', y=top_genes.index, data=top_genes.reset_index(),
        ax=axs[i], palette="vlag"
    )
    axs[i].axvline(0, color='black', linestyle='--')
    axs[i].set_title(f"{dr} | PC{i+1}", fontsize=12)
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
adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))
