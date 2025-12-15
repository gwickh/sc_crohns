import os

import numpy as np
import scanpy as sc


def run_pca(hvgs, adata, method_name, dr_name, PCA_OUTPUT_PATH, vfeature_objects):
    """
    Run PCA on set of HVG and store results in adata.
    """
    hv_adata = adata[:, hvgs].copy()

    features = hv_adata.var_names.tolist()
    vfeature_objects[method_name + "_vfeatures"] = features

    # Save features
    with open(os.path.join(PCA_OUTPUT_PATH, f"{method_name}_features.txt"), "w") as f:
        f.write("\n".join(features))

    # PCA
    sc.tl.pca(hv_adata, svd_solver="auto", n_comps=50, random_state=42)
    adata.obsm[f"X_{dr_name}"] = hv_adata.obsm["X_pca"]
    adata.uns[f"{dr_name}_pca"] = hv_adata.uns["pca"]

    n_genes_total = adata.n_vars
    n_pcs = hv_adata.varm["PCs"].shape[1]
    full_pcs = np.full((n_genes_total, n_pcs), np.nan)

    # Get indices of HVGs in adata.var_names
    hvg_mask = adata.var_names.isin(hv_adata.var_names)
    full_pcs[hvg_mask, :] = hv_adata.varm["PCs"]

    # Store in adata.varm
    adata.varm[f"{dr_name}_PCs"] = full_pcs


def pca_dispersion(
    adata, PCA_OUTPUT_PATH, disps, xmin, xmax, dr_list, vfeature_objects
):
    """
    # Loop over dispersion cutoffs using mean.var.plot equivalent
    """
    for disp in disps:
        method_name = f"mean.var.plot_disp_{disp}"
        dr_name = f"pca_{method_name}"
        dr_list.append(dr_name)

        sc.pp.highly_variable_genes(
            adata, flavor="seurat", min_mean=xmin, max_mean=xmax, min_disp=disp
        )

        # Select top 2000 by dispersion
        hvgs = (
            adata.var[adata.var["highly_variable"]]
            .sort_values(by="dispersions_norm", ascending=False)
            .head(2000)
            .index
        )
        run_pca(hvgs, adata, method_name, dr_name, PCA_OUTPUT_PATH, vfeature_objects)


def pca_varfeatures(
    adata, PCA_OUTPUT_PATH, n_features, dr_list, vfeature_objects
) -> None:
    for n in n_features:
        method_name = f"vst_top_{n}"
        dr_name = f"pca_{method_name}"
        dr_list.append(dr_name)

        sc.pp.highly_variable_genes(adata, n_top_genes=n, batch_key="sample_id")

        hvgs = adata.var.highly_variable

        run_pca(hvgs, adata, method_name, dr_name, PCA_OUTPUT_PATH, vfeature_objects)
