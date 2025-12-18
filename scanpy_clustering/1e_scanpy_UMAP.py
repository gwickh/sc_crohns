import os

import anndata as an
import matplotlib
import scanpy as sc

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy"
adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_clustered.h5ad"))

UMAP_PATH = os.path.join(SCANPY_OBJECT_PATH, "UMAP_plots")
os.makedirs(UMAP_PATH) if not os.path.exists(UMAP_PATH) else None


def plot_umaps_for_grid(
    adata,
    disps,
    n_features,
    neighbors,
    res,
    UMAP_PATH,
    random_state=0,
):
    """
    Loop over the exact same DR/k/res grid as run_clustering and save UMAPs
    colored by cluster, sample_id
    """

    # backup global X_umap since sc.tl.umap writes there
    had_umap = "X_umap" in adata.obsm
    umap_backup = adata.obsm["X_umap"].copy() if had_umap else None

    def _do_one_dr(dr):
        dr_out = os.path.join(UMAP_PATH, dr)
        os.makedirs(dr_out, exist_ok=True)

        for k in neighbors:
            neighbors_key = f"neighbors_{dr}_k{k}"
            if neighbors_key not in adata.uns:
                print(f"[skip] Missing neighbors: {neighbors_key}")
                continue

            # Compute one UMAP per (dr, k)
            basis = f"umap_{dr}_k{k}"
            umap_key = f"X_{basis}"

            if umap_key not in adata.obsm:
                sc.tl.umap(
                    adata, neighbors_key=neighbors_key, random_state=random_state
                )
                adata.obsm[umap_key] = adata.obsm["X_umap"].copy()

            # Plot for each resolution (cluster labels differ)
            for r in res:
                cl_key = f"clusters_{dr}_k{k}_r{r}"
                if cl_key not in adata.obs:
                    print(f"[skip] Missing clusters: {cl_key}")
                    continue

                out_png = os.path.join(dr_out, f"{basis}__r{r}.png")
                title = f"{dr} | k={k} | r={r}"

                sc.pl.embedding(
                    adata,
                    basis=basis,
                    color=[cl_key, "sample_id"],
                    title=title,
                    show=False,
                )
                plt.savefig(out_png, dpi=150, bbox_inches="tight")
                plt.close()

    # mean.var.plot_disp_* DRs
    for disp in disps:
        dr = f"pca_mean.var.plot_disp_{disp}"
        _do_one_dr(dr)

    # vst_top_* DRs
    for n in n_features:
        dr = f"pca_vst_top_{n}"
        _do_one_dr(dr)

    # restore original X_umap if it existed
    if had_umap:
        adata.obsm["X_umap"] = umap_backup

    return adata


disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
neighbors = [10, 20, 30, 50]
res = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

adata = plot_umaps_for_grid(
    adata,
    disps=disps,
    n_features=n_features,
    neighbors=neighbors,
    res=res,
    UMAP_PATH=UMAP_PATH,
    random_state=0,
)


adata.write(os.path.join(UMAP_PATH, "adata_umap.h5ad"))
