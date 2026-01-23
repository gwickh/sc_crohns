import os

import anndata as an
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

print("SCRIPT STARTED", flush=True)

SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy"

adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_clustered.h5ad"))

UMAP_PATH = os.path.join(SCANPY_OBJECT_PATH, "UMAP_plots")
os.makedirs(UMAP_PATH) if not os.path.exists(UMAP_PATH) else None


sid = adata.obs["sample_id"].astype(str)
status = np.where(
    sid.str.contains("crohns", case=False, na=False),
    "crohns",
    np.where(sid.str.contains("normal", case=False, na=False), "normal", "other"),
)
adata.obs["crohns_or_normal"] = pd.Categorical(
    status, categories=["normal", "crohns", "other"]
)


def plot_umaps_for_grid(
    adata,
    disps,
    n_features,
    neighbors,
    res,
    UMAP_PATH,
    random_state=0,
):
    had_umap = "X_umap" in adata.obsm
    umap_backup = adata.obsm["X_umap"].copy() if had_umap else None

    def _do_one_dr(dr) -> None:
        dr_out = os.path.join(UMAP_PATH, dr)
        os.makedirs(dr_out, exist_ok=True)

        for k in neighbors:
            neighbors_key = f"neighbors_{dr}_k{k}"
            if neighbors_key not in adata.uns:
                print(f"[skip] Missing neighbors: {neighbors_key}")
                continue

            basis = f"umap_{dr}_k{k}"
            umap_key = f"X_{basis}"

            if umap_key not in adata.obsm:
                sc.tl.umap(
                    adata, neighbors_key=neighbors_key, random_state=random_state
                )
                adata.obsm[umap_key] = adata.obsm["X_umap"].copy()

            for r in res:
                cl_key = f"clusters_{dr}_k{k}_r{r}"
                if cl_key not in adata.obs:
                    print(f"[skip] Missing clusters: {cl_key}")
                    continue

                out_png = os.path.join(dr_out, f"{basis}__r{r}.png")
                base_title = f"{dr} | k={k} | r={r}"

                colors = [cl_key, "sample_id", "crohns_or_normal", "predicted_doublet"]
                titles = [
                    f"{base_title} | clusters",
                    f"{base_title} | sample_id",
                    f"{base_title} | crohns/normal",
                    f"{base_title} | predicted doublet",
                ]

                fig = sc.pl.embedding(
                    adata,
                    basis=basis,
                    color=colors,
                    title=titles,
                    ncols=4,
                    legend_loc="right margin",
                    legend_fontsize=7,
                    show=False,
                    return_fig=True,
                )

                fig.set_size_inches(7.0 * 4, 6.0)  # figure size for 4 panels
                fig.subplots_adjust(wspace=0.8)

                fig.savefig(out_png, dpi=150, bbox_inches="tight", pad_inches=0.4)
                plt.close(fig)

    for disp in disps:
        _do_one_dr(f"pca_mean.var.plot_disp_{disp}")

    for n in n_features:
        _do_one_dr(f"pca_vst_top_{n}")

    if had_umap:
        adata.obsm["X_umap"] = umap_backup

    return adata


disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
neighbors = [10, 20, 30, 50]
res = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]


s = adata.obs["predicted_doublet"]
print("unique:", pd.unique(s)[:10], flush=True)
print("value_counts:\n", s.value_counts(dropna=False), flush=True)

adata = plot_umaps_for_grid(
    adata,
    disps=disps,
    n_features=n_features,
    neighbors=neighbors,
    res=res,
    UMAP_PATH=UMAP_PATH,
    random_state=0,
)


adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_umap.h5ad"))
