import os

import anndata as an
import numpy as np
import scanpy as sc

# Load existing AnnData object
SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy"
adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_reduced.h5ad"))

# Create output directory
CLUSTERING_OUTPUT_PATH = os.path.join(SCANPY_OBJECT_PATH, "clustering_stats")
os.makedirs(CLUSTERING_OUTPUT_PATH) if not os.path.exists(
    CLUSTERING_OUTPUT_PATH
) else None


def clust_dr_scanpy(
    adata,
    dr,
    pcs=50,
    non_rand_sd_frac=0.5,
    K=(15, 30),
    res=(0.2, 0.5, 1.0),
    file="n_pcs.txt",
    random_state=0,
) -> an.AnnData:
    """
    Perform clustering on dimensionality reduced data
    """
    # Compute minimum std a PC must have to be considered “non-random” assuming
    #   last 10 PCs approximate random noise to exclude SVD artifacts
    stdevs = np.sqrt(np.asarray(adata.uns[f"{dr}_pca"]["variance"]))

    pcs = min(pcs, stdevs.shape[0])
    last10_idx = np.arange(max(0, pcs - 10), pcs)  # last 10 PCs (0-based)
    mean_stdev_last = stdevs[last10_idx].mean()
    min_stdev = (1.0 + non_rand_sd_frac) * mean_stdev_last

    # Loops through PCs to determine nPCs (number of informative PCs that
    #   explains > 0.05 variance with respect to the next)
    n_pcs = 0
    for i in range(pcs):
        if stdevs[i] < min_stdev:
            break
        n_pcs = i + 1

    with open(file, "w") as f:
        f.write(str(n_pcs) + "\n")

    # For each k build a NN graph using the top nPCs components
    for k in K:
        neighbors_key = f"neighbors_{dr}_k{k}"
        print(f"Computing neighbors for {dr} with k={k}...")
        sc.pp.neighbors(
            adata,
            n_neighbors=int(k),
            n_pcs=int(n_pcs),
            use_rep=f"X_{dr}",
            key_added=neighbors_key,
        )

        # For each resolution, run leiden clustering
        for r in res:
            cl_key = f"clusters_{dr}_k{k}_r{r}"
            print(f"Clustering {dr} with k={k}, r={r}...")
            sc.tl.leiden(
                adata,
                resolution=float(r),
                neighbors_key=neighbors_key,
                key_added=cl_key,
                random_state=random_state,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
    return adata


def run_clustering(
    adata,
    disps,
    n_features,
    neighbors,
    res,
    CLUSTERING_OUTPUT_PATH,
    pcs=50,
    non_rand_sd_frac=0.5,
    random_state=0,
) -> an.AnnData:
    """
    Run clustering on multiple PCA reductions
    """
    # mean.var.plot loop
    for disp in disps:
        filename = f"mean.var.plot_disp_{disp}"
        dr_name = f"pca_{filename}"
        out_file = os.path.join(CLUSTERING_OUTPUT_PATH, f"{dr_name}_num_PCs.txt")

        print(f"Clustering for {dr_name}...")
        clust_dr_scanpy(
            adata,
            dr=dr_name,
            pcs=pcs,
            non_rand_sd_frac=non_rand_sd_frac,
            K=neighbors,
            res=res,
            file=out_file,
            random_state=random_state,
        )

    # vst_top_* loop
    for n in n_features:
        filename = f"vst_top_{n}"
        dr_name = f"pca_{filename}"
        out_file = os.path.join(CLUSTERING_OUTPUT_PATH, f"{dr_name}_num_PCs.txt")

        print(f"Clustering for {dr_name}...")
        clust_dr_scanpy(
            adata,
            dr=dr_name,
            pcs=pcs,
            non_rand_sd_frac=non_rand_sd_frac,
            K=neighbors,
            res=res,
            file=out_file,
            random_state=random_state,
        )

    return adata


# Parameters for variable gene selection
disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
xmin = 0.1
xmax = 10
neighbors = [10, 20, 30, 50]
res = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

adata = run_clustering(
    adata,
    disps=disps,
    n_features=n_features,
    neighbors=neighbors,
    res=res,
    CLUSTERING_OUTPUT_PATH=CLUSTERING_OUTPUT_PATH,
)

adata.write(os.path.join(CLUSTERING_OUTPUT_PATH, "adata_merged_clustered.h5ad"))
