import os

import anndata as ad
import pandas as pd
import scanpy as sc
from utils.scanpy_clustering_utils import run_clustering, run_clustering_analysis

# set pandas string handling to use builtin str type, not pyarrow to avoid anndata IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

# Load existing AnnData object
SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy"
adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_reduced.h5ad"))

# Create output directory
CLUSTERING_OUTPUT_PATH = os.path.join(SCANPY_OBJECT_PATH, "clustering_stats")
os.makedirs(CLUSTERING_OUTPUT_PATH, exist_ok=True)


# Parameters for variable gene selection
disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
xmin = 0.1
xmax = 10
neighbors = [10, 20, 30, 50]
res = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]


def main() -> None:
    adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_reduced.h5ad"))

    adata = run_clustering(
        adata,
        disps=disps,
        n_features=n_features,
        neighbors=neighbors,
        res=res,
        CLUSTERING_OUTPUT_PATH=CLUSTERING_OUTPUT_PATH,
    )

    run_clustering_analysis(
        adata,
        CLUSTERING_OUTPUT_PATH,
        disps,
        n_features,
        neighbors,
        res,
    )

    adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_merged_clustered.h5ad"))


if __name__ == "__main__":
    main()
