#!/usr/bin/env python3
"""Reassign labels for cluster-celltypes which have been manually curated."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

SCVI_PATH = Path("project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output")
INPUT_FILE = SCVI_PATH / "query_concat_curated_clustered.h5ad"

candidate_clusters = {
    "B cells": ["8", "20", "21"],
    "ADAMDEC1+ stromal": ["17", "19"],
    "Macrophages": ["34", "14", "2", "3"],
    "Plasma cells": ["4", "16"],
}


def reassign_candidate_clusters_by_centroid(
    adata,
    candidate_clusters,
    label_key="label_spreading_prediction_filtered",
    cluster_key="leiden_n30_r1p0",
    embedding_key="X_embeddings",
    output_key="label_spreading_prediction_filtered",
    distance_metric="cosine",
):
    """
    Split selected cell-type labels by specified Leiden clusters, then assign
    remaining cells from each selected cell type to the nearest candidate
    cluster centroid.
    """

    source_labels = adata.obs[label_key].astype(str)
    leiden = adata.obs[cluster_key].astype(str)
    X = np.asarray(adata.obsm[embedding_key])

    # Start from the existing labels
    new_labels = adata.obs[label_key].astype(str).copy()

    reassignment_records = []

    for cell_type, clusters in candidate_clusters.items():
        clusters = [str(cluster) for cluster in clusters]

        celltype_mask = source_labels.eq(cell_type).to_numpy()

        if celltype_mask.sum() == 0:
            print(f"No cells found for {cell_type!r}; skipping.")
            continue

        # Build centroids from the candidate clusters
        centroids = {}
        centroid_labels = {}

        for cluster in clusters:
            anchor_mask = (source_labels.eq(cell_type) & leiden.eq(cluster)).to_numpy()

            if anchor_mask.sum() == 0:
                print(
                    f"No anchor cells for {cell_type!r}, "
                    f"{cluster_key}={cluster!r}; skipping this centroid.",
                )
                continue

            new_label = f"{cell_type}_cl_{cluster}"

            centroids[cluster] = X[anchor_mask].mean(axis=0)
            centroid_labels[cluster] = new_label

            # Directly relabel cells already in the candidate cluster
            new_labels.loc[anchor_mask] = new_label

            reassignment_records.append(
                {
                    "cell_type": cell_type,
                    "source": "in_candidate_cluster",
                    "assigned_label": new_label,
                    cluster_key: cluster,
                    "n_cells": int(anchor_mask.sum()),
                },
            )

        if len(centroids) == 0:
            print(f"No usable centroids for {cell_type!r}; skipping.")
            continue

        # Assign remaining cells to nearest centroid
        remaining_mask = (
            source_labels.eq(cell_type) & ~leiden.isin(clusters)
        ).to_numpy()

        if remaining_mask.sum() == 0:
            continue

        centroid_cluster_ids = list(centroids.keys())
        centroid_matrix = np.vstack(
            [centroids[cluster] for cluster in centroid_cluster_ids],
        )

        distances = cdist(
            X[remaining_mask],
            centroid_matrix,
            metric=distance_metric,
        )

        nearest_idx = distances.argmin(axis=1)
        nearest_clusters = [centroid_cluster_ids[i] for i in nearest_idx]

        remaining_obs_names = adata.obs_names[remaining_mask]

        for cluster in centroid_cluster_ids:
            assigned_mask = np.array(nearest_clusters) == cluster
            assigned_obs = remaining_obs_names[assigned_mask]

            assigned_label = centroid_labels[cluster]
            new_labels.loc[assigned_obs] = assigned_label

            reassignment_records.append(
                {
                    "cell_type": cell_type,
                    "source": "nearest_centroid",
                    "assigned_label": assigned_label,
                    cluster_key: cluster,
                    "n_cells": len(assigned_obs),
                },
            )

    adata.obs[output_key] = pd.Categorical(new_labels)

    reassignment_summary = pd.DataFrame(reassignment_records)

    return adata, reassignment_summary


def main() -> None:
    """Reassign labels."""
    adata = ad.read_h5ad(INPUT_FILE)

    adata, reassignment_summary = reassign_candidate_clusters_by_centroid(
        adata,
        candidate_clusters,
    )

    reassignment_summary.to_csv(SCVI_PATH / "reassignment_summary.csv", index=False)

    adata.write(SCVI_PATH / "query_concat_curated_clustered_reassigned.h5ad")

    print(f"written to {SCVI_PATH}")


if __name__ == "__main__":
    main()
