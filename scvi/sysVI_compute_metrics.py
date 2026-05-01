#!/usr/bin/env python3
"""Compute integration and biological preservation metrics for sysVI."""

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation

# set pandas string handling to use builtin str type, not pyarrow to avoid IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

ref_dir = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/"
tuning_dir = f"{ref_dir}/sysvi_tuning"


def add_cell_type_annotation(
    query_adata_file: Path,
    cell_type_key: str = "Integrated_05",
) -> sc.AnnData:
    """Add cell type annotation from reference data to query data."""
    query_adata = sc.read_h5ad(query_adata_file)
    print(query_adata.obs["platform"].unique())

    if cell_type_key in query_adata.obs.columns:
        print(f"{cell_type_key} already exists in query data, skipping annotation.")
        return

    print(f"Adding {cell_type_key} annotation to query data.")

    cell_type_ref_adata = sc.read_h5ad(Path(ref_dir) / "query_concat.h5ad")

    if cell_type_key not in cell_type_ref_adata.obs.columns:
        msg = f"Reference AnnData object must have '{cell_type_key}' column in .obs"
        raise ValueError(msg)

    # ensure gene names are unique
    if (
        not cell_type_ref_adata.obs_names.is_unique
        or not query_adata.obs_names.is_unique
    ):
        msg = "Observation names in either reference or query data are not unique."
        raise ValueError(msg)

    missing_ref_obs = query_adata[
        query_adata.obs["platform"] == "10X_Genomics"
    ].obs_names.difference(
        cell_type_ref_adata[
            cell_type_ref_adata.obs["platform"] == "10X_Chromium"
        ].obs_names,
    )

    if len(missing_ref_obs) > 0:
        msg = f"""{len(missing_ref_obs)} query cells are missing from reference.\n
        query: {query_adata.obs_names[:10].tolist()}\n
        ref: {cell_type_ref_adata.obs_names[:10].tolist()}
        """
        pd.DataFrame(missing_ref_obs).to_csv(
            Path(tuning_dir) / "missing_from_reference.csv",
            index=False,
        )
        raise ValueError(msg)

    query_adata.obs[cell_type_key] = (
        cell_type_ref_adata.obs[cell_type_key].reindex(query_adata.obs_names).to_numpy()
    )

    query_adata.write_h5ad(query_adata_file)


def compute_integration_metrics(
    query_adata_file: Path,
) -> pd.DataFrame:
    """Compute integration and biological preservation metrics for sysVI."""
    query_adata = sc.read_h5ad(query_adata_file)

    query_adata.obs["dummy_cell_type"] = "all_cells"

    bm = Benchmarker(
        query_adata,
        batch_key="platform",
        label_key="dummy_cell_type",
        batch_correction_metrics=BatchCorrection(
            bras=False,
            ilisi_knn=True,
            kbet_per_label=False,
            graph_connectivity=False,
            pcr_comparison=True,
        ),
        bio_conservation_metrics=BioConservation(
            isolated_labels=False,
            nmi_ari_cluster_labels_leiden=False,
            nmi_ari_cluster_labels_kmeans=False,
            silhouette_label=False,
            clisi_knn=False,
        ),
        embedding_obsm_keys=["X_embeddings"],
        n_jobs=6,
    )
    bm.benchmark()
    metrics = bm._results.copy()

    print(pd.DataFrame(metrics).transpose())

    return pd.DataFrame(metrics).transpose()


def compute_bio_metrics(
    query_adata_file: Path,
) -> pd.DataFrame:
    """Compute integration and biological preservation metrics for sysVI."""
    query_adata = sc.read_h5ad(query_adata_file)

    query_adata = query_adata[query_adata.obs["platform"] == "10X_Chromium"].copy()

    bm = Benchmarker(
        query_adata,
        batch_key="platform",
        label_key="Integrated_05",
        batch_correction_metrics=BatchCorrection(
            bras=False,
            ilisi_knn=False,
            kbet_per_label=False,
            graph_connectivity=False,
            pcr_comparison=False,
        ),
        bio_conservation_metrics=BioConservation(
            isolated_labels=True,
            nmi_ari_cluster_labels_leiden=False,
            nmi_ari_cluster_labels_kmeans=True,
            silhouette_label=True,
            clisi_knn=True,
        ),
        embedding_obsm_keys=["X_embeddings"],
        n_jobs=6,
    )
    bm.benchmark()
    metrics = bm._results.copy()

    print(list(pd.DataFrame(metrics).transpose().columns.values))

    return pd.DataFrame(metrics).transpose()


def main() -> None:
    """Compute metrics for sysVI."""
    adatas = list(Path(tuning_dir).glob("*_sysvi.h5ad"))
    print(f"Found {len(adatas)} sysVI files.")

    metrics_df = pd.DataFrame(
        columns=[
            "params",
            "ilisi_knn",
            "pcr_comparison",
            "isolated_labels",
            "nmi_ari_cluster_labels_kmeans_nmi",
            "nmi_ari_cluster_labels_kmeans_ari",
            "silhouette_label",
            "clisi_knn",
            "batch_correction_score",
            "bio_conservation_score",
            "aggregate_score",
        ]
    )
    for file in adatas:
        try:
            add_cell_type_annotation(
                query_adata_file=file,
            )
            integration_metrics = compute_integration_metrics(file)

            bio_metrics = compute_bio_metrics(file)

            params = str(file).split("/")[-1].split("_sysvi.h5ad")[0]
            metrics_df.loc[len(metrics_df), "params"] = params

            metrics_df.loc[len(metrics_df) - 1, "ilisi_knn"] = integration_metrics.loc[
                "X_embeddings", "ilisi_knn"
            ]

            metrics_df.loc[len(metrics_df) - 1, "pcr_comparison"] = (
                integration_metrics.loc["X_embeddings", "pcr_comparison"]
            )

            metrics_df.loc[len(metrics_df) - 1, "isolated_labels"] = bio_metrics.loc[
                "X_embeddings", "isolated_labels"
            ]

            metrics_df.loc[len(metrics_df) - 1, "nmi_ari_cluster_labels_kmeans_nmi"] = (
                bio_metrics.loc["X_embeddings", "nmi_ari_cluster_labels_kmeans_nmi"]
            )
            metrics_df.loc[len(metrics_df) - 1, "nmi_ari_cluster_labels_kmeans_ari"] = (
                bio_metrics.loc["X_embeddings", "nmi_ari_cluster_labels_kmeans_ari"]
            )
            metrics_df.loc[len(metrics_df) - 1, "silhouette_label"] = bio_metrics.loc[
                "X_embeddings", "silhouette_label"
            ]

            metrics_df.loc[len(metrics_df) - 1, "batch_correction_score"] = (
                metrics_df.loc[len(metrics_df) - 1, ["ilisi_knn", "pcr_comparison"]]
                .astype(float)
                .mean()
            )

            metrics_df.loc[len(metrics_df) - 1, "bio_conservation_score"] = (
                metrics_df.loc[
                    len(metrics_df) - 1,
                    [
                        "isolated_labels",
                        "nmi_ari_cluster_labels_kmeans_nmi",
                        "nmi_ari_cluster_labels_kmeans_ari",
                        "silhouette_label",
                        "clisi_knn",
                    ],
                ]
                .astype(float)
                .mean()
            )

            metrics_df.loc[len(metrics_df) - 1, "aggregate_score"] = 0.4 * float(
                metrics_df.loc[len(metrics_df) - 1, "batch_correction_score"]
            ) + 0.6 * float(
                metrics_df.loc[len(metrics_df) - 1, "bio_conservation_score"]
            )

        except OSError as err:
            print(f"[skip] Could not read {file}: {err}")
        continue

    metrics_df.to_csv(Path(tuning_dir) / "metrics.csv", index=False)


if __name__ == "__main__":
    main()
