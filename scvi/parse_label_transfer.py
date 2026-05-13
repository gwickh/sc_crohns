#!/usr/bin/env python3
"""Semi-supervised learning of Parse labels."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import jensenshannon
from sklearn.preprocessing import LabelEncoder
from sklearn.semi_supervised import LabelSpreading

# set pandas string handling to use builtin str type, not pyarrow to avoid IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

TUNING_DIR = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/sysvi_tuning",
)

map_curated = {
    "Activated CD4 T": "CD4 T cell",
    "Activated CD8 T": "CD8 T cell",
    "Adult Glia": "Neuronal",
    "arterial capillary": "Endothelial",
    "BEST2+ Goblet cell": "Goblet",
    "BEST4+ epithelial": "BEST4+ epithelial",
    "CD8 Tmem": "CD8 T mem",
    "cDC1": "cDC",
    "cDC2": "cDC",
    "Colonocyte": "Enterocyte",
    "Contractile pericyte (PLN+)": "Pericyte",
    "CX3CR1+ CD8 Tmem": "CD8 T mem",
    "Cycling B cell": "B cells",
    "Cycling plasma cell": "Plasma cells",
    "D cells (SST+)": "Enteroendocrine cell",
    "DZ GC cell": "B cells",
    "EC cells (TAC1+)": "Enteroendocrine cell",
    "EECs": "Enteroendocrine cell",
    "Enterocyte": "Enterocyte",
    "FCRL4+ Memory B": "B cells",
    "FDC": "ADAMDEC1+ stromal",
    "GC B cell": "B cells",
    "gdT": "NK T cell",
    "Goblet cell": "Goblet",
    "IgA plasma cell": "Plasma cells",
    "IgG plasma cell": "Plasma cells",
    "IgM plasma cell": "Plasma cells",
    "ILC3": "NK cell",
    "Immature B": "B cells",
    "L cells (PYY+)": "Enteroendocrine cell",
    "LEC1 (ACKR4+)": "Endothelial",
    "LEC3 (ADGRG3+)": "Endothelial",
    "LEC5 (CLDN11+)": "Endothelial",
    "LEC6 (ADAMTS4+)": "Endothelial",
    "Lymphoid DC": "cDC",
    "LYVE1+ Macrophage": "Macrophages",
    "LZ GC cell": "B cells",
    "Macrophages": "Macrophages",
    "MAIT cell": "CD8 T cell",
    "Mast cell": "Mast cell",
    "Mature arterial EC": "Endothelial",
    "Mature venous EC": "Endothelial",
    "Memory B": "B cells",
    "Mesothelium (PRG4+)": "Mesenchymal",
    "Microfold cell": "Enterocyte",
    "mLN Stroma (FMO2+)": "Mesenchymal",
    "mLTo": "Mesenchymal",
    "MMP9+ Inflammatory macrophage": "Macrophages",
    "Monocytes": "Monocytes",
    "myofibroblast": "Myofibroblast",
    "N cells (NTS+)": "Enteroendocrine cell",
    "Naive B": "B cells",
    "NK cell": "NK cell",
    "NK T cell": "NK T cell",
    "Paneth": "Paneth",
    "pDC": "cDC",
    "Pericyte": "Pericyte",
    "RBC": "Red blood cells",
    "SELL+ CD4 T": "CD4 T cell",
    "SELL+ CD8 T": "Naïve CD8 T cell",
    "STAT1+ Naive B": "B cells",
    "Stem cells": "Stem cell",
    "Stromal 1 (ADAMDEC1+)": "ADAMDEC1+ stromal",
    "Stromal 1 (CCL11+)": "ADAMDEC1+ stromal",
    "Stromal 2 (NPY+)": "NPY+ Stromal",
    "Stromal 3 (C7+)": "ADAMDEC1+ stromal",
    "Stromal 4 (MMP1+)": "ADAMDEC1+ stromal",
    "T reticular": "ADAMDEC1+ stromal",
    "TA": "Stem cell",
    "Tfh": "Tfh",
    "Th1": "CD4 T cell",
    "Th17": "CD4 T cell",
    "Transitional Stromal 3 (C3+)": "Transitional Stromal C3+",
    "Treg": "CD4 T cell",
    "TRGV2 gdT": "NK T cell",
    "Tuft": "Tuft",
    "Fetal venous EC": "Endothelial",
    "TRGV4 gdT": "NK T cell",
    "TRGV5/7 gdT": "NK T cell",
    "myofibroblast (RSPO2+)": "Myofibroblast",
    "LEC4 (STAB2+)": "Endothelial",
    "LEC2 (MADCAM1+)": "Endothelial",
    "Progenitor (NEUROG3+)": "Enteroendocrine cell",
}


def preprocess_adata(
    adata: sc.AnnData,
    embedding_key: str = "X_embeddings",
    platform_key: str = "platform",
    label_key: str = "Integrated_05",
    reference_platform: str = "10X_Chromium",
    query_platform: str = "Parse",
    curated_label_key: str = "curated_labels",
) -> None:
    """Prepare AnnData for diagnosis-stratified label spreading."""
    if embedding_key not in adata.obsm:
        msg = f"{embedding_key!r} not found in adata.obsm."

        raise KeyError(msg)

    required_obs = [platform_key, label_key, "diagnosis"]

    missing = [c for c in required_obs if c not in adata.obs.columns]

    if missing:
        msg = f"Missing columns in adata.obs: {missing}"

        raise KeyError(msg)

    ref_mask = adata.obs[platform_key].eq(reference_platform)

    query_mask = adata.obs[platform_key].eq(query_platform)

    if ref_mask.sum() == 0:
        msg = f"No reference cells found for {platform_key} == {reference_platform!r}"

        raise ValueError(msg)

    if query_mask.sum() == 0:
        msg = f"No query cells found for {platform_key} == {query_platform!r}"

        raise ValueError(msg)

    if adata.obs.loc[ref_mask, label_key].isna().any():
        msg = (
            f"Missing values found in adata.obs[{label_key!r}] "
            f"among {reference_platform!r} reference cells."
        )

        raise ValueError(msg)

    adata.obs[curated_label_key] = (
        adata.obs[label_key].astype("object").map(map_curated)
    )

    # Keep original labels if not found in map_curated

    adata.obs[curated_label_key] = adata.obs[curated_label_key].fillna(
        adata.obs[label_key].astype("object"),
    )

    adata.obs[curated_label_key] = adata.obs[curated_label_key].astype("category")

    print("\nCurated label counts:")
    print(adata.obs.loc[ref_mask, curated_label_key].value_counts(dropna=False))

    print("\nDiagnosis for platform:")
    print(pd.crosstab(adata.obs["diagnosis"], adata.obs[platform_key]))


def train_parse_label_transfer(
    adata: sc.AnnData,
    embedding_key: str = "X_embeddings",
    diagnosis_key: str = "diagnosis",
    kernel: str = "knn",
    n_neighbors: int = 30,
    alpha: float = 0.2,
    confidence_threshold: float = 0.8,
    max_iter: int = 1000,
    tol: float = 1e-3,
    platform_key: str = "platform",
    query_platform: str = "Parse",
    reference_platform: str = "10X_Chromium",
    reference_label_key: str = "curated_labels",
    output_label_key: str = "label_spreading_prediction",
    output_confidence_key: str = "label_spreading_confidence",
    output_unknown_key: str = "label_spreading_prediction_filtered",
) -> pd.DataFrame:
    """Run label spreading separately within each diagnosis and concatenate summaries."""
    if embedding_key not in adata.obsm:
        msg = f"{embedding_key!r} not found in adata.obsm."

        raise KeyError(msg)

    required_obs = [diagnosis_key, platform_key, reference_label_key]

    missing = [c for c in required_obs if c not in adata.obs.columns]

    if missing:
        msg = f"Missing columns in adata.obs: {missing}"

        raise KeyError(msg)

    X_all = np.asarray(adata.obsm[embedding_key])

    if X_all.shape[0] != adata.n_obs:
        msg = f"Embedding has {X_all.shape[0]} rows but adata has {adata.n_obs} cells."

        raise ValueError(msg)

    # Initialise output columns for the full object
    adata.obs[output_label_key] = pd.Series(index=adata.obs_names, dtype="object")
    adata.obs[output_confidence_key] = np.nan
    adata.obs[output_unknown_key] = pd.Series(index=adata.obs_names, dtype="object")

    summary_dfs = []

    diagnosis_values = sorted(adata.obs[diagnosis_key].dropna().astype(str).unique())

    for diagnosis in diagnosis_values:
        print(f"\n{'=' * 80}")
        print(f"Running label spreading for {diagnosis_key} = {diagnosis!r}")
        print(f"{'=' * 80}")
        diagnosis_mask = adata.obs[diagnosis_key].astype(str).eq(diagnosis)

        ref_mask = (
            diagnosis_mask
            & adata.obs[platform_key].eq(reference_platform)
            & adata.obs[reference_label_key].notna()
        )

        query_mask = diagnosis_mask & adata.obs[platform_key].eq(query_platform)
        n_ref = int(ref_mask.sum())
        n_query = int(query_mask.sum())

        print(f"Reference cells: {n_ref:,}")
        print(f"Parse query cells: {n_query:,}")

        if n_ref == 0:
            print(f"[skip] No labelled reference cells for diagnosis {diagnosis!r}.")
            continue

        if n_query == 0:
            print(f"[skip] No Parse query cells for diagnosis {diagnosis!r}.")
            continue

        print("\nReference label counts:")

        print(adata.obs.loc[ref_mask, reference_label_key].value_counts())

        # Subset contains only reference + query cells from this diagnosis
        subset_mask = ref_mask | query_mask
        subset_indices = np.where(subset_mask.to_numpy())[0]
        subset_obs_names = adata.obs_names[subset_indices]
        X = X_all[subset_indices]

        # Encode labels for reference cells in this diagnosis
        le = LabelEncoder()
        ref_labels = adata.obs.loc[ref_mask, reference_label_key].astype(str)
        le.fit(ref_labels)
        y = np.full(subset_indices.shape[0], -1, dtype=int)

        # Positions of reference cells within the subset
        ref_positions = np.where(ref_mask.loc[subset_mask].to_numpy())[0]
        y[ref_positions] = le.transform(ref_labels)

        model = LabelSpreading(
            kernel=kernel,
            n_neighbors=n_neighbors,
            alpha=alpha,
            max_iter=max_iter,
            tol=tol,
            n_jobs=-1,
        )

        model.fit(X, y)
        pred_encoded = model.transduction_
        pred_labels = le.inverse_transform(pred_encoded)
        confidence = model.label_distributions_.max(axis=1)

        # Write predictions back to full adata
        adata.obs.loc[subset_obs_names, output_label_key] = pred_labels
        adata.obs.loc[subset_obs_names, output_confidence_key] = confidence
        filtered = pd.Series(pred_labels, index=subset_obs_names, dtype="object")
        query_subset_mask = adata.obs.loc[subset_obs_names, platform_key].eq(
            query_platform,
        )

        low_conf_parse = query_subset_mask & (
            adata.obs.loc[subset_obs_names, output_confidence_key]
            < confidence_threshold
        )

        filtered.loc[low_conf_parse] = "Unknown"

        adata.obs.loc[subset_obs_names, output_unknown_key] = filtered

        n_low_conf_parse = int(low_conf_parse.sum())

        pct_low_conf_parse = 100 * n_low_conf_parse / n_query if n_query > 0 else np.nan

        print(
            f"{n_low_conf_parse:,}/{n_query:,} Parse cells "
            f"({pct_low_conf_parse:.2f}%) fall below "
            f"confidence threshold {confidence_threshold}.",
        )

        # Summary for diagnosis
        ref_counts = (
            adata.obs.loc[ref_mask, reference_label_key]
            .astype("object")
            .fillna("Unknown")
            .value_counts()
            .rename("n_10X")
        )

        parse_counts = (
            adata.obs.loc[query_mask, output_unknown_key]
            .astype("object")
            .fillna("Unknown")
            .value_counts()
            .rename("n_Parse_predicted")
        )

        summary = pd.concat([ref_counts, parse_counts], axis=1).fillna(0)
        summary["n_10X"] = summary["n_10X"].astype(int)
        summary["n_Parse_predicted"] = summary["n_Parse_predicted"].astype(int)

        summary["pct_10X"] = (
            100 * summary["n_10X"] / (summary["n_10X"].sum() - summary["n_10X"].iloc[0])
            if summary["n_10X"].sum() > 0
            else np.nan
        )

        summary["pct_Parse_predicted"] = (
            100
            * summary["n_Parse_predicted"]
            / (
                summary["n_Parse_predicted"].sum()
                - summary["n_Parse_predicted"].iloc[0]
            )
            if summary["n_Parse_predicted"].sum() > 0
            else np.nan
        )

        summary = (
            summary.rename_axis("cell_type")
            .reset_index()
            .sort_values("n_Parse_predicted", ascending=False)
        )

        summary.insert(0, diagnosis_key, diagnosis)
        summary["n_reference_cells"] = n_ref
        summary["n_parse_cells"] = n_query
        summary["n_low_conf_parse"] = n_low_conf_parse
        summary["pct_low_conf_parse"] = pct_low_conf_parse
        summary["kernel"] = kernel
        summary["n_neighbors"] = n_neighbors
        summary["alpha"] = alpha
        summary["confidence_threshold"] = confidence_threshold
        summary_dfs.append(summary)

    if not summary_dfs:
        raise ValueError("No diagnosis-specific label spreading runs were completed.")

    combined_summary = pd.concat(summary_dfs, ignore_index=True)

    # Clean remaining missing predictions, if any diagnosis was skipped

    adata.obs[output_label_key] = (
        adata.obs[output_label_key].fillna("Unknown").astype("category")
    )

    adata.obs[output_unknown_key] = (
        adata.obs[output_unknown_key].fillna("Unknown").astype("category")
    )

    return combined_summary


def compute_js_divergence(
    summary: pd.DataFrame,
    exclude_labels: tuple[str, ...] = ("Unknown",),
    pseudocount: float = 0.0,
) -> float:
    """
    Compute Jensen-Shannon distance between 10X and Parse cell-type counts.
    """

    summary["unique_cell_type"] = summary["diagnosis"] + "_" + summary["cell_type"]
    counts_10x = summary.set_index("unique_cell_type")["n_10X"]
    counts_parse = summary.set_index("unique_cell_type")["n_Parse_predicted"]

    # align cell types
    all_labels = sorted(set(counts_10x.index) | set(counts_parse.index))
    p = counts_10x.reindex(all_labels, fill_value=0.0)
    q = counts_parse.reindex(all_labels, fill_value=0.0)

    # remove Unknown / low-confidence labels
    if exclude_labels:
        keep = ~p.index.astype(str).str.contains(
            "|".join(map(str, exclude_labels)),
            regex=True,
        )
        p = p.loc[keep]
        q = q.loc[keep]

    # pseudocount smoothing
    if pseudocount > 0:
        p = p + pseudocount
        q = q + pseudocount

    if p.sum() == 0 or q.sum() == 0:
        return np.nan

    p = p / p.sum()
    q = q / q.sum()

    return float(jensenshannon(p, q, base=2))


def main() -> None:
    """Run label transfer."""
    adata = sc.read_h5ad(TUNING_DIR / "c561826c_sysvi.h5ad")

    adata.obs["diagnosis"] = np.where(
        adata.obs["sample_id"].str.contains("crohns", case=False, na=False),
        "Crohn's disease",
        "Normal",
    )

    param_grid = {
        "n_neighbors": [5, 10, 15, 20],
        "alpha": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    }

    preprocess_adata(adata)

    js_d_list = []

    for n_neighbors in param_grid["n_neighbors"]:
        for alpha in param_grid["alpha"]:
            summary = train_parse_label_transfer(
                adata,
                kernel="knn",
                n_neighbors=n_neighbors,
                alpha=alpha,
            )

            summary.to_csv(
                TUNING_DIR
                / f"label_spreading_diagnosis_knn_alpha_{alpha}_n_{n_neighbors}.csv",
                index=False,
            )

            js_d = compute_js_divergence(summary)

            js_d_list.append((n_neighbors, alpha, js_d))

    df = pd.DataFrame(js_d_list, columns=["n_neighbors", "alpha", "js_distance"])
    df.to_csv(TUNING_DIR / "label_spreading_diagnosis_js_distance.csv", index=False)


if __name__ == "__main__":
    main()
