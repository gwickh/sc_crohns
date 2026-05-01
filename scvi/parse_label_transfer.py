#!/usr/bin/env python3
"""Semi-supervised learning of Parse labels."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import LabelEncoder
from sklearn.semi_supervised import LabelSpreading

# set pandas string handling to use builtin str type, not pyarrow to avoid IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

TUNING_DIR = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/old_sysvi_tuning"
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
):
    """Annotate Parse data by label spreading on the latent space."""
    if embedding_key not in adata.obsm:
        msg = f"{embedding_key!r} not found in adata.obsm."
        raise KeyError(msg)

    required_obs = [platform_key, label_key]
    missing = [c for c in required_obs if c not in adata.obs.columns]
    if missing:
        msg = f"Missing columns in adata.obs: {missing}"
        raise KeyError(msg)

    ref_mask = adata.obs[platform_key].eq(reference_platform)
    query_mask = adata.obs[platform_key].eq(query_platform)

    if query_mask.sum() == 0:
        msg = f"No query cells found for {platform_key} == {query_platform!r}"
        raise ValueError(msg)

    if adata.obs.loc[ref_mask, label_key].isna().any():
        msg = (
            f"Missing values found in adata.obs[{label_key!r}] "
            f"among {reference_platform!r} reference cells."
        )
        raise ValueError(msg)

    X = np.asarray(adata.obsm[embedding_key])

    adata.obs[curated_label_key] = (
        adata.obs[label_key].astype("object").map(map_curated)
    )

    adata.obs[curated_label_key] = adata.obs[curated_label_key].astype("category")

    le = LabelEncoder()
    ref_labels = adata.obs.loc[ref_mask, curated_label_key].astype(str)
    encoded_ref_labels = le.fit_transform(ref_labels)

    y = np.full(adata.n_obs, -1, dtype=int)

    y[ref_mask.to_numpy()] = encoded_ref_labels

    print("\nClasses:")
    for i, cls in enumerate(le.classes_):
        print(f"{i}: {cls}")

    return X, y, le


def train_parse_label_transfer(
    adata: sc.AnnData,
    X: np.ndarray,
    y: np.ndarray,
    le: LabelEncoder,
    kernel: str = "knn",
    n_neighbors: int = 30,
    alpha: float = 0.2,
    confidence_threshold: float = 0.8,
    max_iter: int = 1000,
    tol: float = 1e-3,
    platform_key: str = "platform",
    query_platform: str = "Parse",
    output_label_key: str = "label_spreading_prediction",
    output_confidence_key: str = "label_spreading_confidence",
    output_unknown_key: str = "label_spreading_prediction_filtered",
    reference_platform: str = "10X_Chromium",
    reference_label_key: str = "curated_labels",
) -> None:
    """Train label spreading model and save results."""

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

    adata.obs[output_label_key] = pd.Categorical(pred_labels)

    adata.obs[output_confidence_key] = confidence

    filtered = pd.Series(pred_labels, index=adata.obs_names, dtype="object")

    query_mask = adata.obs[platform_key].eq(query_platform)

    low_conf_parse = query_mask & (
        adata.obs[output_confidence_key] < confidence_threshold
    )

    filtered.loc[low_conf_parse] = "Unknown"

    adata.obs[output_unknown_key] = pd.Categorical(filtered)

    n_parse = int(query_mask.sum())
    n_low_conf_parse = int(low_conf_parse.sum())
    pct_low_conf_parse = 100 * n_low_conf_parse / n_parse if n_parse > 0 else np.nan

    print(
        f"{n_low_conf_parse:,}/{n_parse:,} Parse cells "
        f"({pct_low_conf_parse:.2f}%) fall below "
        f"confidence threshold {confidence_threshold}."
    )

    ref_mask = adata.obs[platform_key].eq(reference_platform)

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
        100 * summary["n_10X"] / summary["n_10X"].sum()
        if summary["n_10X"].sum() > 0
        else np.nan
    )

    summary["pct_Parse_predicted"] = (
        100 * summary["n_Parse_predicted"] / summary["n_Parse_predicted"].sum()
        if summary["n_Parse_predicted"].sum() > 0
        else np.nan
    )

    summary = (
        summary.rename_axis("cell_type")
        .reset_index()
        .sort_values("n_Parse_predicted", ascending=False)
    )

    summary.to_csv(
        TUNING_DIR
        / f"label_spreading_cell_type_summary_{kernel}_alpha_{alpha}_n_{n_neighbors}.csv",
        index=False,
    )

    print(f"""
            \nLabel spreading summary:
            kernel: {kernel}
            n_neighbors: {n_neighbors}
            alpha: {alpha}
            pct failed Parse predictions: {pct_low_conf_parse:.2f}%
          """)
    # adata.write_h5ad(TUNING_DIR / "label_spreading_predictions.h5ad")


def main() -> None:
    """Run label transfer."""
    adata = sc.read_h5ad(TUNING_DIR / "760bb93b_sysvi.h5ad")

    param_grid = {
        "n_neighbors": [5, 10, 15, 20],
        "alpha": [0.01, 0.025, 0.05, 0.1, 0.2],
    }

    X, y, le = preprocess_adata(adata)

    for n_neighbors in param_grid["n_neighbors"]:
        for alpha in param_grid["alpha"]:
            train_parse_label_transfer(
                adata,
                X,
                y,
                le,
                kernel="rbf",
                n_neighbors=n_neighbors,
                alpha=alpha,
            )


if __name__ == "__main__":
    main()
