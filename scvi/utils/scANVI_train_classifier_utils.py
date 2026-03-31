#!/usr/bin/env python3
"""Utils for training scANVI classifier on reference and query objects."""

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

import scvi

scvi.settings.seed = 0
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

SCVI_PATH = Path("project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output")
SCVI_PATH.mkdir(parents=True, exist_ok=True)


# train model on reference object
def train_scanvi_ref(label: str = "Integrated_05") -> None:
    """Train scANVI model on reference object."""
    adata_ref = sc.read_h5ad(SCVI_PATH / "gca_ref_scvi.h5ad")

    # add ensembl_id column to var, using gene_ids or gene_id column if available
    if "gene_ids" in adata_ref.var.columns:
        adata_ref.var["ensembl_id"] = adata_ref.var["gene_ids"]
    elif "gene_id" in adata_ref.var.columns:
        adata_ref.var["ensembl_id"] = adata_ref.var["gene_id"]
    else:
        error_msg = (
            "Neither 'gene_ids' nor 'gene_id' column found in adata_ref.var",
            f"for sample {adata_ref.obs['sample_id'].iloc[0]}",
        )
        raise ValueError(error_msg)

    keep = adata_ref.var["ensembl_id"].notna()
    adata_ref = adata_ref[:, keep].copy()
    adata_ref.var_names = adata_ref.var["ensembl_id"].astype(str).to_numpy()

    vae_ref = scvi.model.SCVI.load(
        SCVI_PATH / "scvi_model",
        adata=adata_ref,
    )

    print(adata_ref.obs[label].values)
    adata_ref.obs["labels_scanvi"] = adata_ref.obs[label].astype(str)
    adata_ref.obs["labels_scanvi"] = adata_ref.obs["labels_scanvi"].astype("category")

    vae_ref_scanvi = scvi.model.SCANVI.from_scvi_model(
        vae_ref,
        unlabeled_category="Unknown",
        labels_key="labels_scanvi",
    )

    vae_ref_scanvi.train(
        max_epochs=100,
        n_samples_per_label=100,
        plan_kwargs={"classification_ratio": 10, "lr": 1e-3},
        early_stopping=True,
        early_stopping_patience=20,
        check_val_every_n_epoch=1,
    )
    vae_ref_scanvi.save(SCVI_PATH / "scanvi_model_ref", overwrite=True)

    adata_ref.obsm["X_embeddings_" + label] = vae_ref_scanvi.get_latent_representation()
    adata_ref.write_h5ad(SCVI_PATH / "gca_ref_scanvi.h5ad")

    labels = adata_ref.obs["labels_scanvi"]
    labels.to_csv(SCVI_PATH / "gca_ref_scanvi_labels.csv", sep="\t")


# train model on query object
def train_scanvi_query(
    adata: ad.AnnData,
    label: str = "Integrated_05",
) -> None:
    """Transfer scANVI model on query object and predict labels."""
    scvi.model.SCANVI.prepare_query_anndata(
        adata,
        SCVI_PATH / "scanvi_model_ref",
    )
    vae_q = scvi.model.SCANVI.load_query_data(
        adata,
        SCVI_PATH / "scanvi_model_ref",
    )

    vae_q.train(
        max_epochs=50,
        check_val_every_n_epoch=1,
        early_stopping=True,
        early_stopping_monitor="validation_loss",
        early_stopping_patience=10,
        early_stopping_min_delta=0.0,
        plan_kwargs={
            "lr": 5e-4,
            "weight_decay": 1e-6,
            "eps": 0.01,
            "n_epochs_kl_warmup": 10,
        },
        batch_size=512,
    )

    # predict labels
    labels_hard = vae_q.predict()
    labels_soft = vae_q.predict(soft=True)

    adata.obsm["X_embeddings_" + label] = vae_q.get_latent_representation()
    adata.obs[label] = labels_hard

    adata.write_h5ad(SCVI_PATH / "query_concat.h5ad")
    labels_soft.to_csv(
        SCVI_PATH / "query_soft_labels" + label + ".csv",
        sep="\t",
    )
