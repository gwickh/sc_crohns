#!/usr/bin/env python3
"""Utils for training scVI model on reference object."""

from dataclasses import dataclass
from pathlib import Path

import scanpy as sc

import scvi


def load_ref_obj(gca_obj_path: Path, ref_obj_path: Path) -> sc.AnnData:
    """Load reference object if it exists, otherwise create reference object."""
    if ref_obj_path.exists():
        print(f"Loading existing reference object {ref_obj_path}")
        adata = sc.read_h5ad(ref_obj_path)
    elif gca_obj_path.exists():
        print(
            f"Loading GCA object from {gca_obj_path}",
            "and subsetting to healthy TIL samples",
        )

        adata = sc.read_h5ad(gca_obj_path)
        adata = adata[
            adata.obs["Diagnosis"].isin(["Healthy adult", "Pediatric healthy"])
        ]
        adata = adata[adata.obs["Region code"].isin(["TIL"])]

        adata.write_h5ad(ref_obj_path)
    else:
        error_msg = f"{gca_obj_path} not found."
        raise FileNotFoundError(error_msg)
    return adata


@dataclass
class SCVIConfig:
    """Configuration for scVI training."""

    n_hidden: int = 128
    n_latent: int = 10
    n_layers: int = 1
    dropout_rate: float = 0.1
    max_epochs: int = 400
    lr: float = 1e-3
    weight_decay: float = 1e-6
    eps: float = 0.01


config = SCVIConfig()


def scvi_train(
    adata: sc.AnnData,
    scvi_path: str,
    config: SCVIConfig,
) -> scvi.model.SCVI:
    """Train scVI model on adata with specified hyperparameters."""
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key="batch",
        continuous_covariate_keys=["pct_counts_mt"],
    )

    arches_params = {
        "use_layer_norm": "both",
        "use_batch_norm": "none",
        "encode_covariates": True,
    }

    model = scvi.model.SCVI(
        adata,
        n_hidden=config.n_hidden,
        n_latent=config.n_latent,
        n_layers=config.n_layers,
        dropout_rate=config.dropout_rate,
        **arches_params,
    )

    trainer_kwargs = {
        "early_stopping_monitor": "validation_loss",
        "early_stopping_patience": 20,
        "check_val_every_n_epoch": 1,
    }

    model.train(
        max_epochs=config.max_epochs,
        early_stopping=True,
        plan_kwargs={
            "lr": config.lr,
            "weight_decay": config.weight_decay,
            "eps": config.eps,
            "n_epochs_kl_warmup": 20,
        },
        **trainer_kwargs,
    )

    model.save(scvi_path / "scvi_model", overwrite=True)

    return model


def scvi_get_embeddings_and_normalized_expression(
    adata: sc.AnnData,
    model: scvi.model.SCVI,
    scvi_path: Path,
    outfile: str,
) -> sc.AnnData:
    """Add scVI latent embeddings to adata."""
    adata.obsm["X_embeddings"] = model.get_latent_representation(adata=adata)
    adata.write_h5ad(scvi_path / outfile)
    return adata
