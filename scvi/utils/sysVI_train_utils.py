#!/usr/bin/env python3
"""Utility functions for training of sysVI."""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scvi.external import SysVI


class MissingAnnDataMetadataError(ValueError):
    """Raised when a required column is missing from an AnnData object."""

    def __init__(self, column_name: str, sample_id: str, table: str = "obs"):
        """
        Initialize the error.

        Args:
            column_name: The name of the missing column (or description).
            sample_id: The ID of the sample being processed.
            table: The AnnData attribute where it's missing ('obs', 'var', etc.).

        """
        self.column_name = column_name
        self.sample_id = sample_id
        self.table = table

        # Construct the message inside the class to keep the call-site clean
        self.message = (
            f"Required metadata '{self.column_name}' not found in adata.{self.table} "
            f"for sample: {self.sample_id}"
        )
        super().__init__(self.message)


def sample_sysvi_init() -> dict:
    """Sample a single random sysVI hyperparameter initialisation."""
    rng = np.random.default_rng()

    return {
        # architecture
        "n_latent": int(rng.choice([10, 15, 20, 30, 40])),
        "n_prior_components": int(rng.choice([5, 10, 20, 30, 40])),
        "dropout_rate": float(rng.uniform(0.0, 0.2)),
        "n_hidden": int(rng.choice([64, 128, 256])),
        "n_layers": int(rng.choice([1, 2, 3])),
        # training / optimizer
        "lr": float(np.exp(rng.uniform(np.log(1e-4), np.log(3e-3)))),
        "weight_decay": float(np.exp(rng.uniform(np.log(1e-8), np.log(1e-4)))),
        "eps": float(np.exp(rng.uniform(np.log(1e-8), np.log(1e-4)))),
        "kl_weight": float(rng.uniform(0.1, 1.0)),
        "z_distance_cycle_weight": float(rng.uniform(0.0, 5.0)),
    }


def train_sysvi(
    adata: sc.AnnData,
    output_path: Path,
    run_id: str,
    *,
    n_latent: int,
    n_hidden: int,
    n_layers: int,
    n_prior_components: int,
    dropout_rate: float,
    lr: float,
    weight_decay: float,
    eps: float,
    kl_weight: float,
    z_distance_cycle_weight: float,
    max_epochs: int = 100,
    batch_size: int = 256,
):
    """Train sysVI on the given AnnData and save the model."""
    SysVI.setup_anndata(
        adata,
        batch_key="platform",
        continuous_covariate_keys=["pct_counts_mt"],
    )

    model = SysVI(
        adata=adata,
        n_latent=n_latent,
        n_prior_components=n_prior_components,
        dropout_rate=dropout_rate,
        n_hidden=n_hidden,
        n_layers=n_layers,
    )

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        early_stopping=True,
        early_stopping_monitor="elbo_validation",
        early_stopping_patience=20,
        check_val_every_n_epoch=1,
        plan_kwargs={
            "lr": lr,
            "weight_decay": weight_decay,
            "eps": eps,
            "kl_weight": kl_weight,
            "z_distance_cycle_weight": z_distance_cycle_weight,
        },
    )

    model.save(output_path / f"{run_id}_sysvi_model", overwrite=True)

    params = {
        "run_id": run_id,
        "n_latent": n_latent,
        "n_prior_components": n_prior_components,
        "dropout_rate": dropout_rate,
        "lr": lr,
        "weight_decay": weight_decay,
        "eps": eps,
        "kl_weight": kl_weight,
        "z_distance_cycle_weight": z_distance_cycle_weight,
        "max_epochs": max_epochs,
        "batch_size": batch_size,
    }

    pd.DataFrame([params]).to_csv(
        output_path / f"{run_id}_sysvi_params.csv",
        index=False,
    )

    history_df = pd.DataFrame({k: v.squeeze() for k, v in model.history.items()})
    history_df.to_csv(
        output_path / f"{run_id}_sysvi_losses.csv",
        index=True,
    )

    return model
