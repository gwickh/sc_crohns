#!/usr/bin/env python3
"""Utility functions for training of sysVI."""

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from ray import train, tune
from ray.train import Checkpoint, CheckpointConfig, RunConfig
from ray.tune import TuneConfig
from scvi.external import SysVI
from scvi.train import SaveCheckpoint


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


def load_adata(adata_path: str) -> sc.AnnData:
    adata = sc.read_h5ad(adata_path)

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=5000,
        batch_key="platform",
        subset=True,
    )

    col = "gene_id/gene_ids"
    # add ensembl_id column to var, using gene_ids or gene_id column if available
    sample_id = adata.obs["sample_id"].iloc[0]
    if "gene_ids" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var["gene_ids"]
    elif "gene_id" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var["gene_id"]
    else:
        sample_id = adata.obs["sample_id"].iloc[0]
        raise MissingAnnDataMetadataError(col, sample_id, table="var")

    platform = "platform"
    if platform not in adata.obs.columns:
        raise MissingAnnDataMetadataError(platform, sample_id, table="obs")

    keep = adata.var["ensembl_id"].notna()
    adata = adata[:, keep].copy()
    adata.var_names = adata.var["ensembl_id"].astype(str).to_numpy()

    return adata[adata.obs["predicted_doublet"] != "doublet"].copy()


def train_sysvi_tune(config: dict, adata: sc.AnnData, output_path: str) -> None:
    "Saves the best validation checkpoint and reports it back to Ray."

    output_path = Path(output_path)

    run_id = config.get("run_id", "trial")

    outdir = output_path / run_id
    outdir.mkdir(parents=True, exist_ok=True)

    SysVI.setup_anndata(
        adata,
        batch_key="platform",
        continuous_covariate_keys=["pct_counts_mt"],
    )

    model = SysVI(
        adata=adata,
        n_latent=config["n_latent"],
        n_prior_components=config["n_prior_components"],
        dropout_rate=config["dropout_rate"],
        n_hidden=config["n_hidden"],
        n_layers=config["n_layers"],
    )

    # scvi-tools callback: save best model by validation ELBO
    best_ckpt_callback = SaveCheckpoint(
        dirpath=outdir / "lightning_checkpoints",
        monitor="elbo_validation",
        mode="min",
        save_top_k=1,
    )

    model.train(
        max_epochs=config.get("max_epochs", 100),
        batch_size=config.get("batch_size", 256),
        early_stopping=True,
        early_stopping_monitor="elbo_validation",
        early_stopping_patience=10,
        check_val_every_n_epoch=1,
        plan_kwargs={
            "lr": config["lr"],
            "weight_decay": config["weight_decay"],
            "eps": config["eps"],
            "kl_weight": config["kl_weight"],
            "z_distance_cycle_weight": config["z_distance_cycle_weight"],
        },
        callbacks=[best_ckpt_callback],
    )

    model.save(outdir / "final_model", overwrite=True)

    history_df = pd.DataFrame({k: v.squeeze() for k, v in model.history.items()})
    history_df.to_csv(outdir / f"{run_id}_sysvi_losses.csv", index=True)

    params = {
        "run_id": run_id,
        "n_latent": config["n_latent"],
        "n_prior_components": config["n_prior_components"],
        "dropout_rate": config["dropout_rate"],
        "n_hidden": config["n_hidden"],
        "n_layers": config["n_layers"],
        "lr": config["lr"],
        "weight_decay": config["weight_decay"],
        "eps": config["eps"],
        "kl_weight": config["kl_weight"],
        "z_distance_cycle_weight": config["z_distance_cycle_weight"],
        "max_epochs": config.get("max_epochs", 100),
        "batch_size": config.get("batch_size", 256),
    }
    pd.DataFrame([params]).to_csv(outdir / f"{run_id}_sysvi_params.csv", index=False)

    best_elbo = float(history_df["elbo_validation"].min())
    best_epoch = int(history_df["elbo_validation"].idxmin())

    # Save metadata
    with open(outdir / "summary.json", "w") as f:
        json.dump(
            {
                "run_id": run_id,
                "best_elbo_validation": best_elbo,
                "best_epoch": best_epoch,
                "best_model_path": str(best_ckpt_callback.best_model_path),
                "best_model_score": (
                    float(best_ckpt_callback.best_model_score)
                    if best_ckpt_callback.best_model_score is not None
                    else None
                ),
            },
            f,
            indent=2,
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        best_scvi_dir = tmpdir / "best_scvi_checkpoint"
        best_scvi_dir.mkdir(parents=True, exist_ok=True)

        best_model_path = Path(best_ckpt_callback.best_model_path)
        if best_model_path.exists():
            shutil.copy2(best_model_path, best_scvi_dir / best_model_path.name)

        shutil.copy2(
            outdir / f"{run_id}_sysvi_params.csv", tmpdir / f"{run_id}_sysvi_params.csv"
        )
        shutil.copy2(
            outdir / f"{run_id}_sysvi_losses.csv", tmpdir / f"{run_id}_sysvi_losses.csv"
        )
        shutil.copy2(outdir / "summary.json", tmpdir / "summary.json")

        ray_ckpt = Checkpoint.from_directory(str(tmpdir))
        train.report(
            {
                "run_id": run_id,
                "elbo_validation": best_elbo,
                "best_epoch": best_epoch,
            },
            checkpoint=ray_ckpt,
        )


def make_search_space(num_samples: int) -> list[dict]:
    return [
        sample_sysvi_init() | {"run_id": f"trial_{i:03d}"} for i in range(num_samples)
    ]


def run_sysvi_tuning(
    adata: sc.AnnData,
    output_path: str,
    num_samples: int = 50,
):
    search_space = {
        "n_latent": tune.choice([10, 15, 20, 30, 40]),
        "n_prior_components": tune.choice([5, 10, 20, 30, 40]),
        "dropout_rate": tune.uniform(0.0, 0.2),
        "n_hidden": tune.choice([64, 128, 256]),
        "n_layers": tune.choice([1, 2, 3]),
        "lr": tune.loguniform(1e-4, 3e-3),
        "weight_decay": tune.loguniform(1e-8, 1e-4),
        "eps": tune.loguniform(1e-8, 1e-4),
        "kl_weight": tune.uniform(0.1, 1.0),
        "z_distance_cycle_weight": tune.uniform(0.0, 5.0),
        "max_epochs": 100,
        "batch_size": 256,
    }

    tuner = tune.Tuner(
        tune.with_parameters(
            train_sysvi_tune,
            adata=adata,
            output_path=output_path,
        ),
        param_space=search_space,
        tune_config=TuneConfig(
            metric="elbo_validation",
            mode="min",
            num_samples=num_samples,
        ),
        run_config=RunConfig(
            name="ray_tuning",
            storage_path=str(Path(output_path).resolve()),
            checkpoint_config=CheckpointConfig(
                num_to_keep=1,
                checkpoint_score_attribute="elbo_validation",
                checkpoint_score_order="min",
            ),
        ),
    )

    results = tuner.fit()
    best_result = results.get_best_result(metric="elbo_validation", mode="min")

    print("Best config:")
    print(best_result.config)
    print("Best metrics:")
    print(best_result.metrics)
    print("Best checkpoint:")
    print(best_result.checkpoint)

    return results


def main():
    adata_path = Path(
        "/hpc-home/yep25yan/project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy/adata_umap.h5ad"
    ).resolve()
    output_path = Path(
        "/hpc-home/yep25yan/project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/sysvi_tuning"
    ).resolve()

    adata = load_adata(adata_path)

    run_sysvi_tuning(
        adata=adata,
        output_path=output_path,
        num_samples=75,
    )


if __name__ == "__main__":
    main()
