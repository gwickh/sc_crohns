#!/usr/bin/env python3
"""Train latent space model on reference data and project query data."""

import uuid
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from ray import tune
from utils.scVI_hyperparameter_search_utils import (
    plot_learning_curves,
    scvi_hyperparameter_search,
)
from utils.scVI_train_utils import (
    load_ref_obj,
    scvi_get_embeddings_and_normalized_expression,
    scvi_train,
)
from utils.sysVI_train_utils import (
    MissingAnnDataMetadataError,
    sample_sysvi_init,
    train_sysvi,
)

import scvi

# set pandas string handling to use builtin str type, not pyarrow to avoid IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True


# Define paths
SCVI_PATH = Path("project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output")
SCVI_PATH.mkdir(parents=True, exist_ok=True)

LOG_PATH = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/scvi_autotune_log",
).resolve()
LOG_PATH.mkdir(parents=True, exist_ok=True)

GCA_OBJ_PATH = SCVI_PATH / "Full_obj_raw_counts_nosoupx_v2.h5ad"
REF_OBJ_PATH = SCVI_PATH / "obj_healthy_adult_pediatric_TIL.h5ad"

# define hyperparameters
search_space = {
    "model_params": {
        "n_hidden": tune.choice([64, 128, 256]),
        "n_layers": tune.choice([1, 2, 3]),
        "n_latent": tune.choice([10, 20, 30]),
        "dropout_rate": tune.choice([0.05, 0.1, 0.2]),
    },
    "train_params": {
        "max_epochs": 250,
        "check_val_every_n_epoch": 1,
        "early_stopping": True,
        "early_stopping_monitor": "validation_loss",
        "early_stopping_patience": 20,
        "early_stopping_min_delta": 0.0,
        "plan_kwargs": {
            "lr": tune.loguniform(3e-4, 3e-3),
            "weight_decay": tune.loguniform(1e-8, 1e-4),
            "eps": 1e-2,
            "n_epochs_kl_warmup": 20,
        },
    },
}

scheduler_kwargs = {
    "max_t": 250,
    "grace_period": 25,
    "reduction_factor": 2,
}


def generate_latent_space(adata, method, run_id) -> scvi.model.SCVI:
    """Train latent space model on reference data and project query data."""
    if method.lower() == "scvi":
        print("Training scVI model...")
        # Run hyperparameter search
        scvi_hyperparameter_search(
            adata,
            LOG_PATH,
            search_space,
            scheduler_kwargs,
        )

        plot_learning_curves(LOG_PATH, SCVI_PATH)

        # Train a final model with the best hyperparameters and get embeddings
        model = scvi_train(adata=adata)

    elif method.lower() == "sysvi":
        print("Training SysVI model...")
        cfg = sample_sysvi_init()

        model = train_sysvi(
            adata=adata,
            output_path=SCVI_PATH,
            run_id=run_id,
            **cfg,
        )

    else:
        raise NotImplementedError

    return model


def check_adata(adata: sc.AnnData) -> None:
    """Check that the adata object has the required columns."""
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

    return adata


def main() -> None:
    """Run scVI training and get embeddings."""
    adata = load_ref_obj(GCA_OBJ_PATH, REF_OBJ_PATH)
    adata = sc.read_h5ad(
        "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy/adata_umap.h5ad",
    )
    adata = check_adata(adata)

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=5000,
        batch_key="platform",
        subset=True,
    )

    run_id = uuid.uuid4().hex[:8]
    model = generate_latent_space(adata=adata, method="sysvi", run_id=run_id)

    adata_full = adata.copy()

    scvi_get_embeddings_and_normalized_expression(
        adata=adata_full,
        model=model,
        scvi_path=SCVI_PATH,
        outfile=f"{run_id}_gca_sysvi.h5ad",
    )


if __name__ == "__main__":
    main()
