import os
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from ray import tune

from scvi.utils.scVI_hyperparameter_search_utils import (
    plot_learning_curves,
    scvi_hyperparameter_search,
)
from scvi.utils.scVI_train_utils import (
    load_ref_obj,
    scvi_get_embeddings_and_normalized_expression,
    scvi_train,
)

# set pandas string handling to use builtin str type, not pyarrow to avoid anndata IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

SCVI_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output"
os.makedirs(SCVI_PATH, exist_ok=True)

LOG_PATH = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/scvi_autotune_log"
).resolve()
LOG_PATH.mkdir(parents=True, exist_ok=True)

GCA_OBJ_PATH = os.path.join(SCVI_PATH, "Full_obj_raw_counts_nosoupx_v2.h5ad")
REF_OBJ_PATH = os.path.join(SCVI_PATH, "obj_healthy_adult_pediatric_TIL.h5ad")

adata = load_ref_obj(GCA_OBJ_PATH, REF_OBJ_PATH)

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=5000,
    batch_key="batch",
    subset=True,
)

adata_full = adata.copy()

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

# # Run hyperparameter search
# scvi_hyperparameter_search(
#     adata,
#     LOG_PATH,
#     search_space,
#     scheduler_kwargs,
# )

# plot_learning_curves(LOG_PATH, SCVI_PATH)

# Train a final model with the best hyperparameters and get embeddings for downstream analysis
model = scvi_train(
    adata=adata,
    SCVI_PATH=SCVI_PATH,
    n_hidden=128,
    n_latent=20,
    n_layers=3,
    dropout_rate=0.05,
    max_epochs=100,
    lr=1e-3,
    weight_decay=5.441928187220108e-07,
    eps=1e-2,
)

scvi_get_embeddings_and_normalized_expression(
    adata=adata_full, model=model, SCVI_PATH=SCVI_PATH
)
