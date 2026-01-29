import os
from pathlib import Path

import scanpy as sc
from ray import tune

import scvi
from scvi import autotune

SCVI_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output"
os.makedirs(SCVI_PATH, exist_ok=True)

LOG_PATH = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/scvi_autotune_log"
).resolve()
LOG_PATH.mkdir(parents=True, exist_ok=True)

GCA_OBJ_PATH = os.path.join(SCVI_PATH, "Full_obj_raw_counts_nosoupx_v2.h5ad")
REF_OBJ_PATH = os.path.join(SCVI_PATH, "obj_healthy_adult_pediatric_TIL.h5ad")

if os.path.exists(REF_OBJ_PATH):
    print(f"Loading existing reference object {REF_OBJ_PATH}")
    adata = sc.read_h5ad(REF_OBJ_PATH)
elif os.path.exists(GCA_OBJ_PATH):
    print(
        f"Loading GCA object from {GCA_OBJ_PATH} and subsetting to healthy TIL samples"
    )

    adata = sc.read_h5ad(GCA_OBJ_PATH)
    adata = adata[adata.obs["Diagnosis"].isin(["Healthy adult", "Pediatric healthy"])]
    adata = adata[adata.obs["Region code"].isin(["TIL"])]

    adata.write_h5ad(REF_OBJ_PATH)
else:
    raise FileNotFoundError(f"{GCA_OBJ_PATH} not found.")

# select highly variable genes
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=5000,
    batch_key="batch",
    subset=True,
)

# setup scvi model class
model_cls = scvi.model.SCVI
model_cls.setup_anndata(
    adata,
    batch_key="batch",
    continuous_covariate_keys=["pct_counts_mt"],
    categorical_covariate_keys=["10X"],
)

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

# train models
results = autotune.run_autotune(
    model_cls,
    adata,
    metrics="validation_loss",
    mode="min",
    search_space=search_space,
    scheduler="asha",
    scheduler_kwargs=scheduler_kwargs,
    num_samples=50,
    logging_dir=str(LOG_PATH),
)

print(results)

df = results.get_dataframe()
df.to_csv(os.path.join(LOG_PATH, "scvi_autotune_results.csv"), index=False)

# scvi.model.SCVI.setup_anndata(ref_adata, layer="counts", batch_key="sample_id")

# arches_params = {  # scArches-safe VAE parameters
#     "use_layer_norm": "both",
#     "use_batch_norm": "none",
#     "encode_covariates": True,
#     "dropout_rate": 0.2,
#     "n_layers": 2,
# }

# model = scvi.model.SCVI(ref_adata, **arches_params)

# history = model.history

# model.save(
#     "project-area/data/crohns_scrnaseq/scvi_tools_output/scvi_model_ref", overwrite=True
# )

# # get embeddings and normalized expression
# ref_adata.obsm["X_embeddings"] = model.get_latent_representation()
# ref_adata.layers["normalized_expression"] = model.get_normalized_expression(
#     library_size=1e6
# )

# ref_adata.write_h5ad(
#     "project-area/data/crohns_scrnaseq/scvi_tools_output/gca_ref_scvi.h5ad"
# )
