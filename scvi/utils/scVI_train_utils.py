import os

import scanpy as sc

import scvi


def load_ref_obj(GCA_OBJ_PATH, REF_OBJ_PATH) -> sc.AnnData:
    """
    Load reference object if it exists, otherwise load GCA object and subset to healthy TIL samples to create reference object.
    """
    if os.path.exists(REF_OBJ_PATH):
        print(f"Loading existing reference object {REF_OBJ_PATH}")
        adata = sc.read_h5ad(REF_OBJ_PATH)
    elif os.path.exists(GCA_OBJ_PATH):
        print(
            f"Loading GCA object from {GCA_OBJ_PATH} and subsetting to healthy TIL samples"
        )

        adata = sc.read_h5ad(GCA_OBJ_PATH)
        adata = adata[
            adata.obs["Diagnosis"].isin(["Healthy adult", "Pediatric healthy"])
        ]
        adata = adata[adata.obs["Region code"].isin(["TIL"])]

        adata.write_h5ad(REF_OBJ_PATH)
    else:
        raise FileNotFoundError(f"{GCA_OBJ_PATH} not found.")
    return adata


def scvi_train(
    adata: sc.AnnData,
    SCVI_PATH: str,
    n_hidden: int,
    n_latent: int,
    n_layers: int,
    dropout_rate: float,
    max_epochs: int,
    lr: float,
    weight_decay: float,
    eps: float,
):
    """
    Train scVI model on adata with specified hyperparameters and save model to SCVI_PATH.
    """

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
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
        **arches_params,
    )

    trainer_kwargs = dict(
        early_stopping_monitor="validation_loss",
        early_stopping_patience=20,
        check_val_every_n_epoch=1,
    )

    model.train(
        max_epochs=max_epochs,
        early_stopping=True,
        plan_kwargs=dict(
            lr=lr,
            weight_decay=weight_decay,
            eps=eps,
            n_epochs_kl_warmup=20,
        ),
        **trainer_kwargs,
    )

    model.save(f"{SCVI_PATH}/scvi_model", overwrite=True)

    return model


def scvi_get_embeddings_and_normalized_expression(
    adata: sc.AnnData, model: scvi.model.SCVI, SCVI_PATH: str
) -> sc.AnnData:
    """
    Get scVI latent embeddings and normalized expression and add to adata.
    """
    adata.obsm["X_embeddings"] = model.get_latent_representation(adata=adata)
    adata.write_h5ad(f"{SCVI_PATH}/gca_ref_scvi.h5ad")
    return adata
