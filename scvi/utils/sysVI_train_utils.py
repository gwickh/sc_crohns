import os

import matplotlib.pyplot as plt
import scanpy as sc
from scvi.external import SysVI

import scvi


def train_sysvi(adata_path, output_path) -> None:
    adata = sc.read_h5ad(adata_path)

    # add ensembl_id column to var, using gene_ids or gene_id column if available
    if "gene_ids" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var["gene_ids"]
    elif "gene_id" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var["gene_id"]
    else:
        raise ValueError(
            f"Neither 'gene_ids' nor 'gene_id' column found in adata.var for sample {adata.obs['sample_id'].iloc[0]}"
        )

    if "platform" not in adata.obs.columns:
        raise ValueError(
            f"'platform' column not found in adata.obs for sample {adata.obs['sample_id'].iloc[0]}"
        )

    keep = adata.var["ensembl_id"].notna()
    adata = adata[:, keep].copy()
    adata.var_names = adata.var["ensembl_id"].astype(str).values

    SysVI.model.SCVI.setup_anndata(
        adata,
        batch_key="platform",
        continuous_covariate_keys=["pct_counts_mt"],
        categorical_covariate_keys=["batch"],
    )

    model = SysVI(
        adata=adata,
        embed_categorical_covariates=True,
    )

    model.train(
        plan_kwargs={"kl_weight": 1, "z_distance_cycle_weight": 5},
        max_epochs=200,
        check_val_every_n_epoch=1,
    )

    epochs_detail_plot = 100

    # Losses to plot
    losses = [
        "reconstruction_loss_train",
        "kl_local_train",
        "cycle_loss_train",
    ]
    fig, axs = plt.subplots(2, len(losses), figsize=(len(losses) * 3, 4))
    for ax_i, l_train in enumerate(losses):
        l_val = l_train.replace("_train", "_validation")
        l_name = l_train.replace("_train", "")
        # Change idx of epochs to start with 1
        l_val_values = model.trainer.logger.history[l_val].copy()
        l_val_values.index = l_val_values.index + 1
        l_train_values = model.trainer.logger.history[l_train].copy()
        l_train_values.index = l_train_values.index + 1
        for l_values, c, alpha, dp in [
            (l_train_values, "tab:blue", 1, epochs_detail_plot),
            (l_val_values, "tab:orange", 0.5, epochs_detail_plot),
        ]:
            axs[0, ax_i].plot(l_values.index, l_values.values.ravel(), c=c, alpha=alpha)
            axs[0, ax_i].set_title(l_name)
            axs[1, ax_i].plot(
                l_values.index[dp:], l_values.values.ravel()[dp:], c=c, alpha=alpha
            )

    fig.tight_layout()
    fig.savefig(os.path.join(output_path, "sysvi_losses.png"), dpi=200)
