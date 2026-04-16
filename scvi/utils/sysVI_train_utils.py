#!/usr/bin/env python3
"""Utility functions for training of sysVI."""

from pathlib import Path

import matplotlib.pyplot as plt
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


def train_sysvi(adata: sc.AnnData, output_path: Path) -> None:
    """Train sysVI on the given adata and save the losses plot."""
    SysVI.setup_anndata(
        adata,
        batch_key="platform",
        continuous_covariate_keys=["pct_counts_mt"],
    )

    model = SysVI(
        adata=adata,
    )

    model.train(
        plan_kwargs={"kl_weight": 1, "z_distance_cycle_weight": 5},
        max_epochs=200,
        check_val_every_n_epoch=1,
    )

    model.save(output_path / "sysvi_model", overwrite=True)

    return model


def plot_learning_curves(output_path: Path, model, epochs_detail_plot=100) -> None:
    """Plot learning curves from the training logs."""
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
                l_values.index[dp:],
                l_values.values.ravel()[dp:],
                c=c,
                alpha=alpha,
            )

    fig.tight_layout()
    fig.savefig(Path(output_path) / "sysvi_losses.png", dpi=200)
