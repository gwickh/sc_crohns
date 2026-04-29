#!/usr/bin/env python3
"""Plot learning curves for training sysVI model."""

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

TUNING_DIR = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/sysvi_tuning",
)


def concatenate_sysvi_losses(
    input_dir: str = TUNING_DIR,
) -> None:
    """Concatenate sysVI loss CSV files into a single CSV."""
    dfs = []
    for file in sorted(input_dir.glob("*_sysvi_losses.csv")):
        params = file.name.split("_sysvi_losses.csv")[0]

        df = pd.read_csv(file)
        df.insert(0, "params", params)
        dfs.append(df)

    if not dfs:
        raise FileNotFoundError

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(input_dir / "sysvi_losses.csv", index=False)


def plot_sysvi_learning_curves(
    csv_file: str,
    tuning_dir: Path = TUNING_DIR,
    alpha: float = 0.7,
    linewidth: float = 1.5,
    output_file: str | None = None,
) -> None:
    """Plot learning curves from a sysVI training history CSV file."""
    df = pd.read_csv(tuning_dir / csv_file)

    # pick x-axis
    if "epoch" in df.columns:
        x_col = "epoch"
    else:
        # if epoch was saved as index and then read back unnamed
        unnamed = [c for c in df.columns if c.lower().startswith("unnamed")]
        if unnamed:
            df = df.rename(columns={unnamed[0]: "epoch"})
            x_col = "epoch"
        else:
            msg = f"""Could not find an epoch column in {csv_file}."""
            raise ValueError(msg)

    if "params" not in df.columns:
        msg = f"Expected a 'params' column in {csv_file} to differentiate trials."
        raise ValueError(msg)

    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
        "cycle_loss_validation",
        "kl_local_validation",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    axes = axes.ravel()

    for ax, metric in zip(axes, metrics, strict=False):
        if metric not in df.columns:
            ax.set_visible(False)
            continue

        for params, subdf in df.groupby("params"):
            sorted_subdf = subdf.sort_values(x_col)
            ax.plot(
                sorted_subdf[x_col],
                sorted_subdf[metric],
                label=str(params),
                alpha=alpha,
                linewidth=linewidth,
            )

        ax.set_title(metric)
        ax.set_xlabel("Epoch")
        ax.set_ylabel(metric)
        ax.grid(visible=True, alpha=0.3)

    # one shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            title="params",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
        )

    fig.suptitle("sysVI learning curves across all trials")
    fig.tight_layout()

    if output_file:
        fig.savefig(Path(csv_file).parent / output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def make_sysvi_final_epoch_summary(
    tuning_dir: Path = TUNING_DIR,
    input_csv: str = "sysvi_losses.csv",
    output_csv: str = "sysvi_final_epoch_summary.csv",
) -> None:
    """Create a summary CSV with the final epoch metrics for each sysVI trial."""
    df = pd.read_csv(tuning_dir / input_csv)

    # recover epoch column if it was saved as unnamed index
    if "epoch" not in df.columns:
        unnamed = [c for c in df.columns if c.lower().startswith("unnamed")]
        if unnamed:
            df = df.rename(columns={unnamed[0]: "epoch"})
        else:
            msg = f"Could not find an epoch column in {input_csv}."
            raise ValueError(msg)

    if "params" not in df.columns:
        msg = f"Expected a 'params' column in {input_csv}."
        raise ValueError(msg)

    metric_cols = [
        "kl_weight",
        "validation_loss",
        "elbo_validation",
        "reconstruction_loss_validation",
        "kl_local_validation",
        "kl_global_validation",
        "cycle_loss_validation",
        "train_loss_epoch",
        "elbo_train",
        "reconstruction_loss_train",
        "kl_local_train",
        "kl_global_train",
        "cycle_loss_train",
    ]

    metric_cols = [c for c in metric_cols if c in df.columns]

    # take the final epoch row for each trial
    summary = (
        df.sort_values(["params", "epoch"])
        .groupby("params", as_index=False)
        .tail(1)
        .loc[:, ["params", "epoch", *metric_cols]]
        .sort_values("elbo_validation", ascending="elbo_validation" not in metric_cols)
        .reset_index(drop=True)
    )

    summary.to_csv(tuning_dir / output_csv, index=False)


def rank_sysvi_runs(
    input_csv: str,
    tuning_dir: Path = TUNING_DIR,
) -> None:
    """Rank sysVI runs based on a composite score of validation metrics."""
    df = pd.read_csv(tuning_dir / input_csv)

    required = [
        "params",
        "elbo_validation",
        "reconstruction_loss_validation",
        "cycle_loss_validation",
        "kl_local_validation",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        msg = (
            f"Missing required columns in {input_csv} for ranking: {', '.join(missing)}"
        )
        raise ValueError(msg)

    df["rank_elbo"] = df["elbo_validation"].rank(method="min", ascending=True)
    df["rank_recon"] = df["reconstruction_loss_validation"].rank(
        method="min",
        ascending=True,
    )
    df["rank_cycle"] = df["cycle_loss_validation"].rank(method="min", ascending=False)
    df["rank_kl"] = df["kl_local_validation"].rank(method="min", ascending=True)

    df["overall_rank_score"] = (
        0.5 * df["rank_elbo"]
        + 0.2 * df["rank_recon"]
        + 0.2 * df["rank_cycle"]
        + 0.1 * df["rank_kl"]
    )

    df = df.sort_values(
        ["overall_rank_score", "rank_elbo", "rank_recon"],
        ascending=[True, True, True],
    ).reset_index(drop=True)

    df["overall_rank"] = range(1, len(df) + 1)

    # Put the most useful columns first
    first_cols = [
        "overall_rank",
        "params",
        "overall_rank_score",
        "rank_elbo",
        "rank_recon",
        "rank_cycle",
        "rank_kl",
        "elbo_validation",
        "reconstruction_loss_validation",
        "cycle_loss_validation",
        "kl_local_validation",
    ]
    other_cols = [c for c in df.columns if c not in first_cols]
    df = df[first_cols + other_cols]

    df.to_csv(tuning_dir / input_csv, index=False)


def get_top10_sysvi_learning_curves(
    tuning_dir: Path = TUNING_DIR,
    input_csv: str = "sysvi_final_epoch_summary.csv",
    output_csv: str = "sysvi_top10_learning_curves.csv",
):
    """Plot learning curves for the top 10 ranked sysVI runs."""
    summary = pd.read_csv(tuning_dir / input_csv)
    if "params" not in summary.columns:
        msg = f"'params' column not found in {tuning_dir / input_csv}"
        raise ValueError(msg)

    top10_params = summary["params"].astype(str).head(10).tolist()

    input_dir = tuning_dir
    pattern = "*_sysvi_losses.csv"
    files = input_dir.glob(pattern)

    dfs = []
    found = set()

    for file in sorted(files):
        name = file.name
        if not name.endswith("_sysvi_losses.csv"):
            continue

        params = name.split("_sysvi_losses.csv")[0]
        if params not in top10_params:
            continue

        df = pd.read_csv(file)

        # recover epoch column if it was saved as index
        if "epoch" not in df.columns:
            unnamed = [c for c in df.columns if str(c).lower().startswith("unnamed")]
            if unnamed:
                df = df.rename(columns={unnamed[0]: "epoch"})
            else:
                df = df.reset_index(names="epoch")

        df.insert(0, "params", params)
        dfs.append(df)
        found.add(params)

    (tuning_dir / "_notselected").mkdir(exist_ok=True)

    missing = set(top10_params) - found
    for f in missing:
        for m in list(tuning_dir.glob(f"{f}*")):
            if m.is_file():
                shutil.copy2(m, (tuning_dir / "_notselected" / m.name))

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(tuning_dir / output_csv, index=False)


def plot_top10_sysvi_learning_curves(
    input_csv: str,
    output_plot: str,
    tuning_dir: Path = TUNING_DIR,
) -> None:
    """Plot learning curves for the top 10 ranked sysVI runs."""
    combined = pd.read_csv(tuning_dir / input_csv)

    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
        "cycle_loss_validation",
        "kl_local_validation",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    axes = axes.ravel()

    for ax, metric in zip(axes, metrics, strict=False):
        if metric not in combined.columns:
            ax.set_visible(False)
            continue

        for params, subdf in combined.groupby("params"):
            sorted_subdf = subdf.sort_values("epoch")
            ax.plot(
                sorted_subdf["epoch"],
                sorted_subdf[metric],
                label=str(params),
                alpha=0.9,
                linewidth=1.5,
            )

        ax.set_title(metric)
        ax.set_xlabel("Epoch")
        ax.set_ylabel(metric)
        ax.grid(visible=True, alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            title="params",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
        )

    fig.suptitle("sysVI learning curves for top 10 ranked runs")
    fig.tight_layout()
    fig.savefig(tuning_dir / output_plot, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_top10_sysvi_train_val_curves(
    input_csv: str,
    output_dir: Path,
    tuning_dir: Path = TUNING_DIR,
) -> None:
    """Plot training and validation curves for the top 10 ranked sysVI runs."""
    combined = pd.read_csv(tuning_dir / input_csv)

    metric_pairs = [
        ("elbo_train", "elbo_validation"),
        ("reconstruction_loss_train", "reconstruction_loss_validation"),
        ("cycle_loss_train", "cycle_loss_validation"),
        ("kl_local_train", "kl_local_validation"),
    ]

    for _, row in combined.iterrows():
        run_id = row["run_id"]
        metrics_file = tuning_dir / f"{run_id}_sysvi_losses.csv"
        if not metrics_file.exists():
            print(f"[skip] Missing metrics file: {metrics_file}")
            continue
        df = pd.read_csv(metrics_file)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes = axes.ravel()

        for ax, (train_metric, val_metric) in zip(axes, metric_pairs, strict=False):
            if train_metric not in df.columns or val_metric not in df.columns:
                ax.text(
                    0.5,
                    0.5,
                    f"Missing:\n{train_metric}\n{val_metric}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_axis_off()
                continue

            ax.plot(df.index, df[train_metric], label="train")
            ax.plot(df.index, df[val_metric], label="validation")

            ax.set_title(train_metric.replace("_train", ""))
            ax.set_xlabel("epoch")
            ax.set_ylabel("loss")
            ax.legend()

        fig.suptitle(f"sysVI learning curves: {run_id}", fontsize=14)
        fig.tight_layout()

        out_file = output_dir / f"{run_id}_train_val_learning_curves.png"
        fig.savefig(out_file, dpi=300, bbox_inches="tight")

        plt.close(fig)


def main() -> None:
    """Collect top 10 trials and generate learning curves for sysVI tuning."""
    concatenate_sysvi_losses()

    plot_sysvi_learning_curves(
        input_csv="sysvi_losses.csv",
        output_file="sysvi_trials.png",
    )

    make_sysvi_final_epoch_summary(
        input_csv="sysvi_losses.csv",
        output_csv="sysvi_final_epoch_summary.csv",
    )

    rank_sysvi_runs(
        input_csv="sysvi_final_epoch_summary.csv",
    )

    get_top10_sysvi_learning_curves(
        input_csv="sysvi_final_epoch_summary.csv",
        output_csv="sysvi_top10_learning_curves.csv",
    )

    plot_top10_sysvi_learning_curves(
        input_csv="sysvi_top10_learning_curves.csv",
        output_plot="sysvi_top10_learning_curves.png",
    )

    plot_top10_sysvi_train_val_curves(
        input_csv="sysvi_top10_learning_curves.csv",
        output_dir="top10_train_val_curves",
    )


if __name__ == "__main__":
    main()
