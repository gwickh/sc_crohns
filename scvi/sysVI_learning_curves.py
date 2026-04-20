#!/usr/bin/env python3
"""Plot learning curves for training sysVI model."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

tuning_dir = (
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/sysvi_tuning"
)


def concatenate_sysvi_losses(
    pattern: str,
    output_file: str,
    input_dir: str = ".",
) -> None:
    """Concatenate sysVI loss CSV files into a single CSV."""
    input_dir = Path(input_dir)

    dfs = []
    for file in sorted(input_dir.glob(pattern)):
        # take the 8-char prefix before "_sysvi_losses.csv"
        params = file.name.split("_sysvi_losses.csv")[0]

        df = pd.read_csv(file)
        df.insert(0, "params", params)
        dfs.append(df)

    if not dfs:
        raise FileNotFoundError

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(input_dir / output_file, index=False)


def plot_sysvi_learning_curves(
    csv_file,
    alpha=0.7,
    linewidth=1.5,
    output_file=None,
) -> None:
    """Plot learning curves from a sysVI training history CSV file."""
    df = pd.read_csv(csv_file)

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
            raise ValueError("Could not find an epoch column in the CSV.")

    if "params" not in df.columns:
        raise ValueError("CSV must contain a 'params' column.")

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
            subdf = subdf.sort_values(x_col)
            ax.plot(
                subdf[x_col],
                subdf[metric],
                label=str(params),
                alpha=alpha,
                linewidth=linewidth,
            )

        ax.set_title(metric)
        ax.set_xlabel("Epoch")
        ax.set_ylabel(metric)
        ax.grid(True, alpha=0.3)

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
    input_csv="all_sysvi_losses.csv",
    output_csv="all_sysvi_final_epoch_summary.csv",
):
    df = pd.read_csv(input_csv)

    # recover epoch column if it was saved as unnamed index
    if "epoch" not in df.columns:
        unnamed = [c for c in df.columns if c.lower().startswith("unnamed")]
        if unnamed:
            df = df.rename(columns={unnamed[0]: "epoch"})
        else:
            raise ValueError("Could not find an epoch column.")

    if "params" not in df.columns:
        raise ValueError("Expected a 'params' column.")

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
        .loc[:, ["params", "epoch"] + metric_cols]
        .sort_values(
            "elbo_validation",
            ascending=False if "elbo_validation" in metric_cols else True,
        )
        .reset_index(drop=True)
    )

    summary.to_csv(output_csv, index=False)


def rank_sysvi_runs(
    input_csv,
    output_csv,
):
    df = pd.read_csv(input_csv)

    required = [
        "params",
        "elbo_validation",
        "reconstruction_loss_validation",
        "cycle_loss_validation",
        "kl_local_validation",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Lower/more negative is better for ELBO and reconstruction losses
    df["rank_elbo"] = df["elbo_validation"].rank(method="min", ascending=True)
    df["rank_recon"] = df["reconstruction_loss_validation"].rank(
        method="min", ascending=True
    )

    # Higher cycle loss is preferred here as a soft guard against overcorrection
    df["rank_cycle"] = df["cycle_loss_validation"].rank(method="min", ascending=False)

    # Lower KL is preferred as a soft guard against overcompression
    df["rank_kl"] = df["kl_local_validation"].rank(method="min", ascending=True)

    # Weighted composite score: lower is better
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

    df.to_csv(output_csv, index=False)


def make_top10_sysvi_learning_curves(
    summary_csv="sysvi_final_epoch_summary_ranked.csv",
    input_dir=".",
    output_csv="sysvi_top10_learning_curves.csv",
    output_plot="sysvi_top10_learning_curves.png",
    recursive=True,
):
    # read ranked summary and get top 10 params
    summary = pd.read_csv(summary_csv)
    if "params" not in summary.columns:
        raise ValueError(f"'params' column not found in {summary_csv}")

    top10_params = summary["params"].astype(str).head(10).tolist()

    input_dir = Path(input_dir)
    pattern = "*_sysvi_losses.csv"
    files = input_dir.rglob(pattern) if recursive else input_dir.glob(pattern)

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

    if not dfs:
        raise FileNotFoundError(
            "No matching *_sysvi_losses.csv files found for the top 10 params."
        )

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(input_dir / output_csv, index=False)

    metrics = [
        "elbo_validation",
        "reconstruction_loss_validation",
        "cycle_loss_validation",
        "kl_local_validation",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    axes = axes.ravel()

    for ax, metric in zip(axes, metrics):
        if metric not in combined.columns:
            ax.set_visible(False)
            continue

        for params, subdf in combined.groupby("params"):
            subdf = subdf.sort_values("epoch")
            ax.plot(
                subdf["epoch"],
                subdf[metric],
                label=str(params),
                alpha=0.9,
                linewidth=1.5,
            )

        ax.set_title(metric)
        ax.set_xlabel("Epoch")
        ax.set_ylabel(metric)
        ax.grid(True, alpha=0.3)

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
    fig.savefig(input_dir / output_plot, dpi=300, bbox_inches="tight")
    plt.close(fig)

    missing = sorted(set(top10_params) - found)
    print(f"Saved combined CSV to: {input_dir / output_csv}")
    print(f"Saved plot to: {input_dir / output_plot}")
    if missing:
        print("Missing params with no matching losses file:")
        print(", ".join(missing))


if __name__ == "__main__":
    concatenate_sysvi_losses(
        input_dir=tuning_dir,
        pattern="*_sysvi_losses.csv",
        output_file="sysvi_losses.csv",
    )
    plot_sysvi_learning_curves(
        f"{tuning_dir}/sysvi_losses.csv",
        output_file="sysvi_trials.png",
    )

    make_sysvi_final_epoch_summary(
        input_csv=f"{tuning_dir}/sysvi_losses.csv",
        output_csv=f"{tuning_dir}/sysvi_final_epoch_summary.csv",
    )

    rank_sysvi_runs(
        input_csv=f"{tuning_dir}/sysvi_final_epoch_summary.csv",
        output_csv=f"{tuning_dir}/sysvi_final_epoch_summary_ranked.csv",
    )

    make_top10_sysvi_learning_curves(
        summary_csv=f"{tuning_dir}/sysvi_final_epoch_summary_ranked.csv",
        input_dir=tuning_dir,
        output_csv="sysvi_top10_learning_curves.csv",
        output_plot="sysvi_top10_learning_curves.png",
        recursive=True,
    )
