from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
from ray.tune import ExperimentAnalysis

import scvi
from scvi import autotune


def scvi_hyperparameter_search(
    adata: sc.AnnData,
    LOG_PATH: Path,
    search_space: dict,
    scheduler_kwargs: dict,
):
    """
    Perform hyperparameter search for scVI model on adata and log results to LOG_PATH.
    """
    # setup scvi model class
    model_cls = scvi.model.SCVI
    model_cls.setup_anndata(
        adata,
        batch_key="batch",
        continuous_covariate_keys=["pct_counts_mt"],
        categorical_covariate_keys=["10X"],
    )

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


def plot_learning_curves(
    LOG_PATH: Path,
    SCVI_PATH: Path,
    top_n: int = 10,
) -> None:
    """Plot learning curves for top hyperparameter search trials from LOG_PATH."""
    base = Path(LOG_PATH).resolve()
    exp_dir = max(
        [p for p in base.iterdir() if p.is_dir() and p.name.startswith("scvi_")],
        key=lambda p: p.stat().st_mtime,
    )

    ea = ExperimentAnalysis(str(exp_dir))
    df = ea.dataframe(metric="validation_loss", mode="min")

    cols = [
        c
        for c in df.columns
        if ("config/" in c)
        or (c in ["validation_loss", "training_iteration", "time_total_s"])
    ]
    df_out = df[cols].sort_values("validation_loss", ascending=True)

    csv_path = exp_dir / "autotune_results_final.csv"
    df_out.to_csv(csv_path, index=False)

    trial_dfs = ea.trial_dataframes

    # helper to get final val loss for ranking
    def final_val_loss(df):
        if "validation_loss" not in df.columns:
            return float("inf")
        return (
            df["validation_loss"].dropna().iloc[-1]
            if df["validation_loss"].notna().any()
            else float("inf")
        )

    # rank trials by final validation loss
    ranked = sorted(
        [(tid, df) for tid, df in trial_dfs.items()], key=lambda x: final_val_loss(x[1])
    )

    top_n = 10
    top = ranked[:top_n]

    # Combined plot (top N)
    plt.figure()
    for tid, df in top:
        if "validation_loss" not in df.columns:
            continue
        x = df["training_iteration"] if "training_iteration" in df.columns else df.index
        y = df["validation_loss"]
        plt.plot(x, y, label=f"{tid} (final {final_val_loss(df):.2f})")
    plt.xlabel("epoch")
    plt.ylabel("validation_loss")
    plt.title(f"Top {top_n} learning curves (validation_loss)")
    plt.legend(fontsize=6, loc="best")
    plt.tight_layout()
    out = SCVI_PATH / f"learning_curves_top{top_n}.png"
    plt.savefig(out, dpi=200)
    plt.close()
