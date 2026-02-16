from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
from ray import tune
from ray.tune import ExperimentAnalysis

import scvi
from scvi import autotune

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

base = Path(LOG_PATH).resolve()
exp_dir = max(
    [p for p in base.iterdir() if p.is_dir() and p.name.startswith("scvi_")],
    key=lambda p: p.stat().st_mtime,
)

ea = ExperimentAnalysis(str(exp_dir))
df = ea.dataframe(metric="validation_loss", mode="min")

# Keep the most useful columns
keep = [
    c
    for c in df.columns
    if ("config/" in c)
    or (c in ["validation_loss", "training_iteration", "time_total_s"])
]
df_out = df[keep].sort_values("validation_loss", ascending=True)

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
out = exp_dir / f"learning_curves_top{top_n}.png"
plt.savefig(out, dpi=200)
plt.close()
print("Wrote:", out)

# Individual plots for top N
for tid, df in top:
    if "validation_loss" not in df.columns:
        continue
    x = df["training_iteration"] if "training_iteration" in df.columns else df.index
    y = df["validation_loss"]

    plt.figure()
    plt.plot(x, y)
    plt.xlabel("epoch")
    plt.ylabel("validation_loss")
    plt.title(f"{tid} (final {final_val_loss(df):.2f})")
    plt.tight_layout()
    out = exp_dir / f"learning_curve_{tid}.png"
    plt.savefig(out, dpi=200)
    plt.close()

print(f"Wrote {top_n} individual learning-curve PNGs into:", exp_dir)
