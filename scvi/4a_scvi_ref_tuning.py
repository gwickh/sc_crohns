import os
from pathlib import Path

import scanpy as sc
from ray import tune
from ray.tune import ExperimentAnalysis

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

# if os.path.exists(REF_OBJ_PATH):
#     print(f"Loading existing reference object {REF_OBJ_PATH}")
#     adata = sc.read_h5ad(REF_OBJ_PATH)
# elif os.path.exists(GCA_OBJ_PATH):
#     print(
#         f"Loading GCA object from {GCA_OBJ_PATH} and subsetting to healthy TIL samples"
#     )

#     adata = sc.read_h5ad(GCA_OBJ_PATH)
#     adata = adata[adata.obs["Diagnosis"].isin(["Healthy adult", "Pediatric healthy"])]
#     adata = adata[adata.obs["Region code"].isin(["TIL"])]

#     adata.write_h5ad(REF_OBJ_PATH)
# else:
#     raise FileNotFoundError(f"{GCA_OBJ_PATH} not found.")

# # select highly variable genes
# sc.pp.highly_variable_genes(
#     adata,
#     flavor="seurat_v3",
#     n_top_genes=5000,
#     batch_key="batch",
#     subset=True,
# )

# # setup scvi model class
# model_cls = scvi.model.SCVI
# model_cls.setup_anndata(
#     adata,
#     batch_key="batch",
#     continuous_covariate_keys=["pct_counts_mt"],
#     categorical_covariate_keys=["10X"],
# )

# # define hyperparameters
# search_space = {
#     "model_params": {
#         "n_hidden": tune.choice([64, 128, 256]),
#         "n_layers": tune.choice([1, 2, 3]),
#         "n_latent": tune.choice([10, 20, 30]),
#         "dropout_rate": tune.choice([0.05, 0.1, 0.2]),
#     },
#     "train_params": {
#         "max_epochs": 250,
#         "check_val_every_n_epoch": 1,
#         "early_stopping": True,
#         "early_stopping_monitor": "validation_loss",
#         "early_stopping_patience": 20,
#         "early_stopping_min_delta": 0.0,
#         "plan_kwargs": {
#             "lr": tune.loguniform(3e-4, 3e-3),
#             "weight_decay": tune.loguniform(1e-8, 1e-4),
#             "eps": 1e-2,
#             "n_epochs_kl_warmup": 20,
#         },
#     },
# }

# scheduler_kwargs = {
#     "max_t": 250,
#     "grace_period": 25,
#     "reduction_factor": 2,
# }

# # train models
# results = autotune.run_autotune(
#     model_cls,
#     adata,
#     metrics="validation_loss",
#     mode="min",
#     search_space=search_space,
#     scheduler="asha",
#     scheduler_kwargs=scheduler_kwargs,
#     num_samples=50,
#     logging_dir=str(LOG_PATH),
# )

# print(results)

base = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/scvi_autotune_log"
).resolve()
exp_dir = max(
    [p for p in base.iterdir() if p.is_dir() and p.name.startswith("scvi_")],
    key=lambda p: p.stat().st_mtime,
)

print("Using experiment:", exp_dir)

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
print("Wrote:", csv_path)
print(df_out.head(10))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# adjust these if your columns differ
H = "config/model_params/n_hidden"
L = "config/model_params/n_layers"
Z = "config/model_params/n_latent"
D = "config/model_params/dropout_rate"
LR = "config/train_params/plan_kwargs/lr"
WD = "config/train_params/plan_kwargs/weight_decay"

# (A) LR vs loss
plt.figure()
plt.scatter(df_out[LR], df_out["validation_loss"], s=15)
plt.xscale("log")
plt.xlabel("lr")
plt.ylabel("validation_loss")
plt.title("LR vs validation_loss")
plt.tight_layout()
plt.savefig(exp_dir / "lr_vs_loss.png", dpi=200)
plt.close()

# (B) WD vs loss
plt.figure()
plt.scatter(df_out[WD], df_out["validation_loss"], s=15)
plt.xscale("log")
plt.xlabel("weight_decay")
plt.ylabel("validation_loss")
plt.title("Weight decay vs validation_loss")
plt.tight_layout()
plt.savefig(exp_dir / "wd_vs_loss.png", dpi=200)
plt.close()

# (C) Mean loss heatmap by (n_layers, n_hidden)
pivot = df_out.pivot_table(index=L, columns=H, values="validation_loss", aggfunc="mean")
plt.figure()
plt.imshow(pivot.values, aspect="auto")
plt.xticks(range(len(pivot.columns)), pivot.columns)
plt.yticks(range(len(pivot.index)), pivot.index)
plt.xlabel("n_hidden")
plt.ylabel("n_layers")
plt.title("Mean validation_loss by (n_layers, n_hidden)")
plt.colorbar(label="mean validation_loss")
plt.tight_layout()
plt.savefig(exp_dir / "heatmap_layers_hidden.png", dpi=200)
plt.close()

print("Wrote: lr_vs_loss.png, wd_vs_loss.png, heatmap_layers_hidden.png")


# per-trial time series as dataframes
trial_dfs = ea.trial_dataframes  # dict: trial_id -> dataframe of reported metrics


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
