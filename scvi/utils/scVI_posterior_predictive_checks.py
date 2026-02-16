import os

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

import scvi

ADATA_PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/gca_ref_scvi.h5ad"
OUTDIR = "project-area/data/crohns_scrnaseq/scvi_tools_output/ppc_checks"
COUNTS_LAYER = "counts"
BATCH_KEY = "sample_id"
MODEL_DIR = "path/to/scvi_model"
N_SAMPLES = 25
BATCH_SIZE = 1024
MAX_BATCHES_TO_CHECK = 8
SEED = 0

os.makedirs(OUTDIR, exist_ok=True)
np.random.seed(SEED)

adata = sc.read_h5ad(ADATA_PATH)

if os.path.exists(MODEL_DIR):
    model = scvi.model.SCVI.load(MODEL_DIR, adata=adata)
else:
    layer = COUNTS_LAYER if COUNTS_LAYER in adata.layers else None
    scvi.model.SCVI.setup_anndata(
        adata, layer=layer, batch_key=BATCH_KEY if BATCH_KEY in adata.obs else None
    )
    model = scvi.model.SCVI(adata)
    model.train(max_epochs=200)
    os.makedirs(MODEL_DIR, exist_ok=True)
    model.save(MODEL_DIR, overwrite=True)


def _get_counts(a):
    X = a.layers[COUNTS_LAYER] if (COUNTS_LAYER and COUNTS_LAYER in a.layers) else a.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X)


def _hist_overlay(real, rep, xlabel, outpath, bins=60):
    plt.figure()
    plt.hist(real, bins=bins, density=True, alpha=0.5, label="observed")
    plt.hist(rep, bins=bins, density=True, alpha=0.5, label="posterior_predictive")
    plt.xlabel(xlabel)
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def _scatter(real_x, real_y, rep_x, rep_y, xlabel, ylabel, outpath):
    plt.figure()
    plt.scatter(real_x, real_y, s=6, alpha=0.35, label="observed")
    plt.scatter(rep_x, rep_y, s=6, alpha=0.35, label="posterior_predictive")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def _write_quantile_report(path, name, obs, rep):
    qs = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    q_obs = np.quantile(obs, qs)
    q_rep = np.quantile(rep, qs)
    with open(path, "a", encoding="utf-8") as f:
        f.write(f"\n{name}\n")
        f.write("q\tobs\trep\tdiff\n")
        for q, a, b in zip(qs, q_obs, q_rep):
            f.write(f"{q:.2f}\t{a:.6g}\t{b:.6g}\t{(b - a):.6g}\n")


def _ppc_summaries(X_obs, X_rep):
    lib_obs = X_obs.sum(axis=1)
    det_obs = (X_obs > 0).sum(axis=1)
    zero_obs = (X_obs == 0).mean(axis=1)

    lib_rep = X_rep.sum(axis=2).reshape(-1)
    det_rep = (X_rep > 0).sum(axis=2).reshape(-1)
    zero_rep = (X_rep == 0).mean(axis=2).reshape(-1)

    mean_obs_g = X_obs.mean(axis=0)
    zero_obs_g = (X_obs == 0).mean(axis=0)

    mean_rep_g = X_rep.mean(axis=(0, 1))
    zero_rep_g = (X_rep == 0).mean(axis=(0, 1))

    return {
        "cell": {
            "lib_obs": lib_obs,
            "det_obs": det_obs,
            "zero_obs": zero_obs,
            "lib_rep": lib_rep,
            "det_rep": det_rep,
            "zero_rep": zero_rep,
        },
        "gene": {
            "mean_obs": mean_obs_g,
            "zero_obs": zero_obs_g,
            "mean_rep": mean_rep_g,
            "zero_rep": zero_rep_g,
        },
    }


def run_ppc(a, tag):
    X_obs = _get_counts(a)
    X_rep = model.posterior_predictive_sample(
        adata=a, n_samples=N_SAMPLES, batch_size=BATCH_SIZE
    )
    if hasattr(X_rep, "toarray"):
        X_rep = X_rep.toarray()
    X_rep = np.asarray(X_rep)
    s = _ppc_summaries(X_obs, X_rep)

    _hist_overlay(
        s["cell"]["lib_obs"],
        s["cell"]["lib_rep"],
        f"Library size per cell{tag}",
        os.path.join(OUTDIR, f"ppc_cell_library_size{tag}.png"),
    )
    _hist_overlay(
        s["cell"]["det_obs"],
        s["cell"]["det_rep"],
        f"Detected genes per cell{tag}",
        os.path.join(OUTDIR, f"ppc_cell_detected_genes{tag}.png"),
    )
    _hist_overlay(
        s["cell"]["zero_obs"],
        s["cell"]["zero_rep"],
        f"Zero fraction per cell{tag}",
        os.path.join(OUTDIR, f"ppc_cell_zero_fraction{tag}.png"),
    )

    _scatter(
        s["gene"]["mean_obs"],
        s["gene"]["zero_obs"],
        s["gene"]["mean_rep"],
        s["gene"]["zero_rep"],
        "Gene mean expression",
        "Gene zero fraction",
        os.path.join(OUTDIR, f"ppc_gene_mean_vs_zero{tag}.png"),
    )

    report_path = os.path.join(OUTDIR, "ppc_quantiles.tsv")
    _write_quantile_report(
        report_path,
        f"cell_library_size{tag}",
        s["cell"]["lib_obs"],
        s["cell"]["lib_rep"],
    )
    _write_quantile_report(
        report_path,
        f"cell_detected_genes{tag}",
        s["cell"]["det_obs"],
        s["cell"]["det_rep"],
    )
    _write_quantile_report(
        report_path,
        f"cell_zero_fraction{tag}",
        s["cell"]["zero_obs"],
        s["cell"]["zero_rep"],
    )


run_ppc(adata, tag="")

if BATCH_KEY in adata.obs:
    batches = adata.obs[BATCH_KEY].astype(str).unique().tolist()[:MAX_BATCHES_TO_CHECK]
    for b in batches:
        idx = np.where(adata.obs[BATCH_KEY].astype(str).values == b)[0]
        if idx.size < 50:
            continue
        safe_b = "".join([c if c.isalnum() or c in "-_." else "_" for c in str(b)])
        run_ppc(adata[idx].copy(), tag=f"__batch_{safe_b}")
