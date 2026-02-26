import os

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick

PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output"
LABEL = "Integrated_05"
ADATA_IN = os.path.join(PATH, "query_concat.h5ad")
ADATA_OUT = os.path.join(PATH, "query_concat_umaps.h5ad")

adata = sc.read_h5ad(ADATA_IN)

REDUCT_KEY = f"X_embeddings_{LABEL}"

sid = adata.obs["sample_id"].astype(str)

CROHNS = sorted(sid[sid.str.contains("crohns", case=False, na=False)].unique().tolist())
NORMAL = sorted(sid[sid.str.contains("normal", case=False, na=False)].unique().tolist())

UMAP_MIN_DIST = 0.3
N_NEIGHBORS = 15


if REDUCT_KEY not in adata.obsm:
    alt = f"X_embeddings{LABEL}"
    if alt in adata.obsm:
        REDUCT_KEY = alt
    else:
        raise KeyError(
            f"Missing latent representation: {REDUCT_KEY} (and {alt}) in adata.obsm"
        )


def compute_umap(
    a, use_rep, neighbors_key, umap_key, n_neighbors=N_NEIGHBORS, min_dist=UMAP_MIN_DIST
):
    sc.pp.neighbors(
        a, use_rep=use_rep, n_neighbors=n_neighbors, key_added=neighbors_key
    )
    sc.tl.umap(a, neighbors_key=neighbors_key, min_dist=min_dist)
    a.obsm[umap_key] = a.obsm["X_umap"].copy()
    return a


def save_umap_pdf(a, umap_key, color, out_pdf, ncols=None):
    fig = sc.pl.embedding(
        a,
        basis=umap_key,
        color=color,
        ncols=ncols,
        frameon=True,
        return_fig=True,
        show=False,
    )
    fig.set_size_inches(6, 6)
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


adata = compute_umap(
    adata,
    use_rep=REDUCT_KEY,
    neighbors_key="neighbors_full",
    umap_key="X_umap_full",
)

save_umap_pdf(
    adata,
    umap_key="X_umap_full",
    color=LABEL,
    out_pdf=os.path.join(PATH, f"UMAP_{REDUCT_KEY}_celltypes.pdf"),
)

save_umap_pdf(
    adata,
    umap_key="X_umap_full",
    color="sample_id",
    out_pdf=os.path.join(PATH, f"UMAP_{REDUCT_KEY}_metadata.pdf"),
)

save_umap_pdf(
    adata,
    umap_key="X_umap_full",
    color="platform",
    out_pdf=os.path.join(PATH, f"UMAP_{REDUCT_KEY}_platform.pdf"),
)

adata_c = adata[adata.obs["sample_id"].isin(CROHNS)].copy()
adata_c = compute_umap(
    adata_c,
    use_rep=REDUCT_KEY,
    neighbors_key="neighbors_crohns",
    umap_key="X_umap_crohns",
)

save_umap_pdf(
    adata_c,
    umap_key="X_umap_crohns",
    color=LABEL,
    out_pdf=os.path.join(PATH, f"UMAP_Crohns_{REDUCT_KEY}_celltypes.pdf"),
)

save_umap_pdf(
    adata_c,
    umap_key="X_umap_crohns",
    color="sample_id",
    out_pdf=os.path.join(PATH, f"UMAP_Crohns_{REDUCT_KEY}_metadata.pdf"),
)

adata_n = adata[adata.obs["sample_id"].isin(NORMAL)].copy()
adata_n = compute_umap(
    adata_n,
    use_rep=REDUCT_KEY,
    neighbors_key="neighbors_normal",
    umap_key="X_umap_normal",
)

save_umap_pdf(
    adata_n,
    umap_key="X_umap_normal",
    color=LABEL,
    out_pdf=os.path.join(PATH, f"UMAP_Normal_{REDUCT_KEY}_celltypes.pdf"),
)

save_umap_pdf(
    adata_n,
    umap_key="X_umap_normal",
    color="sample_id",
    out_pdf=os.path.join(PATH, f"UMAP_Normal_{REDUCT_KEY}_metadata.pdf"),
)

adata.obsm["X_umap_crohns"] = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
adata.obsm["X_umap_normal"] = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)

crohns_mask = adata.obs["sample_id"].isin(CROHNS).values
normal_mask = adata.obs["sample_id"].isin(NORMAL).values

adata.obsm["X_umap_crohns"][crohns_mask, :] = adata_c.obsm["X_umap_crohns"].astype(
    np.float32
)
adata.obsm["X_umap_normal"][normal_mask, :] = adata_n.obsm["X_umap_normal"].astype(
    np.float32
)

assert "sample_id" in adata.obs
assert LABEL in adata.obs

df = adata.obs[["sample_id", LABEL]].copy()

prop = (
    df.groupby(["sample_id", LABEL], observed=True)
    .size()
    .groupby(level=0, observed=True)
    .apply(lambda s: s / s.sum())
    .unstack(fill_value=0)
    .sort_index()
)

if isinstance(prop.index, pd.MultiIndex):
    prop = prop.copy()
    prop.index = prop.index.get_level_values(0)

prop.index = prop.index.astype(str)

cts_alpha = sorted(prop.columns)
palette_name = "tab20" if len(cts_alpha) <= 20 else "gist_ncar"
cmap = plt.get_cmap(palette_name, len(cts_alpha))
color_by_ct = {ct: cmap(i) for i, ct in enumerate(cts_alpha)}

subA = prop.reindex([s for s in CROHNS if s in prop.index]).dropna(how="all")
subB = prop.reindex([s for s in NORMAL if s in prop.index]).dropna(how="all")

present_cols = sorted(
    set(subA.columns[subA.sum(axis=0) > 0]).union(
        set(subB.columns[subB.sum(axis=0) > 0])
    )
)
subA = subA.reindex(columns=present_cols)
subB = subB.reindex(columns=present_cols)

width_ratios = [max(1, subA.shape[0]), max(1, subB.shape[0])]
fig_w = max(10, 0.6 * (subA.shape[0] + subB.shape[0]) + 2)

fig, axes = plt.subplots(
    1, 2, figsize=(fig_w, 5), sharey=True, gridspec_kw={"width_ratios": width_ratios}
)

for ax, sub, title in [
    (axes[0], subA, "Crohn's Disease Terminal Ileum"),
    (axes[1], subB, "Normal Terminal Ileum"),
]:
    x = np.arange(sub.shape[0])
    bottom = np.zeros(sub.shape[0])

    for ct in present_cols:
        vals = sub[ct].values if ct in sub.columns else np.zeros(sub.shape[0])
        ax.bar(x, vals, bottom=bottom, label=ct, color=color_by_ct[ct], linewidth=0)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(sub.index.tolist(), rotation=0, ha="center")
    ax.set_xlabel("Sample")
    ax.set_title(title)

    ax.set_ylim(0, 1)
    ax.yaxis.set_major_locator(mtick.MultipleLocator(0.25))
    ax.grid(axis="y", which="major", linestyle="-", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

handles, labels_ = axes[0].get_legend_handles_labels()
order = np.argsort(labels_)
handles = [handles[i] for i in order]
labels_ = [labels_[i] for i in order]
ncol = min(4, max(1, len(labels_) // 10 + 1))

fig.legend(
    handles,
    labels_,
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False,
    ncol=ncol,
    title=LABEL,
)

fig.tight_layout()
fig.savefig(
    os.path.join(PATH, "stacked_bar_celltypes_by_sample_panels.pdf"),
    bbox_inches="tight",
)
plt.close(fig)

adata.write_h5ad(ADATA_OUT)
