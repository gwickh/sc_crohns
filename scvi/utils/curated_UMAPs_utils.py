#!/usr/bin/env python3
"""Utility functions for generating curated UMAPs."""

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


def compute_umap(adata, color, outdir, reduct_name) -> None:
    """Compute and save UMAP plot colored by specified variable."""
    fig = sc.pl.umap(
        adata,
        color=color,
        frameon=True,
        legend_loc="right margin",
        return_fig=True,
        show=False,
    )

    fig.set_size_inches(6, 6)
    fig.savefig(
        outdir / ("joint_UMAP_" + reduct_name + "_diagnosis" + ".pdf"),
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )


def specify_umap_plots(
    adata,
    adata_crohns,
    adata_normal,
    color,
    legends,
    annotations,
) -> tuple:
    """Define UMAP plot parameters."""
    adata.obs[color] = adata.obs[color].astype("category")
    cats = list(adata.obs[color].cat.categories)

    if len(cats) <= 10:
        cmap = plt.get_cmap("tab10", len(cats))
        palette = {sid: cmap(i) for i, sid in enumerate(cats)}
    elif len(cats) <= 20:
        cmap = plt.get_cmap("tab20", len(cats))
        palette = {sid: cmap(i) for i, sid in enumerate(cats)}
    else:
        cmap = sns.color_palette("hls", len(cats))
        palette = {sid: cmap[i] for i, sid in enumerate(cats)}

    adata.uns[f"{color}_colors"] = [palette[sid] for sid in cats]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12, 6),
        gridspec_kw={"width_ratios": [1, 1]},
        constrained_layout=True,
    )

    if legends:
        legend_loc = "right margin"
    elif annotations:
        legend_loc = "on data"
    else:
        legend_loc = None

    sc.pl.umap(
        adata_crohns,
        color=color,
        frameon=True,
        legend_loc=legend_loc,
        title="Crohn's Disease",
        legend_fontsize="small" if annotations else None,
        legend_fontweight="normal" if annotations else None,
        legend_fontoutline=2 if annotations else None,
        ax=axes[0],
        show=False,
    )

    sc.pl.umap(
        adata_normal,
        color=color,
        frameon=True,
        legend_loc=legend_loc,
        title="Normal",
        legend_fontsize="small" if annotations else None,
        legend_fontweight="normal" if annotations else None,
        legend_fontoutline=2 if annotations else None,
        ax=axes[1],
        show=False,
    )

    return fig, axes


def compute_joint_umap(
    adata,
    reduct_name,
    outdir,
    color,
    legends,
    annotations,
) -> None:
    """Compute and save joint UMAPs."""
    xy = adata.obsm["X_umap"]

    x_min, x_max = xy[:, 0].min(), xy[:, 0].max()
    y_min, y_max = xy[:, 1].min(), xy[:, 1].max()

    x_pad = (x_max - x_min) * 0.05
    y_pad = (y_max - y_min) * 0.05

    xlim = (x_min - x_pad, x_max + x_pad)
    ylim = (y_min - y_pad, y_max + y_pad)

    mask = adata.obs["Diagnosis"] == "Crohn's Disease"
    adata_crohns = adata[mask].copy()
    adata_normal = adata[~mask].copy()

    fig, axes = specify_umap_plots(
        adata,
        adata_crohns,
        adata_normal,
        color,
        legends,
        annotations,
    )

    for ax in axes:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("auto")

    fig.savefig(
        outdir / f"joint_UMAP_{reduct_name}_{color}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )


def compute_marginal_umap(adata, reduct_name, color, legends, annotations, outdir):
    """Compute UMAP embedding for each diagnosis and plot."""

    def compute_embeddings(adata, diagnosis, reduct_name, min_dist=0.3) -> ad.AnnData:
        adata_subset = adata[adata.obs["Diagnosis"] == diagnosis].copy()

        sc.pp.neighbors(adata_subset, use_rep=reduct_name)
        sc.tl.umap(adata_subset, min_dist=min_dist, random_state=0)

        return adata_subset

    adata_crohns = compute_embeddings(adata, "Crohn's Disease", reduct_name)
    adata_normal = compute_embeddings(adata, "Normal", reduct_name)

    fig, axes = specify_umap_plots(
        adata,
        adata_crohns,
        adata_normal,
        color,
        legends,
        annotations,
    )

    def set_centered_limits(ax, xy, pad=0.5):
        """Center axis limits around data with equal aspect ratio."""
        x_min, x_max = xy[:, 0].min(), xy[:, 0].max()
        y_min, y_max = xy[:, 1].min(), xy[:, 1].max()
        x_mid, y_mid = (x_min + x_max) / 2, (y_min + y_max) / 2
        span = max(x_max - x_min, y_max - y_min) / 2 + pad
        ax.set_xlim(x_mid - span, x_mid + span)
        ax.set_ylim(y_mid - span, y_mid + span)
        ax.set_aspect("auto")
        return ax

    set_centered_limits(axes[0], adata_crohns.obsm["X_umap"])
    set_centered_limits(axes[1], adata_normal.obsm["X_umap"])

    fig.savefig(
        outdir / f"marginal_UMAP_{reduct_name}_{color}.pdf",
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.2,
    )


def compute_celltype_props(adata, labels, crohns_samples, normal_samples, outdir):
    """Compute proportions of cell types across samples."""
    df = adata.obs[["sample_id", labels]].copy()

    prop = (
        df.groupby(["sample_id", labels], observed=True)
        .size()
        .groupby(level=0, observed=True)
        .apply(lambda s: s / s.sum())
        .unstack(fill_value=0)  # samples x celltypes
        .sort_index()
    )

    if isinstance(prop.index, pd.MultiIndex):
        prop = prop.copy()
        prop.index = prop.index.get_level_values(0)
    prop.index = prop.index.astype(str)

    # color mapping: alphabetical by cell type
    cts_alpha = sorted(prop.columns)

    if len(cts_alpha) <= 10:
        cmap = plt.get_cmap("tab10", len(cts_alpha))
        color_by_ct = {ct: cmap(i) for i, ct in enumerate(cts_alpha)}
    elif len(cts_alpha) <= 20:
        cmap = plt.get_cmap("tab20", len(cts_alpha))
        color_by_ct = {ct: cmap(i) for i, ct in enumerate(cts_alpha)}
    else:
        cmap = sns.color_palette("hls", len(cts_alpha))
        color_by_ct = {ct: cmap[i] for i, ct in enumerate(cts_alpha)}

    # split panels
    crohns_subset = prop.reindex([s for s in crohns_samples if s in prop.index]).dropna(
        how="all",
    )
    normal_subset = prop.reindex([s for s in normal_samples if s in prop.index]).dropna(
        how="all",
    )

    width_ratios = [max(1, crohns_subset.shape[0]), max(1, normal_subset.shape[0])]
    fig_w = max(10, 0.6 * (crohns_subset.shape[0] + normal_subset.shape[0]) + 2)
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(fig_w, 8),
        sharey=True,
        gridspec_kw={"width_ratios": width_ratios},
    )

    for ax, sub, title in [
        (axes[0], crohns_subset, "Crohn's Disease"),
        (axes[1], normal_subset, "Normal"),
    ]:
        x = np.arange(sub.shape[0])
        bottom = np.zeros(sub.shape[0])

        # stacked bars in abundance order
        for ct in prop.columns[::-1]:
            vals = sub[ct].to_numpy() if ct in sub.columns else np.zeros(sub.shape[0])
            ax.bar(
                x,
                vals,
                bottom=bottom,
                alpha=0.75,
                label=ct,
                color=color_by_ct[ct],
                edgecolor="dimgray",
                linewidth=1,
                zorder=2,
            )
            bottom += vals

        ax.set_xticks(x)
        ax.set_xticklabels(sub.index.tolist(), rotation=45, ha="right")
        ax.set_xlabel("Sample")
        ax.set_title(title)

        ax.set_ylim(0, 1)
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.25))
        ax.grid(
            axis="y",
            which="major",
            color="black",
            linestyle="-",
            linewidth=0.75,
            zorder=3,
        )
        ax.set_axisbelow(False)  # grid above bars

        ax.spines["left"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # one legend for both panels, alphabetical order
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
    )

    fig.tight_layout()
    out_pdf = outdir / f"stacked_bar_{labels}.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
