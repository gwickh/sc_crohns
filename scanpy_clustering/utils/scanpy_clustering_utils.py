import os

import anndata as an
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_samples, silhouette_score


def clust_dr_scanpy(
    adata,
    dr,
    pcs=50,
    non_rand_sd_frac=0.5,
    K=(15, 30),
    res=(0.2, 0.5, 1.0),
    file="n_pcs.txt",
    random_state=0,
) -> an.AnnData:
    """
    Perform clustering on dimensionality reduced data
    """
    # Compute minimum std a PC must have to be considered “non-random” assuming
    #   last 10 PCs approximate random noise to exclude SVD artifacts
    stdevs = np.sqrt(np.asarray(adata.uns[f"{dr}_pca"]["variance"]))

    pcs = min(pcs, stdevs.shape[0])
    last10_idx = np.arange(max(0, pcs - 10), pcs)  # last 10 PCs (0-based)
    mean_stdev_last = stdevs[last10_idx].mean()
    min_stdev = (1.0 + non_rand_sd_frac) * mean_stdev_last

    # Loops through PCs to determine nPCs (number of informative PCs that
    #   explains > 0.05 variance with respect to the next)
    n_pcs = 0
    for i in range(pcs):
        if stdevs[i] < min_stdev:
            break
        n_pcs = i + 1

    with open(file, "w") as f:
        f.write(str(n_pcs) + "\n")

    # For each k build a NN graph using the top nPCs components
    for k in K:
        neighbors_key = f"neighbors_{dr}_k{k}"
        print(f"Computing neighbors for {dr} with k={k}...")
        sc.pp.neighbors(
            adata,
            n_neighbors=int(k),
            n_pcs=int(n_pcs),
            use_rep=f"X_{dr}",
            key_added=neighbors_key,
        )

        # For each resolution, run leiden clustering
        for r in res:
            cl_key = f"clusters_{dr}_k{k}_r{r}"
            print(f"Clustering {dr} with k={k}, r={r}...")
            sc.tl.leiden(
                adata,
                resolution=float(r),
                neighbors_key=neighbors_key,
                key_added=cl_key,
                random_state=random_state,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
    return adata


def run_clustering(
    adata,
    disps,
    n_features,
    neighbors,
    res,
    CLUSTERING_OUTPUT_PATH,
    pcs=50,
    non_rand_sd_frac=0.5,
    random_state=0,
) -> an.AnnData:
    """
    Run clustering on multiple PCA reductions
    """
    # mean.var.plot loop
    for disp in disps:
        filename = f"mean.var.plot_disp_{disp}"
        dr_name = f"pca_{filename}"
        out_file = os.path.join(CLUSTERING_OUTPUT_PATH, f"{dr_name}_num_PCs.txt")

        print(f"Clustering for {dr_name}...")
        clust_dr_scanpy(
            adata,
            dr=dr_name,
            pcs=pcs,
            non_rand_sd_frac=non_rand_sd_frac,
            K=neighbors,
            res=res,
            file=out_file,
            random_state=random_state,
        )

    # vst_top_* loop
    for n in n_features:
        filename = f"vst_top_{n}"
        dr_name = f"pca_{filename}"
        out_file = os.path.join(CLUSTERING_OUTPUT_PATH, f"{dr_name}_num_PCs.txt")

        print(f"Clustering for {dr_name}...")
        clust_dr_scanpy(
            adata,
            dr=dr_name,
            pcs=pcs,
            non_rand_sd_frac=non_rand_sd_frac,
            K=neighbors,
            res=res,
            file=out_file,
            random_state=random_state,
        )

    return adata


def compute_cluster_composition(
    adata,
    CLUSTERING_OUTPUT_PATH,
    out_prefix,
    cluster_id,
) -> None:
    """
    Compute cluster composition and save to CSV and plot as 2x2 bar plots
    (absolute and normalized by sample and cluster)
    """

    df = pd.DataFrame(
        {
            "cluster": adata.obs[cluster_id],
            "sample": adata.obs["sample_id"],
            "num": adata.obs[cluster_id].nunique(),
        }
    )

    def plot_bar(data, x, y, fill, ax, norm=False):
        # reshape to wide form for stacked bars
        plot_df = data.pivot_table(
            index=x,
            columns=fill,
            values=y,
            aggfunc="sum",
            fill_value=0,
        )

        # normalize row-wise if requested
        if norm:
            plot_df = plot_df.div(plot_df.sum(axis=1), axis=0).fillna(0)

        # plot horizontal stacked bar chart
        plot_df.plot(
            kind="barh",
            stacked=True,
            ax=ax,
            legend=False,
        )

        ax.set_xlabel("fraction of cells" if norm else "number of cells")
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelsize=15)

    # inner function for 2x2 cluster composition plots
    def _plot_cluster_composition(
        data,
        out_prefix,
        clustering_output_path,
        width=12,
        height=8,
    ):
        fig, axes = plt.subplots(2, 2, figsize=(width, height))

        # top-left
        plot_bar(
            data=data,
            x="sample",
            y="num",
            fill="cluster",
            ax=axes[0, 0],
            norm=False,
        )
        axes[0, 0].set_title("Sample by Cluster (absolute)")

        # top-right
        plot_bar(
            data=data,
            x="sample",
            y="num",
            fill="cluster",
            ax=axes[0, 1],
            norm=True,
        )
        axes[0, 1].set_title("Sample by Cluster (normalized)")

        # bottom-left
        plot_bar(
            data=data,
            x="cluster",
            y="num",
            fill="sample",
            ax=axes[1, 0],
            norm=False,
        )
        axes[1, 0].set_title("Cluster by Sample (absolute)")

        # bottom-right
        plot_bar(
            data=data,
            x="cluster",
            y="num",
            fill="sample",
            ax=axes[1, 1],
            norm=True,
        )
        axes[1, 1].set_title("Cluster by Sample (normalized)")

        # legend
        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.02, 0.5))

        plt.tight_layout(rect=(0, 0, 0.9, 1))

        out_file = os.path.join(
            clustering_output_path,
            f"{out_prefix}_cluster_composition.pdf",
        )
        plt.savefig(out_file, bbox_inches="tight")
        plt.close(fig)

    _plot_cluster_composition(
        data=df,
        out_prefix=out_prefix,
        clustering_output_path=CLUSTERING_OUTPUT_PATH,
    )


def calculate_silhouette_width(
    adata, cluster_id, data_matrix, CLUSTERING_OUTPUT_PATH
) -> tuple:
    silh_mean = silhouette_score(data_matrix, adata.obs[cluster_id])
    silh_width = silhouette_samples(data_matrix, adata.obs[cluster_id])

    fig, ax = plt.subplots(figsize=(8, 6))

    y_lower = 10
    yticks = []
    clusters = pd.Index(adata.obs[cluster_id].unique()).sort_values()

    for cluster in clusters:
        cluster_sil_vals = silh_width[(adata.obs[cluster_id].to_numpy() == cluster)]
        cluster_sil_vals.sort()

        size_cluster = cluster_sil_vals.shape[0]
        y_upper = y_lower + size_cluster

        ax.fill_betweenx(
            y=np.arange(y_lower, y_upper),
            x1=0,
            x2=cluster_sil_vals,
            alpha=0.7,
        )

        yticks.append(y_lower + 0.5 * size_cluster)
        ax.text(-0.05, y_lower + 0.5 * size_cluster, str(cluster), va="center")

        y_lower = y_upper + 10

    ax.axvline(silh_mean, linestyle="--")
    ax.set_xlabel("Silhouette width")
    ax.set_ylabel("Cluster")
    ax.set_yticks(yticks)
    ax.set_yticklabels(adata.obs[cluster_id].unique())
    ax.set_title(f"Silhouette plot (mean = {silh_mean:.3f})")

    plt.savefig(
        os.path.join(CLUSTERING_OUTPUT_PATH, f"{cluster_id}_silhouette.pdf"),
        bbox_inches="tight",
    )
    plt.close(fig)

    return silh_mean, silh_width


def clustering_summary(
    adata,
    silh_mean,
    cluster_id,
    case,
    dim,
    k,
    r,
) -> pd.DataFrame:
    summary_df = pd.DataFrame(
        [
            {
                "case": case,
                "dim": dim,
                "k": k,
                "r": r,
                "n_clusters": adata.obs[cluster_id].nunique(),
                "mean_silhouette_width": silh_mean,
            }
        ]
    )
    return summary_df


def run_clustering_analysis(
    adata,
    CLUSTERING_OUTPUT_PATH,
    disps,
    n_features,
    neighbors,
    res,
) -> None:
    cases = []

    for disp in disps:
        cases.append(f"mean.var.plot_disp_{disp}")
    for n in n_features:
        cases.append(f"vst_top_{n}")

    clustering_output = pd.DataFrame(
        columns=[
            "case",
            "dim",
            "k",
            "r",
            "n_clusters",
            "mean_silhouette_width",
        ]
    )

    for case in cases:
        dr = f"pca_{case}"
        os.makedirs(os.path.join(CLUSTERING_OUTPUT_PATH, dr), exist_ok=True)
        pc_file = os.path.join(CLUSTERING_OUTPUT_PATH, f"{dr}_num_PCs.txt")

        with open(pc_file, "r") as f:
            dim = int(f.read().strip())

        data_matrix = np.asarray(adata.obsm[f"X_{dr}"][:, :dim])

        for k in neighbors:
            for r in res:
                cluster_id = f"clusters_{dr}_k{k}_r{r}"
                if cluster_id not in adata.obs.columns:
                    raise ValueError(f"Cluster ID {cluster_id} not found in adata.obs")

                compute_cluster_composition(
                    adata=adata,
                    CLUSTERING_OUTPUT_PATH=CLUSTERING_OUTPUT_PATH,
                    out_prefix=f"{case}_{dim}_k{k}_r{r}",
                    cluster_id=cluster_id,
                )

                silh_mean, _ = calculate_silhouette_width(
                    adata=adata,
                    cluster_id=cluster_id,
                    data_matrix=data_matrix,
                    CLUSTERING_OUTPUT_PATH=CLUSTERING_OUTPUT_PATH,
                )

                summary_df = clustering_summary(
                    adata=adata,
                    silh_mean=silh_mean,
                    cluster_id=cluster_id,
                    case=case,
                    dim=dim,
                    k=k,
                    r=r,
                )

                clustering_output = pd.concat(
                    [clustering_output, summary_df], ignore_index=True
                )
