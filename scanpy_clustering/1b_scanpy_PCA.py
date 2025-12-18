# from sc_crohns.scanpy_clustering.utils.scanpy_pca_utils import (
#     pca_dispersion,
#     pca_varfeatures,
#     plot_elbow_plots,
#     plot_pca_loadings,
#     plot_variable_feature_plots,
# )
import os

import anndata as an
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

dr_list = []
vfeature_objects = {}


def run_pca(
    hvgs,
    adata,
    method_name,
    dr_name,
    PCA_OUTPUT_PATH,
    vfeature_objects,
    s_genes,
    g2m_genes,
) -> None:
    """
    Run PCA on set of HVG and store results in adata.
    """
    hv_adata = adata[:, hvgs].copy()

    features = hv_adata.var_names.tolist()
    vfeature_objects[method_name + "_vfeatures"] = features

    # Save features
    with open(os.path.join(PCA_OUTPUT_PATH, f"{method_name}_features.txt"), "w") as f:
        f.write("\n".join(features))

    sc.tl.score_genes_cell_cycle(hv_adata, s_genes=s_genes, g2m_genes=g2m_genes)
    # sc.pp.regress_out(adata, ["S_score", "G2M_score"])    # scvi does not need regressing out of cell cycle scores

    sc.pp.scale(hv_adata, max_value=10)

    # PCA
    sc.tl.pca(hv_adata, svd_solver="auto", n_comps=50, random_state=42)
    adata.obsm[f"X_{dr_name}"] = hv_adata.obsm["X_pca"]
    adata.uns[f"{dr_name}_pca"] = hv_adata.uns["pca"]

    n_genes_total = adata.n_vars
    n_pcs = hv_adata.varm["PCs"].shape[1]
    full_pcs = np.full((n_genes_total, n_pcs), np.nan)

    # Get indices of HVGs in adata.var_names
    hvg_mask = adata.var_names.isin(hv_adata.var_names)
    full_pcs[hvg_mask, :] = hv_adata.varm["PCs"]

    # Store in adata.varm
    adata.varm[f"{dr_name}_PCs"] = full_pcs


def pca_dispersion(
    adata,
    PCA_OUTPUT_PATH,
    disps,
    xmin,
    xmax,
    s_genes,
    g2m_genes,
    dr_list=dr_list,
    vfeature_objects=vfeature_objects,
) -> an.AnnData:
    """
    # Loop over dispersion cutoffs using mean.var.plot equivalent
    """
    for disp in disps:
        method_name = f"mean.var.plot_disp_{disp}"
        dr_name = f"pca_{method_name}"
        dr_list.append(dr_name)

        sc.pp.highly_variable_genes(
            adata, flavor="seurat", min_mean=xmin, max_mean=xmax, min_disp=disp
        )

        # Select top 2000 by dispersion
        hvgs = (
            adata.var[adata.var["highly_variable"]]
            .sort_values(by="dispersions_norm", ascending=False)
            .head(2000)
            .index
        )

        run_pca(
            hvgs,
            adata,
            method_name,
            dr_name,
            PCA_OUTPUT_PATH,
            vfeature_objects,
            s_genes,
            g2m_genes,
        )
    return adata


def pca_varfeatures(
    adata,
    PCA_OUTPUT_PATH,
    n_features,
    s_genes,
    g2m_genes,
    dr_list=dr_list,
    vfeature_objects=vfeature_objects,
) -> an.AnnData:
    for n in n_features:
        method_name = f"vst_top_{n}"
        dr_name = f"pca_{method_name}"
        dr_list.append(dr_name)

        sc.pp.highly_variable_genes(adata, n_top_genes=n, batch_key="sample_id")

        hvgs = adata.var.highly_variable

        run_pca(
            hvgs,
            adata,
            method_name,
            dr_name,
            PCA_OUTPUT_PATH,
            vfeature_objects,
            s_genes,
            g2m_genes,
        )
    return adata


def plot_elbow_plots(adata, PCA_OUTPUT_PATH, dr_list=dr_list) -> None:
    """
    # Plot elbow plots
    """
    for dr in dr_list:
        print(f"Plotting elbow plot for {dr}...")
        var_ratio = adata.uns[f"{dr}_pca"]["variance_ratio"]
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(18, 5))
        axs.plot(np.arange(1, len(var_ratio) + 1), var_ratio, marker="o")
        axs.set_xlabel("PC")
        axs.set_ylabel("Variance ratio")
        axs.set_title(dr, fontsize=10)
        axs.set_xticks(np.arange(1, len(var_ratio) + 1))
        plt.tight_layout()
        plt.savefig(os.path.join(PCA_OUTPUT_PATH, f"elbow_plots_{dr}.png"), dpi=300)
        plt.close()


def plot_pca_loadings(adata, PCA_OUTPUT_PATH, dr_list=dr_list) -> None:
    """
    # Plot PC loadings
    """
    for dr in dr_list:
        print(f"Plotting PCA loadings for {dr}...")
        pcs = adata.varm[f"{dr}_PCs"]  # shape: (n_genes, n_pcs)
        genes = adata.var_names

        n_pcs_to_plot = 3
        top_n = 10

        fig, axs = plt.subplots(nrows=n_pcs_to_plot, figsize=(10, 3 * n_pcs_to_plot))

        for i in range(n_pcs_to_plot):
            pc_loadings = pcs[:, i]
            loading_df = pd.DataFrame(
                {"gene": genes, "loading": pc_loadings}
            ).set_index("gene")

        # Get top positive and negative contributing genes
        top_genes = pd.concat(
            [
                loading_df.sort_values("loading", ascending=False).head(top_n),
                loading_df.sort_values("loading", ascending=True).head(top_n),
            ]
        )

        sns.barplot(
            x="loading",
            y=top_genes.index,
            data=top_genes.reset_index(),
            ax=axs[i],
            palette="vlag",
        )
        axs[i].axvline(0, color="black", linestyle="--")
        axs[i].set_title(f"{dr} | PC{i + 1}", fontsize=12)
        axs[i].set_xlabel("Loading weight")

    plt.tight_layout()
    plt.savefig(os.path.join(PCA_OUTPUT_PATH, f"{dr}_pca_loadings.png"), dpi=300)
    plt.close()


def plot_variable_feature_plots(
    adata, PCA_OUTPUT_PATH, vfeature_objects=vfeature_objects
) -> None:
    """
    # Variable feature plots
    """
    for name, features in vfeature_objects.items():
        sc.pl.highly_variable_genes(adata, show=False)
        plt.title(name, fontsize=10)
        plt.savefig(
            os.path.join(PCA_OUTPUT_PATH, f"hvg_{name}.png"),
            dpi=200,
            bbox_inches="tight",
        )
        plt.close()


SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy"

# Load existing AnnData object
adata = sc.read_h5ad(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))

# Create output directory
PCA_OUTPUT_PATH = os.path.join(SCANPY_OBJECT_PATH, "PCA_stats")
os.makedirs(PCA_OUTPUT_PATH) if not os.path.exists(PCA_OUTPUT_PATH) else None

# Cell cycle scoring
s_genes = [
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM7",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "CENPU",
    "HELLS",
    "RFC2",
    "POLR1B",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "MRPL36",
    "E2F8",
]

g2m_genes = [
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "PIMREG",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "CDCA3",
    "JPT1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
]

# Parameters for variable gene selection
disps = [0.25, 0.5, 0.75, 1]
n_features = [500, 1000, 2000, 3000]
xmin = 0.1
xmax = 10

# Loop over dispersion cutoffs using mean.var.plot equivalent
pca_dispersion(adata, PCA_OUTPUT_PATH, disps, xmin, xmax, s_genes, g2m_genes)

# Loop over nfeatures with vst method
pca_varfeatures(
    adata,
    PCA_OUTPUT_PATH,
    n_features,
    s_genes,
    g2m_genes,
)


adata.write(os.path.join(PCA_OUTPUT_PATH, "adata_merged_reduced.h5ad"))

# Plot elbow plots
plot_elbow_plots(adata, PCA_OUTPUT_PATH)

# Plot loadings (PC1-3)
n_pcs_to_plot = 5
top_n = 10

plot_pca_loadings(adata, PCA_OUTPUT_PATH)

# Variable feature plots
plot_variable_feature_plots(adata, PCA_OUTPUT_PATH)
