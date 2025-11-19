import os
from glob import glob
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation
import seaborn as sns
import matplotlib.pyplot as plt

# Set paths
SCRIPT_DIR = "sc_crohns/scanpy_clustering"
MATRIX_DIR = "project-area/data/crohns_scrnaseq/10c_14n_analysis/crohns_samples"

OUTPATH = os.path.join("project-area/data/crohns_scrnaseq/10c_14n_analysis", "scanpy", "qc_stats")
os.makedirs(OUTPATH, exist_ok=True) if not os.path.exists(OUTPATH) else None

# Check if adata_merged already exists, if true then end script
if os.path.exists(os.path.join(OUTPATH, "adata_merged.h5ad")):
    print("adata_merged already created, skipping")
    exit()

# Load matrix paths and sample names
matrix_paths = sorted(
    [p for p in glob(os.path.join(MATRIX_DIR, "**"), recursive=True) if not os.path.isdir(p)]
)
sample_names = [os.path.basename(p) for p in matrix_paths]

# Create AnnData objects
adata_list = []
for path, sample in zip(matrix_paths, sample_names):
    if ".h5ad" in sample:
        adata = sc.read_h5ad(path)
        adata.obs["platform"] = "Parse"  
    elif ".h5" in sample:
        adata = sc.read_10x_h5(path)
        adata.var["platform"] = "10X_Chromium"  
        adata.var["gene_name"] = adata.var.index
    else:
        raise ValueError(f"Unknown file type: {sample}")
    adata.obs["sample_id"] = sample                                     # add sample_id column
    adata.var_names_make_unique()                                       # make var names unique
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]      # make barcodes unique
    adata.raw = adata                                                   # store raw counts
    adata_list.append(adata)

# Compute per-sample QC metrics and mitochondrial and ribosomal genes
def compute_qc_metrics(adata_list) -> list:
    """
    Compute QC metrics for each sample in adata_list.
    """
    for sample in adata_list:
        
        # Remove empty droplets
        sc.pp.filter_cells(sample, min_counts=1000)

        # Identify mitochondrial and ribosomal genes
        sample.var["ribo"] = sample.var["gene_name"].str.startswith(("RPS", "RPL"))
        if sample.var["ribo"].sum() == 0:
            raise ValueError(f"No ribosomal genes found in {sample.obs['sample_id'][0]}")
        
        sample.var["mt"] = sample.var["gene_name"].str.startswith("MT-")
        if sample.var["mt"].sum() == 0:
            raise ValueError(f"No mitochondrial genes found in {sample.obs['sample_id'][0]}")

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(
            sample,
            qc_vars=["mt", "ribo"],
            inplace=True,
            log1p=True
        )
    return adata_list

# Create dataframe of QC metrics from anndata object list
def obtain_metrics(adata_list, qc_dict) -> pd.DataFrame:
    """
    Return long-format QC metrics for all samples.
    """
    
    rows = []

    for label, metric in qc_dict.items():
        for sample in adata_list:
            
            sample_id = sample.obs["sample_id"].iloc[0]
            diagnosis = "Crohn's Disease" if "Crohns" in sample_id else "Normal"
            
            values = sample.obs[metric].values

            rows.append(pd.DataFrame({
                "Sample ID": sample_id,
                "Diagnosis": diagnosis,
                "Metric": label,
                "Value": values,
                "Log10(cell count)": np.log10(sample.n_obs)
            }))

    return pd.concat(rows, ignore_index=True)

# Plot QC metrics
def qc_plots(metrics_df, qc_dict, stage="") -> None:
    """
    Plot QC metrics across samples
    """

    for label, metric in qc_dict.items():
        metric_df = metrics_df[metrics_df["Metric"] == label]
        # Plot violins over each metric
        plt.figure(figsize=(12, 5)) 
        sns.violinplot(
            data=metric_df,
            x="Sample ID", 
            y="Value", 
            hue="Diagnosis",
            cut=0, 
            density_norm="width",
            inner="quart",
            fill=False
        )
        plt.ylabel(metric)
        plt.title(f"{stage} {label} per cell")
        plt.xticks(rotation=90)

        plt.savefig(os.path.join(OUTPATH, f"qc_plot_{metric}_{stage}.png"), dpi=300, bbox_inches="tight")

    # Plot log10 cell counts per sample
    plt.figure(figsize=(3, 6)) 
    sns.stripplot(
        data=metrics_df.drop_duplicates(subset="Sample ID"), 
        x="Diagnosis", 
        y="Log10(cell count)",
        hue="Diagnosis",
        dodge=True, 
        alpha=0.5, 
        legend=False
    )
    sns.pointplot(    
        data=metrics_df.drop_duplicates(subset="Sample ID"), 
        x="Diagnosis", 
        y="Log10(cell count)",
        hue="Diagnosis",
        errorbar="sd", 
        capsize=0.2
    )
    plt.title(f"{stage} cell counts")
    plt.savefig(
        os.path.join(OUTPATH, f"qc_plot_log_cell_counts_per_diagnosis_{stage}.png"), 
        dpi=300, 
        bbox_inches="tight"
    )
    plt.clf()

    # Plot log10 cell counts per sample
    plt.figure(figsize=(12, 5)) 
    sns.barplot(
        data=metrics_df, 
        x="Sample ID", 
        y="Log10(cell count)",
        hue="Diagnosis",
    )
    plt.xticks(rotation=90)
    plt.title(f"{stage} cell counts")
    plt.savefig(
        os.path.join(OUTPATH, f"qc_plot_log_cell_counts_per_sample_{stage}.png"), 
        dpi=300, 
        bbox_inches="tight"
    )

# Utility to create MAD mask for a given metric
def create_mad_mask(adata_list, sample, metric, nmads=3) -> np.ndarray:

    platform = sample.obs["platform"].iloc[0]
    platform_values = []

    for s in adata_list:
        if s.obs["platform"].iloc[0] == platform:
            platform_values.append(s.obs[metric].values)

    platform_values = np.concatenate(platform_values)
    sample_values = sample.obs[metric].values

    med = np.median(platform_values)
    mad = median_abs_deviation(platform_values, scale=0.6745)

    # guard for collapse of bounds and dropping all cells when mad is near zero
    if mad < 1e-8:
        mask = np.ones(sample.n_obs, dtype=bool)    
    else:
        mask = (sample_values > med - nmads * mad) & (sample_values < med + nmads * mad)
    return mask

# Apply MAD filtering across samples
def mad_filter(adata_list, qc_dict, nmads=3) -> list:
    """
    Compute MAD boolean masks and filter per sample
    """
    adata_list_filtered = []

    for sample in adata_list:

        combined_mask = np.ones(sample.n_obs, dtype=bool)   # init combined mask to all True

        for label, metric in qc_dict.items():
            mask = create_mad_mask(adata_list, sample, metric)
            combined_mask &= mask                            # logical AND to combine masks across metrics 

        adata_filtered = sample[combined_mask, :].copy()
        adata_list_filtered.append(adata_filtered)

    return adata_list_filtered

# Get QC stats and save to CSV
def obtain_qc_stats(adata_list, adata_list_filtered, qc_dict, outpath) -> None:
    """
    Save QC stats to CSV.
    """
    adata_stats = []

    for sample, filtered in zip(adata_list, adata_list_filtered):

        sample_id = sample.obs["sample_id"].iloc[0]
        diagnosis = "Crohn's Disease" if "Crohns" in sample_id else "Normal"

        row = {
            "Sample ID": sample_id,
            "Diagnosis": diagnosis,
            "Initial number of cells": sample.n_obs,
            "Final number of cells": filtered.n_obs,
            "Proportion retained": filtered.n_obs / sample.n_obs,
        }
        for label, metric in qc_dict.items():

            values = sample.obs[metric].values
            med = np.median(values)
            mad = median_abs_deviation(values, scale=0.6745)
            mask = create_mad_mask(adata_list, sample, metric)

            row.update({
                f"{label} median": med,
                f"{label} MAD": mad,
                f"{label} passed": mask.sum(),
                f"{label} passed fraction": mask.sum() / sample.n_obs,
            })

        adata_stats.append(row)

    adata_df = pd.DataFrame(adata_stats)
    adata_df.to_csv(os.path.join(outpath, "qc_stats_mad_filtering.csv"), index=False)


# Run preprocessing steps
qc_dict = {
        "Number of genes": "n_genes_by_counts", 
        "Log10 total read count": "log1p_total_counts",
        "Percent mitochondrial reads": "pct_counts_mt",
        "Percent ribosomal reads": "pct_counts_ribo",
    }

compute_qc_metrics(adata_list)

metrics_df = obtain_metrics(adata_list, qc_dict)
qc_plots(metrics_df, qc_dict, stage="Raw")

adata_list_filtered = mad_filter(adata_list, qc_dict)

metrics_df_filtered = obtain_metrics(adata_list_filtered, qc_dict)
qc_plots(metrics_df_filtered, qc_dict, stage="Filtered")

obtain_qc_stats(adata_list, adata_list_filtered, qc_dict, OUTPATH)



# Normalize
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# # Cell cycle scoring
# # Load cell cycle gene lists (use Seurat's list)
# s_genes = [
#     "MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG", 
#     "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", 
#     "RFC2", "POLR1B", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", 
#     "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", 
#     "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", 
#     "MRPL36", "E2F8"
# ]

# g2m_genes = [
#     "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", 
#     "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", 
#     "TACC3", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", 
#     "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", 
#     "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", 
#     "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", 
#     "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", 
#     "G2E3", "GAS2L3", "CBX5", "CENPA"
# ]
# sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
# adata.obs['Phase'] = adata.obs['phase']

# sc.pp.regress_out(adata, keys=["S_score", "G2M_score"])
# sc.pp.scale(adata, max_value=10)

# # Save the object
# adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))
