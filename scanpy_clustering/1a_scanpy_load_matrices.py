import os

from utils.scanpy_qc_utils import (
    load_count_matrices,
    compute_qc_metrics,
    filter_low_count_cells,
    mad_filter,
    obtain_qc_stats,
    qc_plots,
)

from utils.scanpy_doublet_utils import run_scrublet, calculate_doublet_threshold

os.chdir("/users/yep25yan/dev")

# Set paths
MATRIX_DIR = "project-area/data/crohns_scrnaseq/10c_14n_analysis/crohns_samples"

OUTPATH = os.path.join(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis", "scanpy", "qc_stats"
)

os.makedirs(OUTPATH) if not os.path.exists(OUTPATH) else None

# Check if adata_merged already exists, if true then end script
if os.path.exists(os.path.join(OUTPATH, "adata_merged.h5ad")):
    print("adata_merged already created, skipping")
    exit()

# Load matrix paths and sample names
adata_list = load_count_matrices(MATRIX_DIR)

adata_list_raw = [adata.copy() for adata in adata_list]

# Run preprocessing steps
qc_dict = {
    "Number of genes": "n_genes_by_counts",
    "Log10 total read count": "log1p_total_counts",
    "Percent mitochondrial reads": "pct_counts_mt",
    "Percent ribosomal reads": "pct_counts_ribo",
}

adata_list = filter_low_count_cells(adata_list)

adata_list_raw = compute_qc_metrics(adata_list_raw)
adata_list = compute_qc_metrics(adata_list)

adata_list_filtered = mad_filter(adata_list, qc_dict)

adata_list_filtered = run_scrublet(adata_list_filtered)

transformation = ["probit"]

calculate_doublet_threshold(
    adata_list_filtered,
    outpath=OUTPATH,
    transformation=transformation,
    bins=50
)

obtain_qc_stats(adata_list, adata_list_filtered, qc_dict, OUTPATH)

stages_dict = {
    "Raw": adata_list_raw,
    "Background filtering": adata_list,
    "Platform-level MAD filtering": adata_list_filtered,
}
qc_plots(stages_dict, qc_dict, OUTPATH)

# # Normalize to median total counts of each sample and log1p transform
# for stages in stages_dict.values():
#     for sample in adata_list:
#         sc.pp.normalize_total(sample, target_sum=1e4)
#         sc.pp.log1p(sample)

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
