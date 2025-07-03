import scanpy as sc
import os
from glob import glob

# Set static variables
MIN_CELLS = 3
MIN_FEATURES = 200
VARS_TO_REGRESS = ["S_score", "G2M_score"]  # Default should be empty if not used

# Set paths
SCRIPT_DIR = "sc_crohns/scanpy_clustering"
SCANPY_OBJECT_PATH = "project-area/data/crohns_scrnaseq/scanpy_clustering_output"
MATRIX_DIR = "project-area/data/crohns_scrnaseq/crohns_samples"
os.makedirs(SCANPY_OBJECT_PATH, exist_ok=True)

if os.path.exists(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad")):
    print("adata_merged already created, skipping")
    exit()

# Load .h5 matrices and sample names
matrix_paths = sorted(glob(os.path.join(MATRIX_DIR, "**", "filtered_feature_bc_matrix.h5"), recursive=True))
sample_names = [os.path.basename(os.path.dirname(os.path.dirname(p))) for p in matrix_paths]

# Load matrices into AnnData objects
adatas = []
for path, sample in zip(matrix_paths, sample_names):
    ad = sc.read_10x_h5(path)
    ad.var_names_make_unique()
    ad.obs["sample_id"] = sample
    ad.obs_names = [f"{sample}_{bc}" for bc in ad.obs_names]  # emulate Seurat’s add.cell.ids
    ad.raw = ad  # store raw counts
    adatas.append(ad)

# Merge into one object
adata = adatas[0].concatenate(*adatas[1:], batch_key="sample_id", batch_categories=sample_names)

# Annotate mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Filter low-quality cells
adata = adata[(adata.obs.n_genes_by_counts > 500) & 
              (adata.obs.n_genes_by_counts < 6000) & 
              (adata.obs.pct_counts_mt < 10), :]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Cell cycle scoring
# Load cell cycle gene lists (use Seurat's list)
s_genes = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS",
    "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP",
    "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
    "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"
]

g2m_genes = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
    "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA",
    "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2",
    "G2E3", "GAS2L3", "CBX5", "CENPA"
]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
print(adata.obs.columns)
# Regress out variables
if VARS_TO_REGRESS:
    sc.pp.scale(adata, max_value=10)
    sc.pp.regress_out(adata, keys=VARS_TO_REGRESS)

# Save the object
adata.write(os.path.join(SCANPY_OBJECT_PATH, "adata_merged.h5ad"))
