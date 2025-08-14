import scanpy as sc
import scvi
import glob
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logging

print(scvi.__version__)
# Set up logging
logging.basicConfig(
    filename="scanvi.log",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.info

# Also print to stdout
def print_and_log(msg):
    print(msg)
    log(msg)

# Step 1: Load reference
print_and_log("🔹 Step 1: Loading reference dataset...")
ref = sc.read_h5ad("project-area/data/crohns_scrnaseq/scvi_tools_output/gca_ref.h5ad")
ref.var_names_make_unique()
ref.obs["source"] = "ref"
print_and_log(f"✅ Loaded reference: {ref.shape}")

ref.obs["cell_type"] = ref.obs["Integrated_05"]

# Step 2: Load query datasets
print_and_log("🔹 Step 2: Loading Crohn's query datasets...")
query_paths = glob.glob("project-area/data/crohns_scrnaseq/crohns_samples/*/outs/filtered_feature_bc_matrix.h5")
query_adatas = []

for path in query_paths:
    sample_id = Path(path).parts[2]  # crohns_samples/SAMPLE/outs/...
    print_and_log(f"   → Reading {sample_id} from {path}")
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    adata.obs["batch"] = sample_id
    adata.obs["cell_type"] = "Unknown"
    adata.obs["source"] = "query"
    query_adatas.append(adata)

query = sc.concat(query_adatas, label="batch", join="outer", fill_value=0)
print_and_log(f"✅ Merged query: {query.shape}")

# Step 3: Match gene space
print_and_log("🔹 Step 3: Matching genes between reference and query...")
common_genes = ref.var_names.intersection(query.var_names)
print_and_log(f"✅ Common genes: {len(common_genes)}")
ref = ref[:, common_genes].copy()
query = query[:, common_genes].copy()

# Step 4: Combine reference + query
print_and_log("🔹 Step 4: Concatenating reference and query...")
adata = ref.concatenate(query, batch_key="source", uns_merge="unique")
adata.var_names_make_unique()
print_and_log(f"✅ Combined AnnData shape: {adata.shape}")

# Step 4.5: Ensure unique names and valid categories
adata.var_names_make_unique()
adata.obs_names_make_unique()

# Fix categorical type for cell_type
import pandas as pd
all_cell_types = ref.obs["cell_type"].unique().tolist()
if "Unknown" not in all_cell_types:
    all_cell_types.append("Unknown")

adata.obs["cell_type"] = adata.obs["cell_type"].astype(pd.CategoricalDtype(categories=all_cell_types))

# Step 5: Setup scANVI
print_and_log("🔹 Step 5: Setting up SCANVI...")
scvi.model.SCANVI.setup_anndata(
    adata,
    batch_key="source",
    labels_key="cell_type",
    unlabeled_category="Unknown"
)
print_and_log("✅ SCANVI anndata setup complete.")

# Step 6: Train scANVI directly
print_and_log("🔹 Step 6: Training SCANVI model...")
scanvi = scvi.model.SCANVI(adata)
scanvi.train(max_epochs=400, train_size=0.95)
scanvi.save("project-area/data/crohns_scrnaseq/scvi_tools_output/scanvi_model", overwrite=True)
print_and_log("✅ SCANVI model training complete.")
# scanvi = scvi.model.SCANVI.load("project-area/data/crohns_scrnaseq/scvi_tools_output/scanvi_model", adata=adata)

# Step 7: Predict labels
print_and_log("🔹 Step 7: Predicting labels...")
adata.obs["scanvi_predicted_labels"] = scanvi.predict()
print_and_log("✅ Cell type prediction complete.")

# Step 8: Extract query and save
print_and_log("🔹 Step 8: Extracting and saving annotated query data...")
annotated_query = adata[adata.obs["source"] == "query"].copy()
annotated_query.write("project-area/data/crohns_scrnaseq/scvi_tools_output/scanvi_annotated_crohns.h5ad")
print_and_log("✅ Annotated Crohn's samples saved to scvi_tools_output/scanvi_annotated_crohns.h5ad")

# Step 9: Compute UMAP on latent space

print_and_log("🔹 Step 9: Computing latent space and UMAP...")
adata.obsm["X_scanvi"] = scanvi.get_latent_representation()
print_and_log("✅ Latent space extracted.")

sc.pp.neighbors(adata, use_rep="X_scanvi")
print_and_log("✅ Neighbors computed.")

sc.tl.umap(adata)
print_and_log("✅ UMAP embedding computed.")

# Save UMAP plot
fig_path = "project-area/data/crohns_scrnaseq/scvi_tools_output/scanvi_umap.png"
import matplotlib.pyplot as plt
import seaborn as sns

adata.obs["scanvi_predicted_labels"] = adata.obs["scanvi_predicted_labels"].astype("category")
unique_labels = adata.obs["scanvi_predicted_labels"].cat.categories
palette = sns.color_palette("tab20", len(unique_labels))
sc.pl.umap(
    adata,
    color="scanvi_predicted_labels",
    palette=palette,
    legend_loc="right margin",
    wspace=0.4,
    show=False
)
plt.gcf().set_size_inches(16, 8)
plt.savefig(fig_path, dpi=300)
print_and_log(f"✅ UMAP plot saved to {fig_path}")