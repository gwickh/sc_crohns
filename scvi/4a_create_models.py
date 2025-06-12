import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import pickle
import scvi


METADATA = "Sample name,Diagnosis,Region code,category,Integrated_05"
IN_OBJ = "filtered/filtered.h5ad"
OUT_OBJ = "filtered/filtered.h5ad"

adata = sc.read_h5ad(IN_OBJ)
meta = METADATA.split(",")
adata.obs = adata.obs[meta]
adata.obs.to_csv("filtered/filtered_metadata.csv", sep = "\t")

# train scvi model
scvi.settings.seed = 0 
adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000, 
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key=BATCH_KEY
)

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
)

arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)

model = scvi.model.SCVI(adata, **arches_params)
model.train()

model.save("scvi_model", overwrite = True)

adata.obsm["X_" + "scvi"] = model.get_latent_representation()
denoised = model.get_normalized_expression(adata, library_size=1e6)
adata.layers["normalized_" + "scvi"] = model.get_normalized_expression(library_size=10e6)

adata.write_h5ad(OUT_OBJ)

# train scanvi model
vae_ref = scvi.model.SCVI.load("scvi_model", adata=adata)

adata_ref.obs["labels_scanvi"] = adata_ref.obs[META].values

vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    unlabeled_category="Unknown",
    labels_key="labels_scanvi"
)
vae_ref_scan.train(max_epochs=20, n_samples_per_label=100)

vae_ref_scan.save(R_MOD_OUT, overwrite = True)

adata_ref.obsm[LATENT_ID] = vae_ref.get_latent_representation()
adata_ref.write_h5ad(R_OBJ)

labels = adata_ref.obs["labels_scanvi"]
labels.to_csv(R_LAB, sep = "\t")

exit()