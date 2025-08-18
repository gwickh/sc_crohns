import scanpy as sc
import scvi
import os
from pathlib import Path

scvi.settings.seed = 0 

PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/"
label = "category"

# create query object from concatenated .h5ad files
query_path = sorted(Path("project-area/data/crohns_samples/").glob("*/outs/filtered_feature_bc_matrix.h5"))

print(f"Found {len(matrix_paths)} matrices:")
for q in matrix_paths:
    print(" -", q)

adata_query = []
for p in matrix_paths:
    ad = sc.read_10x_h5(p) 
    ad.var_names_make_unique()
    
    # Store sample ID from parent dir of outs
    sample_id = p.parent.parent.name
    ad.obs["sample_id"] = sample_id
    
    adatas.append(ad)

adata_query = adatas[0].concatenate(
    *adatas[1:],
    batch_key="batch",                      # new column in obs
    batch_categories=[a.obs["sample_id"][0] for a in adatas]
)

adata_query.layers["counts"] = adata_query.X.copy()
adata_query.write_h5ad(os.path.join(PATH, "query_concat.h5ad"))

# set query object for scanvi labels 
if "labels_scanvi" in adata_query.obs:
    labels = adata_query.obs["labels_scanvi"] 
adata_query.obs["labels_scanvi"] = "Unknown"

# train model on query object
scvi.model.SCANVI.prepare_query_anndata(
    adata_query, 
    os.path.join(PATH, "scanvi_model_ref"),
) 
vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    os.path.join(PATH, "scanvi_model_ref"),
)

vae_q.train(
    max_epochs=200,
    plan_kwargs={"weight_decay": 0.0}, 
    check_val_every_n_epoch=10
)

# predict labels
labels_hard = vae_q.predict()
labels_soft = vae_q.predict(soft = True)

adata_query.obsm["X_embeddings"+label] = vae_q.get_latent_representation()
adata_query.obs[label] = labels_hard

if lab_exists:
    adata_query.obs["labels_scanvi"] = labels

adata_query.write_h5ad(os.path.join(PATH, "query_concat.h5ad"))
labels_soft.to_csv(os.path.join(PATH, "query_soft_labels"+label+".csv"), sep = "\t")