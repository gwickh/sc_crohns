import scanpy as sc
import scvi
import os
from pathlib import Path

scvi.settings.seed = 0 

PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/"
label = "Integrated_05" # "category"

# create query object from concatenated .h5ad files
query_path = sorted(Path("project-area/data/crohns_scrnaseq/crohns_samples/").glob("*/outs/filtered_feature_bc_matrix.h5"))

print(f"Found {len(query_path)} matrices:")
for q in query_path:
    print(" -", q)

adata = []
for p in query_path:
    ad = sc.read_10x_h5(p) 
    ad.var_names_make_unique()
    
    # Store sample ID from parent dir of outs
    sample_id = p.parent.parent.name
    ad.obs["sample_id"] = sample_id
    
    adata.append(ad)

adata = adata[0].concatenate(
    *adata[1:],
    batch_key="batch",                      # new column in obs
    batch_categories=[a.obs["sample_id"][0] for a in adata]
)

adata.layers["counts"] = adata.X.copy()
adata.write_h5ad(os.path.join(PATH, "query_concat.h5ad"))

# set query object for scanvi labels 
if "labels_scanvi" in adata.obs:
    labels = adata.obs["labels_scanvi"] 
adata.obs["labels_scanvi"] = "Unknown"

# train model on query object
scvi.model.SCANVI.prepare_query_anndata(
    adata, 
    os.path.join(PATH, "scanvi_model_ref"),
) 
vae_q = scvi.model.SCANVI.load_query_data(
    adata,
    os.path.join(PATH, "scanvi_model_ref"),
)

vae_q.train(
    max_epochs=100,
    plan_kwargs={"weight_decay": 0.0}, 
    check_val_every_n_epoch=10
)

# predict labels
labels_hard = vae_q.predict()
labels_soft = vae_q.predict(soft = True)

adata.obsm["X_embeddings"+label] = vae_q.get_latent_representation()
adata.obs[label] = labels_hard

if "labels_scanvi" in adata.obs:
    adata.obs["labels_scanvi"] = label

adata.write_h5ad(os.path.join(PATH, "query_concat.h5ad"))
labels_soft.to_csv(os.path.join(PATH, "query_soft_labels"+label+".csv"), sep = "\t")