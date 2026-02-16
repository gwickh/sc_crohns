import scanpy as sc
import scvi
import os

scvi.settings.seed = 0 

SCVI_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output"
os.makedirs(SCVI_PATH, exist_ok=True)

adata_ref = sc.read_h5ad(os.path.join(SCVI_PATH, "gca_ref_scvi.h5ad"))
vae_ref = scvi.model.SCVI.load(os.path.join(SCVI_PATH, "scvi_model_ref"), adata=adata_ref)

label = "Integrated_05" # "category" 
print(adata_ref.obs[label].values)
adata_ref.obs["labels_scanvi"] = adata_ref.obs[label].values

vae_ref_scanvi = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    unlabeled_category="Unknown",
    labels_key="labels_scanvi"
)

vae_ref_scanvi.train(max_epochs=100, n_samples_per_label=100)
vae_ref_scanvi.save(os.path.join(SCVI_PATH, "scanvi_model_ref"), overwrite = True)

adata_ref.obsm["X_embeddings_"+label] = vae_ref_scanvi.get_latent_representation()
adata_ref.write_h5ad(os.path.join(SCVI_PATH, "gca_ref_scvi.h5ad"))

labels = adata_ref.obs["labels_scanvi"]
labels.to_csv(os.path.join(SCVI_PATH, "gca_ref_scvi_labels.csv"), sep = "\t")