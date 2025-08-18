import scanpy as sc
import scvi

scvi.settings.seed = 0

ref_adata = sc.read_h5ad("project-area/data/crohns_scrnaseq/scvi_tools_output/Full_obj_raw_counts_nosoupx_v2.h5ad")

# normalize and log transform
ref_adata.layers["counts"] = ref_adata.X.copy() 
sc.pp.normalize_total(ref_adata, target_sum=1e6)
sc.pp.log1p(ref_adata)
ref_adata.raw = ref_adata

# select highly variable genes
sc.pp.highly_variable_genes(
    ref_adata,
    flavor="seurat_v3",
    n_top_genes=5000,
    subset=True,
    layer="counts",
    batch_key="Diagnosis"
)

# setup scvi model
scvi.model.SCVI.setup_anndata(
    ref_adata,
    layer="counts",
)

arches_params = {               # scArches-safe VAE parameters 
    "use_layer_norm": "both",
    "use_batch_norm": "none",
    "encode_covariates": True,
    "dropout_rate": 0.2,
    "n_layers": 2,
}

model = scvi.model.SCVI(ref_adata, **arches_params)
model.train(    
    max_epochs=200,                 
    batch_size=2048,               
    early_stopping=True,
    check_val_every_n_epoch=1)

model.save("project-area/data/crohns_scrnaseq/scvi_tools_output/scvi_model_ref", overwrite=True)

# get embeddings and normalized expression
ref_adata.obsm["X_embeddings"] = model.get_latent_representation()
ref_adata.layers["normalized_expression"] = model.get_normalized_expression(library_size=1e6)

ref_adata.write_h5ad("project-area/data/crohns_scrnaseq/scvi_tools_output/gca_ref_scvi.h5ad")
