import scanpy as sc
from matplotlib import pyplot as plt
import os

PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/"
label = "Integrated_05" # "category"
adata = sc.read_h5ad(os.path.join(PATH, "query_concat.h5ad"))
REDUCT_NAME = "X_embeddings"+label

# use scVI latent space for UMAP generation of full dataset coloured by high granularity labels
sc.pp.neighbors(adata, use_rep=REDUCT_NAME)
sc.tl.umap(adata, min_dist=0.3)

plt = sc.pl.umap(
    adata,
    color=label,
    ncols=2,
    frameon=True, 
    return_fig=True,
    show=False
)
plt.set_size_inches(8, 4)
plt.savefig(
    os.path.join(PATH, "UMAP_" + REDUCT_NAME + "_celltypes.pdf"),
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

plt = sc.pl.umap(
    adata,
    color="sample_id",
    ncols=4,
    frameon=True, 
    return_fig=True,
    show=False
)
plt.set_size_inches(12, 4)
plt.savefig(
    os.path.join(PATH, "UMAP_" + REDUCT_NAME + "_metadata.pdf"), 
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

# use scVI latent space for UMAP generation of Crohn's samples
selected = ["C_03", "C_08", "C_13"]
mask = adata.obs["sample_id"].isin(selected)

sc.pp.neighbors(adata[mask], use_rep=REDUCT_NAME)
sc.tl.umap(adata[mask], min_dist=0.3)

plt = sc.pl.umap(
    adata[mask],
    color=label,
    ncols=2,
    frameon=True, 
    return_fig=True,
    show=False
)
plt.set_size_inches(8, 4)
plt.savefig(
    os.path.join(PATH, "UMAP_Crohns" + REDUCT_NAME + "_celltypes.pdf"), 
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

plt = sc.pl.umap(
    adata[mask],
    color="sample_id",
    ncols=4,
    frameon=True, 
    return_fig=True,
    show=False
)
plt.set_size_inches(12, 4)
plt.savefig(
    os.path.join(PATH, "UMAP_Crohns_" + REDUCT_NAME + "_metadata.pdf"), 
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

# use scVI latent space for UMAP generation of normal samples
selected = ["N_04", "N_06", "N_07", "N_12"]
mask = adata.obs["sample_id"].isin(selected)

sc.pp.neighbors(adata[mask], use_rep=REDUCT_NAME)
sc.tl.umap(adata[mask], min_dist=0.3)

plt = sc.pl.umap(
    adata[mask],
    color=label,
    ncols=2,
    frameon=True, 
    return_fig=True,
    show=False
)
plt.set_size_inches(8, 4)
plt.savefig(
    os.path.join(PATH, "UMAP_Normal" + REDUCT_NAME + "_celltypes.pdf"),
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

plt = sc.pl.umap(
    adata[mask],
    color="sample_id",
    ncols=4,
    frameon=True, 
    return_fig=True,
    show=False
)
plt.set_size_inches(12, 4)
plt.savefig(
    os.path.join(PATH, "UMAP_Normal" + REDUCT_NAME + "_metadata.pdf"), 
    format="pdf",     
    bbox_inches="tight",    
    pad_inches=0.2
)

