
import scanpy as sc
from matplotlib import pyplot as plt

adata = sc.read_h5ad("filtered/filtered.h5ad")

# run PCA then generate UMAP plots
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)

plt = sc.pl.umap(
    adata,
    color=["category","Integrated_05"],
    ncols=2,
    frameon=False,
    return_fig=True,
    show=False
)
plt.set_size_inches(8, 4)
plt.savefig("filtered/UMAP_PCA_celltypes.pdf", format="pdf")

plt = sc.pl.umap(
    adata,
    color=["Diagnosis", "Sample name", "Region code"],
    ncols=4,
    frameon=False, 
    return_fig=True,
    show=False
)
plt.set_size_inches(12, 4)
plt.savefig("filtered/UMAP_PCA_metadata.pdf", format="pdf")

# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scvi")
sc.tl.umap(adata, min_dist=0.3)

plt = sc.pl.umap(
    adata,
    color=["category","Integrated_05"],
    ncols=2,
    frameon=False, 
    return_fig=True,
    show=False
)
plt.set_size_inches(8, 4)
plt.savefig("filtered/UMAP_scvi_celltypes.pdf", format="pdf")

plt = sc.pl.umap(
    adata,
    color=["Diagnosis", "Sample name", "Region code"],
    ncols=4,
    frameon=False, 
    return_fig=True,
    show=False
)
plt.set_size_inches(12, 4)
plt.savefig("filtered/UMAP_" + "scvi" + "_metadata.pdf", format="pdf")

import pandas as pd

R_OBJ_PREFIX = "filtered"
META = ["category", "Integrated_05"]

df1_path = f"{R_OBJ_PREFIX}/UMAP_scvi_celltypes.pdf"
df2_path = f"{R_OBJ_PREFIX}/UMAP_scvi_metadata.pdf"

df_cat = pd.read_csv(df1_path, sep="\t")
df_ct = pd.read_csv(df2_path, sep="\t")

df = pd.DataFrame({
    META[0]: df_cat["category"],
    META[1]: df_ct["category"]
})

tab = df.value_counts().reset_index(name="count")

tab = tab[tab["count"] > 0]
tab = tab.sort_values(by=META[0])

tab.to_csv("filtered/UMAP_scvi_celltypes.tsv", sep="\t", index=False)


exit()