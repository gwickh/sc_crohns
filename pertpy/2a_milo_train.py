import pertpy as pt
import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np

# paths
PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/Integrated_05_label" 
OUTPATH = PATH + "/../../pertpy_output/milo"

# load data and init milo
adata = sc.read_h5ad(os.path.join(PATH, "query_concat_curated.h5ad"))
milo = pt.tl.Milo()
mdata = milo.load(adata)

print("mdata:")
print(mdata)

# build KNN graph and construct neighbourhoods
print("building KNN graph and constructing neighbourhoods...")

sc.pp.neighbors(mdata["rna"], use_rep="X_embeddingsIntegrated_05")
milo.make_nhoods(mdata["rna"], prop=0.1)

print("mdata nhoods:")
print(mdata["rna"].obsm["nhoods"])

# plot neighbourhood size distribution
print("plotting neighbourhood size distribution...")

nhood_size = np.array(mdata["rna"].obsm["nhoods"].sum(0)).ravel()
plt.hist(nhood_size, bins=100)
plt.xlabel("cells in neighbourhood")
plt.ylabel("count")

nhood_distribution = os.path.join(OUTPATH, "nhood_distribution.pdf")
plt.savefig(nhood_distribution, bbox_inches="tight")
plt.close()

# count neighbourhoods
print("counting neighbourhoods...")

mdata = milo.count_nhoods(mdata, sample_col="sample_id")
    # "Milo leverages the variation in cell numbers between replicates for the same 
    # experimental condition to test for differential abundance. Therefore we have to count
    # how many cells from each sample (in this case the patient) are in each neighbourhood."

print("mdata milo modality:")
print(mdata["milo"])

# reorder diagnosis categories
mdata["rna"].obs["Diagnosis"] = mdata["rna"]\
    .obs["Diagnosis"]\
    .cat\
    .reorder_categories(["Normal", "Crohn's Disease"])  # ensure Normal is baseline

print("mdata Diagnosis order:")
print(mdata["rna"].obs["Diagnosis"].cat.categories)

# fit milo for differential abundance
print("fitting Milo model...")
milo.da_nhoods(mdata, design="~Diagnosis", solver="pydeseq2")

print("mdata milo modality DA variables:")
print(mdata["milo"].var)

# plot nhood graph
print("calculating UMAP...")
sc.tl.umap(mdata["rna"])

print("plotting nhood graph...")
milo.build_nhood_graph(mdata)
milo.plot_nhood_graph(
    mdata,
    alpha=0.2,
    min_size=1
)

nhood_graph = os.path.join(OUTPATH, "nhood_graph.pdf")
plt.savefig(nhood_graph, bbox_inches="tight")
plt.close()

# annotate nhoods 
print("annotating nhoods...")
milo.annotate_nhoods(mdata, anno_col="category")

# plot nhood annotation fraction
plt.hist(mdata["milo"].var["nhood_annotation_frac"], bins=30)
plt.xlabel("celltype fraction")
plt.ylabel("count")

nhood_annotation_frac = os.path.join(OUTPATH, "nhood_annotation_frac.pdf")
plt.savefig(nhood_annotation_frac, bbox_inches="tight")
plt.close()

# assign ambiguous neighbourhoods to 'Mixed'
print("Assigning ambiguous neighbourhoods to 'Mixed'")
mdata["milo"].var["nhood_annotation"] = mdata["milo"]\
    .var["nhood_annotation"]\
    .cat.add_categories("Mixed")

mdata["milo"].var.loc[mdata["milo"]\
    .var["nhood_annotation_frac"] < 0.6, "nhood_annotation"] = "Mixed"

plt.hist(mdata["milo"].var["nhood_annotation_frac"], bins=30)
plt.xlabel("celltype fraction")
plt.ylabel("count")

nhood_annotation_frac_mixed = os.path.join(OUTPATH, "nhood_annotation_frac_mixed.pdf")
plt.savefig(nhood_annotation_frac_mixed, bbox_inches="tight")
plt.close()

# plot DA beeswarm
milo.plot_da_beeswarm(mdata, alpha=0.1)

plt.savefig(os.path.join(OUTPATH, "da_beeswarm.pdf"), bbox_inches="tight")
plt.close()
