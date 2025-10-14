import pertpy as pt
import pandas as pd
import scanpy as sc
import numpy as np
import os
import arviz as az
import matplotlib.pyplot as plt

# vars
PATH = "project-area/data/crohns_scrnaseq/scvi_tools_output/Integrated_05_label" 
OUTPATH = PATH + "/../../pertpy_output"

# load data
adata = sc.read_h5ad(os.path.join(PATH, "query_concat_curated.h5ad"))

# create model input on category
adata.obs["Diagnosis"] = adata.obs["Diagnosis"].cat.reorder_categories(
    ["Normal", "Crohn's Disease"], ordered=True
)

scc = pt.tl.Sccoda()
mdata = scc.load(
    adata,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="category",
    sample_identifier="sample_id",
    covariate_obs=["Diagnosis"],
)

# prepare model
mdata = scc.prepare(
    mdata,
    modality_key="coda",
    formula="Diagnosis",      
    reference_cell_type="Mesenchymal"
)

# fit model
kwargs = {
    "target_accept_prob": 0.95, 
    "max_tree_depth": 12, 
    "find_heuristic_step_size": True
}

print("number of categories = ", int(adata.obs["category"].unique().size))

scc.run_nuts(
    mdata,
    modality_key="coda",
    num_samples=1000 * int(adata.obs["category"].nunique()/4),
    num_warmup=4000,
    rng_key=1,             
    **kwargs                 
)

scc.set_fdr(mdata, est_fdr=0.05, modality_key="coda")

# MCMC diagnostics
idata = scc.make_arviz(mdata, modality_key="coda", rng_key=1)
ess_bulk = az.ess(idata, method="bulk")
print("Min bulk ESS:", int(ess_bulk.to_array().min()))
n_div = int(idata.sample_stats["diverging"].sum().item()) if "diverging" in idata.sample_stats else 0
print("Divergences:", n_div)

summ = az.summary(idata, round_to=None)
summ.to_csv(os.path.join(OUTPATH, "category_sccoda_mcmc_diagnostics.csv"))

# summary ouput
scc.summary(mdata, modality_key="coda")

params = mdata.mod["coda"].uns.get("scCODA_params", {})
print(params)            # inspect
print(params.get("reference_cell_type", None))

effect_df = scc.get_effect_df(mdata, modality_key="coda").copy()

# save
os.makedirs(OUTPATH, exist_ok=True)
effect_df.to_csv(os.path.join(OUTPATH, "category_effect_df.csv")) 

# calculate credible interval in log2FC space
idata = scc.make_arviz(mdata, modality_key="coda", rng_key=1)
post  = idata.posterior

mod = mdata.mod["coda"] if hasattr(mdata, "mod") else mdata
print(mod.var_names)
cell_types = list(mod.var_names)

alpha = post["alpha"]  # (chain, draw, K-1)

if "beta" in post and "covariate" in post["beta"].dims:
    covs = list(map(str, post["beta"].coords["covariate"].values))
    if "Diagnosis" in covs:
        beta = post["beta"].sel(covariate="Diagnosis")
    else:
        pick = next((c for c in covs if "Diagnosis" in c), covs[0])
        beta = post["beta"].sel(covariate=pick)
else:
    beta = post["beta"]

# Two design points (contrast 0 vs 1)
eta0 = alpha + beta * 0   # (chain, draw, K)
eta1 = alpha + beta * 1   # (chain, draw, K)

# Softmax 
def softmax(x, axis=-1):
    x = np.asarray(x)
    x = x - np.max(x, axis=axis, keepdims=True)
    ex = np.exp(x)
    return ex / ex.sum(axis=axis, keepdims=True)

p0 = softmax(eta0.values, axis=-1)  # (chain, draw, K)
p1 = softmax(eta1.values, axis=-1)  # (chain, draw, K)

# Log2 fold changes
eps = 1e-12
log2fc = np.log2((p1 + eps) / (p0 + eps))  # (chain, draw, K)
S, K = log2fc.shape[0] * log2fc.shape[1], log2fc.shape[2]
log2fc_flat = log2fc.reshape(S, K)         # (samples, K)

# Mean + 95% HDI
def hdi_1d(samples, mass=0.95):
    x = np.sort(np.asarray(samples)); n = x.size
    m = int(np.floor(mass * n))
    if m < 1 or m >= n: return x[0], x[-1]
    j = int(np.argmin(x[m:] - x[:n-m]))
    return x[j], x[j+m]

means  = log2fc_flat.mean(axis=0)
hdi_lo = np.array([hdi_1d(log2fc_flat[:, j])[0] for j in range(K)])
hdi_hi = np.array([hdi_1d(log2fc_flat[:, j])[1] for j in range(K)])

cell_types = list(map(str, mod.var_names))
assert len(cell_types) == K, f"Expected {K} names, got {len(cell_types)}"

log2fc_df = pd.DataFrame({
    "Cell Type": cell_types,
    "log2FC_mean": means,
    "log2FC_HDI_low": hdi_lo,
    "log2FC_HDI_high": hdi_hi,
}).sort_values("log2FC_mean", ascending=False).reset_index(drop=True)

log2fc_df.to_csv(os.path.join(OUTPATH, "category_log2fc_df.csv"), index=False)

# create model input on curated
scc = pt.tl.Sccoda()
mdata = scc.load(
    adata,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="curated",
    sample_identifier="sample_id",
    covariate_obs=["Diagnosis"],
)

# prepare model
mdata = scc.prepare(
    mdata,
    modality_key="coda",
    formula="Diagnosis",      
    reference_cell_type="Mast cell"
)

# fit model
kwargs = {
    "target_accept_prob": 0.95, 
    "max_tree_depth": 12, 
    "find_heuristic_step_size": True
}

print("number of curated = ", int(adata.obs["curated"].unique().size))

scc.run_nuts(
    mdata,
    modality_key="coda",
    num_samples=1000 * int(adata.obs["curated"].nunique()/4),
    num_warmup=4000,
    rng_key=1,             
    **kwargs                 
)

scc.set_fdr(mdata, est_fdr=0.1, modality_key="coda")

# MCMC diagnostics
idata = scc.make_arviz(mdata, modality_key="coda", rng_key=1)
ess_bulk = az.ess(idata, method="bulk")
print("Min bulk ESS:", int(ess_bulk.to_array().min()))
n_div = int(idata.sample_stats["diverging"].sum().item()) if "diverging" in idata.sample_stats else 0
print("Divergences:", n_div)

summ = az.summary(idata, round_to=None)
summ.to_csv(os.path.join(OUTPATH, "curated_sccoda_mcmc_diagnostics.csv"))

# summary ouput
scc.summary(mdata, modality_key="coda")

params = mdata.mod["coda"].uns.get("scCODA_params", {})
print(params)       
print(params.get("reference_cell_type", None))

effect_df = scc.get_effect_df(mdata, modality_key="coda").copy()

# save
os.makedirs(OUTPATH, exist_ok=True)
effect_df.to_csv(os.path.join(OUTPATH, "curated_effect_df.csv")) 

# calculate credible intervals in log2FC space
idata = scc.make_arviz(mdata, modality_key="coda", rng_key=1)
post  = idata.posterior

mod = mdata.mod["coda"] if hasattr(mdata, "mod") else mdata
print(mod.var_names)
cell_types = list(mod.var_names)

alpha = post["alpha"]  # (chain, draw, K-1)

if "beta" in post and "covariate" in post["beta"].dims:
    covs = list(map(str, post["beta"].coords["covariate"].values))
    if "Diagnosis" in covs:
        beta = post["beta"].sel(covariate="Diagnosis")
    else:
        pick = next((c for c in covs if "Diagnosis" in c), covs[0])
        beta = post["beta"].sel(covariate=pick)
else:
    beta = post["beta"]

# contrast 0 vs 1
eta0 = alpha + beta * 0   # (chain, draw, K)
eta1 = alpha + beta * 1   # (chain, draw, K)

# Softmax
def softmax(x, axis=-1):
    x = np.asarray(x)
    x = x - np.max(x, axis=axis, keepdims=True)
    ex = np.exp(x)
    return ex / ex.sum(axis=axis, keepdims=True)

p0 = softmax(eta0.values, axis=-1)  # (chain, draw, K)
p1 = softmax(eta1.values, axis=-1)  # (chain, draw, K)

# Log2 fold changes
eps = 1e-12
log2fc = np.log2((p1 + eps) / (p0 + eps))  # (chain, draw, K)
S, K = log2fc.shape[0] * log2fc.shape[1], log2fc.shape[2]
log2fc_flat = log2fc.reshape(S, K)         # (samples, K)

# Mean + 95% HDI
def hdi_1d(samples, mass=0.95):
    x = np.sort(np.asarray(samples)); n = x.size
    m = int(np.floor(mass * n))
    if m < 1 or m >= n: return x[0], x[-1]
    j = int(np.argmin(x[m:] - x[:n-m]))
    return x[j], x[j+m]

means  = log2fc_flat.mean(axis=0)
hdi_lo = np.array([hdi_1d(log2fc_flat[:, j])[0] for j in range(K)])
hdi_hi = np.array([hdi_1d(log2fc_flat[:, j])[1] for j in range(K)])

cell_types = list(map(str, mod.var_names))
assert len(cell_types) == K, f"Expected {K} names, got {len(cell_types)}"

log2fc_df = pd.DataFrame({
    "Cell Type": cell_types,
    "log2FC_mean": means,
    "log2FC_HDI_low": hdi_lo,
    "log2FC_HDI_high": hdi_hi,
}).sort_values("log2FC_mean", ascending=False).reset_index(drop=True)

log2fc_df.to_csv(os.path.join(OUTPATH, "curated_log2fc_df.csv"), index=False)