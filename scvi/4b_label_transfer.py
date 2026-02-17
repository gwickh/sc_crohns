import os

import anndata as ad
import pandas as pd
import scanpy as sc

from scvi.utils.scANVI_train_classifier_utils import (
    train_scanvi_query,
    train_scanvi_ref,
)

# set pandas string handling to use builtin str type, not pyarrow to avoid anndata IO issues
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

SCVI_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output"
os.makedirs(SCVI_PATH, exist_ok=True)

ADATA_OBJ_PATH = (
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy/adata_umap.h5ad"
)

adata = sc.read_h5ad(ADATA_OBJ_PATH)
adata = adata[adata.obs["predicted_doublet"] != "doublet"].copy()

train_scanvi_ref(SCVI_PATH)
train_scanvi_query(adata, SCVI_PATH)
