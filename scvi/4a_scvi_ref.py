import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# from utils.scVI_hyperparameter_search_utils import func
from utils.scVI_train import (
    load_ref_obj,
    scvi_get_embeddings_and_normalized_expression,
    scvi_train,
)

import scvi

# set pandas string handling to use builtin str type, not pyarrow to avoid anndata IO issues
pd.options.mode.string_storage = "python"

SCVI_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output"
os.makedirs(SCVI_PATH, exist_ok=True)

LOG_PATH = Path(
    "project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/scvi_autotune_log"
).resolve()
LOG_PATH.mkdir(parents=True, exist_ok=True)

GCA_OBJ_PATH = os.path.join(SCVI_PATH, "Full_obj_raw_counts_nosoupx_v2.h5ad")
REF_OBJ_PATH = os.path.join(SCVI_PATH, "obj_healthy_adult_pediatric_TIL.h5ad")

adata = load_ref_obj(GCA_OBJ_PATH, REF_OBJ_PATH)

adata = sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=5000,
    batch_key="batch",
    subset=True,
)

model = scvi_train(
    adata=adata,
    SCVI_PATH=SCVI_PATH,
    n_hidden=128,
    n_latent=20,
    n_layers=3,
    dropout_rate=0.05,
    max_epochs=100,
    lr=1e-3,
    weight_decay=5.441928187220108e-07,
    eps=1e-2,
)

scvi_get_embeddings_and_normalized_expression(
    adata=adata, model=model, SCVI_PATH=SCVI_PATH
)
