import numpy as np
import scanpy as sc

ADATA_OBJ_PATH = "project-area/data/crohns_scrnaseq/10c_14n_analysis/scanpy/adata_merged_clustered.h5ad"


def main():
    adata = sc.read_h5ad(ADATA_OBJ_PATH)

    print(adata.var[["gene_ids", "gene_id", "feature_types"]].isna().sum())
    print(adata.var[["gene_ids", "gene_id", "feature_types"]].head(20))
    print(adata.var[["gene_ids", "gene_id", "feature_types"]].tail(20))

    adata_chromium = adata[adata.obs["platform"] == "10X_Chromium"]
    adata_chromium = adata_chromium[
        :, np.asarray(adata_chromium.X.sum(axis=0)).ravel() > 0
    ].copy()
    chromium_gene_set = set(adata_chromium.var["gene_id"])

    adata_parse = adata[adata.obs["platform"] == "Parse"]

    parse_gene_set = set(adata_parse.var["gene_id"])

    print(f"Parse gene set size: {len(parse_gene_set)}")
    print(f"Chromium gene set size: {len(chromium_gene_set)}")
    print(
        f"Intersection gene set size: {len(parse_gene_set.intersection(chromium_gene_set))}"
    )


if __name__ == "__main__":
    main()
