import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy.stats import spearmanr

def compare_clusterings(seurat_labels, scanpy_labels):
    ari = adjusted_rand_score(seurat_labels, scanpy_labels)
    nmi = normalized_mutual_info_score(seurat_labels, scanpy_labels)
    print(f"Adjusted Rand Index (ARI): {ari:.3f}")
    print(f"Normalized Mutual Information (NMI): {nmi:.3f}")

def compare_embeddings(seurat_emb, scanpy_emb, emb_name="Embedding"):
    for i in range(seurat_emb.shape[1]):
        corr = np.corrcoef(seurat_emb[:, i], scanpy_emb[:, i])[0, 1]
        spearman_corr, _ = spearmanr(seurat_emb[:, i], scanpy_emb[:, i])
        print(f"{emb_name} dim {i+1}: Pearson corr = {corr:.3f}, Spearman corr = {spearman_corr:.3f}")

def plot_cell_cycle_scores(seurat_scores, scanpy_scores, score_name="S Score"):
    plt.figure(figsize=(6,6))
    plt.scatter(seurat_scores, scanpy_scores, alpha=0.5)
    plt.xlabel(f"Seurat {score_name}")
    plt.ylabel(f"Scanpy {score_name}")
    plt.title(f"Comparison of {score_name}")
    plt.grid(True)
    plt.show()

def hvg_jaccard(seurat_hvgs, scanpy_hvgs):
    seurat_set = set(seurat_hvgs)
    scanpy_set = set(scanpy_hvgs)
    jaccard_index = len(seurat_set & scanpy_set) / len(seurat_set | scanpy_set)
    print(f"HVG Jaccard index: {jaccard_index:.3f}")

# === Example usage ===

# 1. Load your data and results
# Assume you have pandas DataFrames or arrays:
# seurat_meta: DataFrame with 'barcode' and 'seurat_clusters', optionally 'S.Score' etc.
# scanpy_adata: AnnData object with .obs including 'scanpy_clusters', cell cycle scores, etc.

# Example: Align data by barcode
common_barcodes = seurat_meta['barcode'].isin(scanpy_adata.obs_names)
seurat_common = seurat_meta[common_barcodes].set_index('barcode')
scanpy_common = scanpy_adata[seurat_common.index]

# 2. Compare clustering labels
compare_clusterings(seurat_common['seurat_clusters'], scanpy_common.obs['scanpy_clusters'])

# 3. Compare embeddings (e.g., PCA stored in .obsm['X_pca'])
compare_embeddings(
    seurat_emb=seurat_pca_array,    # numpy array, shape (cells, PCs)
    scanpy_emb=scanpy_common.obsm['X_pca'],
    emb_name="PCA"
)

# 4. Compare HVGs (lists of gene names)
hvg_jaccard(seurat_hvgs_list, scanpy_hvgs_list)

# 5. Compare cell cycle scores (if available)
plot_cell_cycle_scores(
    seurat_scores=seurat_common['S.Score'], 
    scanpy_scores=scanpy_common.obs['S_score'], 
    score_name="S Phase Score"
)

plot_cell_cycle_scores(
    seurat_scores=seurat_common['G2M.Score'], 
    scanpy_scores=scanpy_common.obs['G2M_score'], 
    score_name="G2M Phase Score"
)
