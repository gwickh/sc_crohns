#!/bin/bash
#SBATCH --job-name=seurat_cluster
#SBATCH --output=seurat_load_%j.out
#SBATCH --error=seurat_load_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

# Arguments
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

echo "Create the Seurat object"
Rscript "${SCRIPT_DIR}/1_seurat_load_matrices.R" \
        1> "${OUT_DIR}/1_seurat_load_matrices.STDOUT" \
        2> "${OUT_DIR}/1_seurat_load_matrices.STDERR"

exit 0