#!/bin/bash

# SLURM directives
#SBATCH --job-name=seurat_cluster
#SBATCH --output=1_seurat_load_matrices_%j.SLURMout
#SBATCH --error=1_seurat_load_matrices_%j.SLURMerr
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk

# Arguments
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

echo "Create the Seurat object"

source /hpc-home/yep25yan/mamba/etc/profile.d/conda.sh
conda activate r_env

Rscript "${SCRIPT_DIR}/1_seurat_load_matrices.R" \
        1> "${OUT_DIR}/1_seurat_load_matrices.STDOUT" \
        2> "${OUT_DIR}/1_seurat_load_matrices.STDERR"

exit 0