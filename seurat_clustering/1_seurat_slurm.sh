#!/bin/bash

# SLURM directives
#SBATCH --job-name=seurat_cluster
#SBATCH --output=1_clustering_%j.SLURM.stdout
#SBATCH --error=1_clustering_%j.SLURM.stderr
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk
#SBATCH --partition=ei-medium

# Argumentsß
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

source /hpc-home/yep25yan/mamba/etc/profile.d/conda.sh
conda activate r_env
scripts=(
    "1a_seurat_load_matrices.R"
    "1b_seurat_PCA.R"
    "1c_seurat_clustering.R"
    # "1c_filter_doublets.R"
    "1d_seurat_cluster_comp.R"
    "1e_seurat_UMAP.R"
)

for script in "${scripts[@]}"; do
    echo "Running ${script}"
    Rscript "${SCRIPT_DIR}/${script}" \
            1> "${OUT_DIR}/${script}.stdout" \
            2> "${OUT_DIR}/${script}.stderr"
    echo "${script} completed"
done

exit 0