#!/bin/bash

# SLURM directives
#SBATCH --job-name=scanvi_cluster
#SBATCH --output=4_scanvi_%j.SLURM.stdout
#SBATCH --error=4_scanvi_%j.SLURM.stderr
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk
#SBATCH --partition=ei-medium

# Usage: sbatch sc_crohns/scvi/4_scanvi_slurm.sh sc_crohns/scvi/ project-area/data/crohns_scrnaseq/10c_14n_analysis/scvi_tools_output/

# Arguments
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

source ~/.bashrc
conda activate python_scripts
scripts=(
    "4a_scvi_ref.py"
    # "4b_scanvi_ref.py"
    # "4c_scanvi_query.py"
    # "4d_scanpy_clustering.py"
    # "4e_celltype_similarity_matrix.py"
    # "4f_curated_celltype_clustering.py"
)

for script in "${scripts[@]}"; do
    echo "Running ${script}"
    python "${SCRIPT_DIR}/${script}" \
            1> "${OUT_DIR}/${script}.stdout" \
            2> "${OUT_DIR}/${script}.stderr"
    echo "${script} completed"
done

exit 0