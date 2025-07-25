#!/bin/bash

# SLURM directives
#SBATCH --job-name=scanpy_cluster
#SBATCH --output=1_scanpy_%j.SLURM.stdout
#SBATCH --error=1_scanpy_%j.SLURM.stderr
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk
#SBATCH --partition=ei-short

# Arguments
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

source ~/.bashrc
conda activate python_scripts
scripts=(
    "1a_scanpy_load_matrices.py"
    "1b_scanpy_PCA.py"
)

for script in "${scripts[@]}"; do
    echo "Running ${script}"
    python "${SCRIPT_DIR}/${script}" \
            1> "${OUT_DIR}/${script}.stdout" \
            2> "${OUT_DIR}/${script}.stderr"
    echo "${script} completed"
done

exit 0