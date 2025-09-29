#!/bin/bash

# SLURM directives
#SBATCH --job-name=cellchat
#SBATCH --output=cellchat_%j.SLURM.stdout
#SBATCH --error=cellchat_%j.SLURM.stderr
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk
#SBATCH --partition=ei-medium

# Usage: sbatch sc_crohns/cellchat/cellchat_slurm.sh sc_crohns/cellchat/ project-area/data/crohns_scrnaseq/cellchat_output/

# Arguments
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

source ~/.bashrc
conda activate r_env
scripts=(
    "1a_cellchat_calculate_interactions.R"
    # "1b_cellchat_compare.R"
)

for script in "${scripts[@]}"; do
    echo "Running ${script}"
    Rscript "${SCRIPT_DIR}/${script}" \
            1> "${OUT_DIR}/${script}.stdout" \
            2> "${OUT_DIR}/${script}.stderr"
    echo "${script} completed"
done

exit 0