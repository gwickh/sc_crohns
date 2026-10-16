#!/bin/bash

# SLURM directives
#SBATCH --job-name=analysis_scripts
#SBATCH --output=cellchat_%j.SLURM.stdout
#SBATCH --error=cellchat_%j.SLURM.stderr
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk
#SBATCH --partition=ei-medium

# Usage: sbatch sc_crohns/analysis/analysis_slurm.sh sc_crohns/analysis/ project-area/data/crohns_scrnaseq/10c_14n_analysis/pseudobulk_dea/

# Arguments
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

source ~/.bashrc
conda activate r_scripts
scripts=(
    # 2_pathway_analysis.R
    3_GO_hierarchy_traversal.R
)

for script in "${scripts[@]}"; do
    echo "Running ${script}"
    Rscript "${SCRIPT_DIR}/${script}" \
            1> "${OUT_DIR}/${script}.stdout" \
            2> "${OUT_DIR}/${script}.stderr"
    echo "${script} completed"
done

exit 0