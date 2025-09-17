#!/bin/bash

# SLURM directives
#SBATCH --job-name=pertpy
#SBATCH --output=pertpy_%j.SLURM.stdout
#SBATCH --error=pertpy_%j.SLURM.stderr
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk
#SBATCH --partition=ei-medium

# Usage: sbatch sc_crohns/pertpy/pertpy_slurm.sh sc_crohns/pertpy/ project-area/data/crohns_scrnaseq/pertpy_output/

# Arguments
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

source ~/.bashrc
conda activate python_scripts
scripts=(
    "scCODA.py"
    "scCODA_visualisations.py"
)

for script in "${scripts[@]}"; do
    echo "Running ${script}"
    python "${SCRIPT_DIR}/${script}" \
            1> "${OUT_DIR}/${script}.stdout" \
            2> "${OUT_DIR}/${script}.stderr"
    echo "${script} completed"
done

exit 0