#!/bin/bash

# SLURM directives
#SBATCH --job-name=dea
#SBATCH --output=2_dea_%j.SLURM.stdout
#SBATCH --error=2_dea_%j.SLURM.stderr
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk
#SBATCH --partition=ei-medium

# Arguments
SCRIPT_DIR="$1"
OUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

source /hpc-home/yep25yan/mamba/etc/profile.d/conda.sh
conda activate r_env
scripts=(
    "2a_dea.R"
)

for script in "${scripts[@]}"; do
    echo "Running ${script}"
    Rscript "${SCRIPT_DIR}/${script}" \
            1> "${OUT_DIR}/${script}.stdout" \
            2> "${OUT_DIR}/${script}.stderr"
    echo "${script} completed"
done

exit 0