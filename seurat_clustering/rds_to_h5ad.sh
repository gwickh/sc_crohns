#!/bin/bash

# SLURM directives
#SBATCH --job-name=convert_rds_to_h5ad
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yep25yan@nbi.ac.uk
#SBATCH --partition=ei-medium

# Usage: sbatch PATH/TO/rds_to_h5ad.sh <SCRIPT_PATH> <RDS_FILE_PATH>

# Arguments
SCRIPT_PATH="$1"
RDS_FILE_PATH="$2"

source /hpc-home/yep25yan/mamba/etc/profile.d/conda.sh
conda activate r_env

script_dir=$(dirname "${SCRIPT_PATH}")
script_name=$(basename "${SCRIPT_PATH}")

echo "Running ${script_name}"
Rscript "${SCRIPT_PATH}" "${RDS_FILE_PATH}" \
        1> "${script_dir}/${script_name}.stdout" \
        2> "${script_dir}/${script_name}.stderr"
echo "${script_name} completed"

exit 0