#!/bin/bash
#SBATCH -J ftq_sub-001_ses-01_roi95
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Optional (if you want physical cores only; ignore if unsure):
# SBATCH --hint=nomultithread
# Optional (bind ranks/threads to cores for better cache locality):
# SBATCH --cpu-bind=cores

set -euo pipefail

# Ensure log directory exists
mkdir -p logs

# Optional: activate conda env if you use one for R
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate base || true
fi

# Align threaded libs with the Slurm allocation
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

DATA_ROOT="$HOME/Projects/FibrTraQ-data"
SUB="sub-001"
SES="ses-01"
FDESC="roi95"
DISCRET="any"     # any|round|floor|ceil

# Launch under srun so cgroups & CPU/MEM limits are enforced for the task
srun --cpu-bind=cores Rscript --vanilla bench_fibrtraq.R \
  "$DATA_ROOT" "$SUB" "$SES" "$FDESC" "$DISCRET"
