#!/bin/bash
#SBATCH -J ftq_sub-001_ses-01_roi95
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-logs/%x_%j.out
#SBATCH --error=slurm-logs/%x_%j.err
# Optional: physical cores only (no SMT/HT)
#SBATCH --hint=nomultithread

set -euo pipefail

# --- preflight: print environment for provenance ---
echo "[$(date -Is)] Host: $(hostname)"
echo "[$(date -Is)] CWD:  $PWD"
echo "[$(date -Is)] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset}"

# If you use modules instead of conda, do 'module load R' here.
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate base || true
fi

# Align threaded libs to Slurm allocation
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Inputs
DATA_ROOT="$HOME/Projects/FibrTraQ-data"
SUB="sub-001"
SES="ses-01"
FDESC="roi95"
DISCRET="any"     # any|round|floor|ceil

# Ensure log dir exists BEFORE submission (do this in shell: mkdir -p slurm-logs)
# Submit-time directives above will write logs to slurm-logs/

# Launch under srun so Slurm enforces cgroups/CPU binding
srun --cpu-bind=cores Rscript --vanilla bench_fibrtraq.R \
  "${DATA_ROOT}" "${SUB}" "${SES}" "${FDESC}" "${DISCRET}"

echo "[$(date -Is)] Done."
