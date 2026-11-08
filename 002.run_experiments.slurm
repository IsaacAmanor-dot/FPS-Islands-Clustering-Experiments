#!/bin/bash
#SBATCH --job-name=SI_FPS
#SBATCH --partition=rn-long-40core
#SBATCH --array=0-2%3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --time=2-00:00:00
#SBATCH --output=comparison/experiment_%a.out
#SBATCH --error=comparison/experiment_%a.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source 000.config.sh
METHODS=(cartesian multigrid_323 multigrid_20)

# One method per node; concurrent cutoffs are controlled by MAX_PARALLEL.
python3 -u 001.prepare_experiments.py run "${METHODS[SLURM_ARRAY_TASK_ID]}"
