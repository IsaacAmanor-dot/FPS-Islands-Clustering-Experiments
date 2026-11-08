#!/bin/bash
#SBATCH --job-name=SI_grids
#SBATCH --partition=rn-long-40core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --output=comparison/grid_generation.out
#SBATCH --error=comparison/grid_generation.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

if [[ -f inputs/whole_grid/whole.nrg && -f multigrid_323/components.json && -f multigrid_20/components.json ]]; then
    python3 -u 004.generate_multigrids.py validate
else
    python3 -u 004.generate_multigrids.py
fi
