#!/bin/bash
#SBATCH --job-name=cite-seq-wnn-trajectory
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
# Edit for your site: partition/account, and `module load` names.
##SBATCH --partition=compute
##SBATCH --account=YOUR_ACCOUNT
#
# Submit from the repository root:  mkdir -p logs && sbatch hpc/submit_slurm.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Site modules (names vary; failures are ignored if the tool is already on PATH).
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

export RUNTIME=apptainer
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SCRATCH:-$HOME}/.apptainer}"

# Images are Apptainer .sif files (build once, see HPC_AWS.md). Data should be
# fetched beforehand on a login node:  RUNTIME=apptainer ./run_analysis.sh prefetch
./run_analysis.sh run
