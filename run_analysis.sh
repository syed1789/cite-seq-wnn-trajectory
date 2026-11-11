#!/usr/bin/env bash
# Run the whole analysis in containers. Works with Docker (laptop, EC2) or
# Apptainer (HPC, where Docker is not permitted).
#
#   ./run_analysis.sh prefetch   # download input data (needs internet)
#   ./run_analysis.sh run        # run the analysis (default)
#   RUNTIME=apptainer ./run_analysis.sh run
#   DRY_RUN=1 ./run_analysis.sh  # print the commands without running them
set -euo pipefail
cd "$(dirname "$0")"
RUNTIME="${RUNTIME:-docker}"
DRY_RUN="${DRY_RUN:-0}"

run() { if [ "$DRY_RUN" = 1 ]; then echo "+ $*"; else "$@"; fi; }

# build_image <docker-tag> <build-context>
#   docker:    builds the image.
#   apptainer: expects ./<tag>.sif (build once: see HPC_AWS.md); builds it here
#              from the local Docker daemon if Docker happens to be available.
build_image() {
  local tag="$1" ctx="$2" sif="${1%%:*}.sif"
  case "$RUNTIME" in
    docker) run docker build -t "$tag" "$ctx" ;;
    apptainer)
      if [ ! -f "$sif" ]; then
        if command -v docker >/dev/null 2>&1; then
          run docker build -t "$tag" "$ctx"
          run apptainer build "$sif" "docker-daemon://$tag"
        else
          echo "Missing $sif. Build it on a machine with Docker (see HPC_AWS.md)." >&2
          [ "$DRY_RUN" = 1 ] || exit 1
        fi
      fi ;;
    *) echo "RUNTIME must be docker or apptainer" >&2; exit 1 ;;
  esac
}

# ctr <docker-ref-or-tag> <apptainer-source> <command...>
#   apptainer-source: a .sif path, or docker://<ref> for public registry images.
ctr() {
  local ref="$1" src="$2"; shift 2
  case "$RUNTIME" in
    docker) run docker run --rm -e HOME=/tmp -e MPLCONFIGDIR=/tmp -e NUMBA_CACHE_DIR=/tmp \
              -v "$PWD":/work -w /work "$ref" "$@" ;;
    apptainer) run apptainer exec --bind "$PWD":/work --pwd /work "$src" "$@" ;;
  esac
}

step="${1:-run}"

WNN=cite-seq-wnn:latest
TRAJ=monocle3-trajectory:latest
H5=data/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5
fetch() {
  mkdir -p data
  [ -f "$H5" ] || run curl -L --fail -o "$H5" \
    https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5
}
case "$step" in
  prefetch) fetch ;;
  run)
    fetch
    build_image "$WNN" wnn/
    ctr "$WNN" cite-seq-wnn.sif Rscript /work/wnn/run_wnn.R
    build_image "$TRAJ" trajectory/
    ctr "$TRAJ" monocle3-trajectory.sif Rscript /work/trajectory/run_trajectory.R ;;
  *) echo "usage: $0 [prefetch|run]" >&2; exit 1 ;;
esac
