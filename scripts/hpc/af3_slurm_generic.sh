#!/usr/bin/env bash
# Reference Slurm job body for a site-managed AlphaFold 3 wrapper.
# NOT used automatically by ppigFinder; administrator must adapt sbatch directives.
# Do not submit this template without setting the environment and paths.
set -Eeuo pipefail

: "${PPIG_AF3_JSON:?Set the absolute cluster path to a local AF3 JSON input}"
: "${PPIG_AF3_WORKDIR:?Set the absolute cluster path to the output directory}"
: "${PPIG_AF3_JOB_NAME:?Set a safe job name}"
AF3_EXECUTABLE="${PPIG_AF3_EXECUTABLE:-af3}"

if ! command -v "$AF3_EXECUTABLE" >/dev/null 2>&1; then
  echo "ERROR: AF3 executable '$AF3_EXECUTABLE' unavailable. Activate the site's approved environment." >&2
  exit 127
fi

mkdir -p "$PPIG_AF3_WORKDIR"
args=(
  --json-path "$PPIG_AF3_JSON"
  --job-name "$PPIG_AF3_JOB_NAME"
  --workdir "$PPIG_AF3_WORKDIR"
  --stage all
  --executor local
)
[[ -z "${AF3_IMAGE:-}" ]] || args+=(--image "$AF3_IMAGE")
[[ -z "${AF3_MODEL_DIR:-}" ]] || args+=(--model-dir "$AF3_MODEL_DIR")
[[ -z "${AF3_DB_DIR:-}" ]] || args+=(--db-dir "$AF3_DB_DIR")

# Do not add --slurm-* flags: this script itself runs inside a Slurm allocation.
"$AF3_EXECUTABLE" "${args[@]}"
