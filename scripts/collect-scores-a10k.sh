#!/usr/bin/env bash
# collect-scores-a10k.sh
# Aggregates per-replicate STELAR-X run stats from the A10K layout.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/a10k-outputs-dir.sh"

DATA_DIR=""
START_REP=""
END_REP=""
A10K_OUTPUTS_DIR=""
OUTPUTS_MIRROR=true

print_help() {
  cat <<EOF
collect-scores-a10k.sh

Usage: $0 --data-dir <dir> --start-rep <N> --end-rep <M> [options]

Writes <data-dir>/a10k_stelarx_scores_merged.csv and copies it into the
reproducibility mirror at <outputs>/stelarx_outputs/.

Options:
  --a10k-outputs-dir PATH  Mirror root (default: the "outputs" sibling of the
                           dataset, i.e. ".../10k-astral-dataset" ->
                           ".../outputs/10k-astral-dataset")
  --no-outputs-mirror      Do not copy the merged CSV into the mirror
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --start-rep|-sr) START_REP="$2"; shift 2 ;;
    --end-rep|-er) END_REP="$2"; shift 2 ;;
    --a10k-outputs-dir|--outputs-dir) A10K_OUTPUTS_DIR="$2"; shift 2 ;;
    --a10k-outputs-dir=*|--outputs-dir=*) A10K_OUTPUTS_DIR="${1#*=}"; shift ;;
    --no-outputs-mirror) OUTPUTS_MIRROR=false; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ -z "$DATA_DIR" || -z "$START_REP" || -z "$END_REP" ]]; then
  echo "Error: --data-dir, --start-rep and --end-rep are required."
  exit 2
fi

DATA_DIR="$(realpath "$DATA_DIR")"
MERGED_CSV="${DATA_DIR}/a10k_stelarx_scores_merged.csv"
echo "alg,setting,replicate,tree_type,rf-rate,optimal-triplet-score,running-time-s,max-cpu-mb,max-gpu-mb" > "$MERGED_CSV"

for i in $(seq "$START_REP" "$END_REP"); do
  while IFS= read -r -d '' stat_file; do
    tail -n +2 "$stat_file" >> "$MERGED_CSV"
  done < <(find "${DATA_DIR}/10k-simphy/R${i}/stelarx_outputs" -type f -name 'stat-stelarx.csv' -print0 2>/dev/null | sort -z)
done

echo "Merged A10K STELAR-X stats saved to: $MERGED_CSV"

# The merged summary belongs with the mirrored results; it is a stat, not input
# data. A mirror problem is reported but never fails the collection itself.
if [[ "$OUTPUTS_MIRROR" == true ]]; then
  if A10K_OUTPUTS_DIR="$(stelarx_prepare_a10k_outputs_dir "$A10K_OUTPUTS_DIR" "$DATA_DIR")"; then
    MIRROR_METHOD_DIR="${A10K_OUTPUTS_DIR}/stelarx_outputs"
    if mkdir -p -- "$MIRROR_METHOD_DIR" &&
       stelarx__copy_small_file_atomic "$MERGED_CSV" "${MIRROR_METHOD_DIR}/$(basename "$MERGED_CSV")"; then
      stelarx_write_a10k_dataset_record "$DATA_DIR" "$MIRROR_METHOD_DIR" || true
      echo "Mirrored merged stats to: ${MIRROR_METHOD_DIR}/$(basename "$MERGED_CSV")"
    else
      echo "WARNING: could not mirror the merged stats into ${MIRROR_METHOD_DIR}" >&2
    fi
  else
    echo "WARNING: outputs mirror was not updated for $MERGED_CSV" >&2
  fi
fi
