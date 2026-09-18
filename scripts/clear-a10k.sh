#!/usr/bin/env bash
# Remove all STELAR-X results produced by run-a10k.sh while preserving the A10K
# source data, rooted gene trees, and simulation files.

set -euo pipefail

DATA_DIR=""
DRY_RUN=false
ASSUME_YES=false

show_usage() {
  cat <<'EOF'
Usage: ./scripts/clear-a10k.sh --data-dir DIR [options]

Remove all STELAR-X A10K results beneath:
  DIR/10k-simphy/R*/stelarx_outputs

The merged DIR/a10k_stelarx_scores_merged.csv file is also removed when present.
Input gene trees, rooted gene trees, species trees, and all other dataset files
are preserved.

The reproducibility mirror (by default DIR/../outputs/<dataset>) is NOT touched,
so results wiped here survive there. Run ./scripts/sync-a10k-outputs.sh before clearing
if you want the mirror to be complete first.

Required:
  --data-dir DIR   A10K dataset root containing 10k-simphy/

Options:
  --dry-run        List exact targets without deleting anything
  --yes, -y        Delete without interactive confirmation
  --help, -h       Show this help

Examples:
  ./scripts/clear-a10k.sh --data-dir /path/to/10k-astral-dataset --dry-run
  ./scripts/clear-a10k.sh --data-dir /path/to/10k-astral-dataset --yes

Do not run this cleaner while run-a10k.sh is active on the same dataset.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 2; }
      DATA_DIR="$2"
      shift 2
      ;;
    --data-dir=*) DATA_DIR="${1#*=}"; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    --yes|-y) ASSUME_YES=true; shift ;;
    --help|-h) show_usage; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; show_usage >&2; exit 2 ;;
  esac
done

if [[ -z "$DATA_DIR" ]]; then
  echo "Error: --data-dir is required." >&2
  show_usage >&2
  exit 2
fi

if [[ ! -d "$DATA_DIR" ]]; then
  echo "Error: data directory does not exist: $DATA_DIR" >&2
  exit 2
fi

DATA_DIR="$(realpath "$DATA_DIR")"
case "$DATA_DIR" in
  /|"$HOME")
    echo "Error: refusing to clean unsafe data directory: $DATA_DIR" >&2
    exit 2
    ;;
esac

SIMPHY_DIR="${DATA_DIR}/10k-simphy"
if [[ ! -d "$SIMPHY_DIR" ]]; then
  echo "Error: expected A10K dataset directory at $SIMPHY_DIR" >&2
  exit 2
fi

declare -a RESULT_DIRS=()
declare -a TARGETS=()
while IFS= read -r -d '' replicate_dir; do
  result_dir="${replicate_dir}/stelarx_outputs"
  [[ -d "$result_dir" ]] && RESULT_DIRS+=("$result_dir")
done < <(find "$SIMPHY_DIR" -mindepth 1 -maxdepth 1 -type d -name 'R*' -print0 | sort -z -V)

TARGETS=("${RESULT_DIRS[@]}")
MERGED_CSV="${DATA_DIR}/a10k_stelarx_scores_merged.csv"
[[ -f "$MERGED_CSV" ]] && TARGETS+=("$MERGED_CSV")

echo "A10K data: $DATA_DIR"
echo "Results:   $SIMPHY_DIR/R*/stelarx_outputs"
echo

if [[ ${#TARGETS[@]} -eq 0 ]]; then
  echo "No A10K STELAR-X results found. Nothing to remove."
  exit 0
fi

echo "Targets (${#TARGETS[@]}):"
for target in "${TARGETS[@]}"; do
  case "$target" in
    "$SIMPHY_DIR"/R*/stelarx_outputs|"$MERGED_CSV") ;;
    *) echo "Error: unsafe target outside the A10K result layout: $target" >&2; exit 3 ;;
  esac
  printf '  %s\n' "$target"
done
echo

if [[ "$DRY_RUN" == true ]]; then
  echo "Dry run only; nothing was removed."
  exit 0
fi

if [[ "$ASSUME_YES" != true ]]; then
  if [[ ! -t 0 ]]; then
    echo "Refusing non-interactive deletion without --yes. Use --dry-run to preview." >&2
    exit 4
  fi
  read -r -p "Remove all listed A10K result targets? [y/N] " reply
  case "${reply,,}" in
    y|yes) ;;
    *) echo "Cancelled; nothing was removed."; exit 0 ;;
  esac
fi

for result_dir in "${RESULT_DIRS[@]}"; do
  [[ "$(basename "$result_dir")" == "stelarx_outputs" ]] || {
    echo "Error: refusing unexpected result directory: $result_dir" >&2
    exit 3
  }
  [[ "$(basename "$(dirname "$result_dir")")" == R* ]] || {
    echo "Error: result directory is not beneath an R* replicate: $result_dir" >&2
    exit 3
  }
  rm -rf -- "$result_dir"
done

if [[ -f "$MERGED_CSV" ]]; then
  rm -f -- "$MERGED_CSV"
fi

echo "Removed ${#RESULT_DIRS[@]} replicate result director$(
  [[ ${#RESULT_DIRS[@]} -eq 1 ]] && printf 'y' || printf 'ies'
) and the merged CSV if it existed."
echo "All A10K input and simulation files were preserved."
echo "The outputs mirror was not touched; it still holds these results."
