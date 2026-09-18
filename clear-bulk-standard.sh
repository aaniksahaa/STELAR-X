#!/usr/bin/env bash
# Remove generated standard-dataset statistics without touching source data or
# results outside the selected method(s).

set -euo pipefail

BASE_DIR="${HOME}/phylogeny"
DATASET_DIR=""
METHOD=""
DRY_RUN=false
ASSUME_YES=false
ALL_RESULTS=false
SUPPORTED_METHODS=(stelarx stelar aster astral treeqmc wqfmtree supertriplets stp-nni tmc)

show_usage() {
  cat <<'EOF'
Usage: ./clear-bulk-standard.sh --method METHOD [options]

Remove generated files for one method, or all methods, from the bulk-standard
dataset tree.
By default, only method-specific statistics CSVs and lock markers are removed;
output trees and logs are preserved. This is enough to exclude the method from
the next collect-stats-standard.sh run and allow run-bulk-standard.sh to rerun it.

Required:
  --method, -m METHOD   stelarx | stelar | aster | astral | treeqmc |
                        wqfmtree | supertriplets | stp-nni | tmc | all

Paths:
  --base-dir, -b DIR    Base directory (default: $HOME/phylogeny)
  --dataset-dir, -d DIR Standard dataset directory
                        (default: BASE_DIR/datasets/standard)

Modes:
  --dry-run             List matching targets without deleting anything
  --yes, -y             Delete without an interactive confirmation
  --all-results         Remove the selected method(s)' complete *_outputs
                        directories, including trees, CSVs, logs, and locks
  --help, -h            Show this help

Examples:
  ./clear-bulk-standard.sh --method stelarx --dry-run
  ./clear-bulk-standard.sh --method stelar --yes
  ./clear-bulk-standard.sh --method all --dry-run
  ./clear-bulk-standard.sh --method stelarx
  ./clear-bulk-standard.sh --method stelarx --yes
EOF
}

normalize_method() {
  case "${1,,}" in
    stelarx)                      printf 'stelarx' ;;
    stelar|stelar-x)              printf 'stelar' ;;
    aster)                        printf 'aster' ;;
    astral)                       printf 'astral' ;;
    treeqmc|tree-qmc)             printf 'treeqmc' ;;
    wqfm|wqfmtree|wqfm-tree)      printf 'wqfmtree' ;;
    supertriplets|super-triplets) printf 'supertriplets' ;;
    stp-nni|stpnni)               printf 'stp-nni' ;;
    tmc)                          printf 'tmc' ;;
    all)                          printf 'all' ;;
    *) return 1 ;;
  esac
}

is_supported_output_dir_name() {
  local name="$1" method
  for method in "${SUPPORTED_METHODS[@]}"; do
    [[ "$name" == "${method}_outputs" ]] && return 0
  done
  return 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --method|-m)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 2; }
      METHOD="$2"
      shift 2
      ;;
    --method=*) METHOD="${1#*=}"; shift ;;
    --base-dir|-b)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 2; }
      BASE_DIR="$2"
      shift 2
      ;;
    --base-dir=*) BASE_DIR="${1#*=}"; shift ;;
    --dataset-dir|-d)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 2; }
      DATASET_DIR="$2"
      shift 2
      ;;
    --dataset-dir=*) DATASET_DIR="${1#*=}"; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    --yes|-y) ASSUME_YES=true; shift ;;
    --all-results) ALL_RESULTS=true; shift ;;
    --help|-h) show_usage; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; show_usage >&2; exit 2 ;;
  esac
done

if [[ -z "$METHOD" ]]; then
  echo "Error: --method is required." >&2
  show_usage >&2
  exit 2
fi

METHOD_INPUT="$METHOD"
if ! METHOD="$(normalize_method "$METHOD_INPUT")"; then
  echo "Error: unsupported method '$METHOD_INPUT'." >&2
  exit 2
fi

if [[ -z "$DATASET_DIR" ]]; then
  DATASET_DIR="${BASE_DIR%/}/datasets/standard"
fi

if [[ ! -d "$DATASET_DIR" ]]; then
  echo "Error: dataset directory does not exist: $DATASET_DIR" >&2
  exit 2
fi

DATASET_DIR="$(realpath "$DATASET_DIR")"
case "$DATASET_DIR" in
  /|"$HOME")
    echo "Error: refusing to clean unsafe dataset directory: $DATASET_DIR" >&2
    exit 2
    ;;
esac

declare -a METHODS=()
if [[ "$METHOD" == "all" ]]; then
  METHODS=("${SUPPORTED_METHODS[@]}")
else
  METHODS=("$METHOD")
fi

declare -a OUTPUT_DIRS=()
declare -a TARGETS=()
for selected_method in "${METHODS[@]}"; do
  mapfile -d '' -O "${#OUTPUT_DIRS[@]}" OUTPUT_DIRS < <(
    find "$DATASET_DIR" -type d -name "${selected_method}_outputs" -print0 2>/dev/null
  )

  if [[ "$ALL_RESULTS" == true ]]; then
    continue
  fi

  for output_dir in "${OUTPUT_DIRS[@]}"; do
    [[ "$(basename "$output_dir")" == "${selected_method}_outputs" ]] || continue
    mapfile -d '' -O "${#TARGETS[@]}" TARGETS < <(
      find "$output_dir" -type f \
        \( -name "stat-${selected_method}.csv" \
           -o -name "*-${selected_method}_stats.csv" \
           -o -name ".${selected_method}.lock" \) \
        -print0 2>/dev/null
    )
  done
done

if [[ "$ALL_RESULTS" == true ]]; then
  TARGETS=("${OUTPUT_DIRS[@]}")
fi

echo "Dataset: $DATASET_DIR"
echo "Method:  $METHOD"
if [[ "$ALL_RESULTS" == true ]]; then
  echo "Mode:    complete method output directories"
else
  echo "Mode:    statistics CSVs and lock markers (trees/logs preserved)"
fi
echo

if [[ ${#TARGETS[@]} -eq 0 ]]; then
  if [[ "$METHOD" == "all" ]]; then
    echo "No matching results for any supported method. Nothing to remove."
  else
    echo "No matching $METHOD results found. Nothing to remove."
  fi
  exit 0
fi

echo "Targets (${#TARGETS[@]}):"
for target in "${TARGETS[@]}"; do
  case "$target" in
    "$DATASET_DIR"/*) ;;
    *) echo "Error: unsafe target escaped dataset directory: $target" >&2; exit 3 ;;
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
  read -r -p "Remove these ${#TARGETS[@]} target(s)? [y/N] " reply
  case "${reply,,}" in
    y|yes) ;;
    *) echo "Cancelled; nothing was removed."; exit 0 ;;
  esac
fi

if [[ "$ALL_RESULTS" == true ]]; then
  for target in "${TARGETS[@]}"; do
    is_supported_output_dir_name "$(basename "$target")" || {
      echo "Error: refusing unexpected directory target: $target" >&2
      exit 3
    }
    rm -rf -- "$target"
  done
else
  for target in "${TARGETS[@]}"; do
    rm -f -- "$target"
  done
fi

if [[ "$METHOD" == "all" ]]; then
  echo "Removed ${#TARGETS[@]} target(s) across all supported methods."
  echo "Run collect-stats-standard.sh again to rebuild the merged CSV without the cleared rows."
else
  echo "Removed ${#TARGETS[@]} $METHOD target(s)."
  echo "Run collect-stats-standard.sh again to rebuild the merged CSV without $METHOD rows."
fi
