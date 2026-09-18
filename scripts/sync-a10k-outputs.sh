#!/usr/bin/env bash
# Back-fill or refresh the reproducibility mirror of A10K run outputs.
#
# Every results directory found in the A10K data tree,
#   <data>/10k-simphy/<replicate>/<method>_outputs/<tree-type>/<setting>
# is copied to
#   <outputs>/<method>_outputs/<replicate>/<tree-type>/<setting>
# together with the dataset provenance record and the per-replicate input
# fingerprints. Gene trees and true species trees are never copied.
#
# New runs mirror themselves automatically (run-a10k.sh); this tool exists for
# results produced before the mirror existed and for verifying that the mirror
# is complete before uploading it.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/a10k-outputs-dir.sh"

DATA_DIR=""
OUTPUTS_DIR=""
METHODS_RAW=""
DRY_RUN=false
QUIET=false

print_help() {
  cat <<EOF
sync-a10k-outputs.sh

Copies every 10k-simphy/<replicate>/<method>_outputs/<tree-type>/<setting>
results directory from the A10K data tree into the outputs mirror, alongside the
dataset record and each replicate's input fingerprints. Existing mirror leaves
are replaced so they equal the current results exactly. A10K gene trees and
species trees are never copied.

Options:
  --data-dir PATH          A10K dataset root containing 10k-simphy/ (required
                           unless \$PHYLOGENY_DATA_DIR/10k-astral-dataset exists)
  --a10k-outputs-dir PATH  Mirror root to write (default: the "outputs" sibling
                           of the dataset, i.e. ".../10k-astral-dataset"
                           mirrors into ".../outputs/10k-astral-dataset")
  --methods LIST           Only mirror these methods, comma/space separated
                           (e.g. "stelarx" or "stelarx_outputs"; default: all)
  --dry-run                Show what would be mirrored without writing
  --quiet, -q              Print only the summary and problems
  --help, -h               Show this message

Examples:
  ./scripts/sync-a10k-outputs.sh --dry-run
  ./scripts/sync-a10k-outputs.sh --data-dir \$PHYLOGENY_DATA_DIR/10k-astral-dataset
  ./scripts/sync-a10k-outputs.sh --methods stelarx
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir|--a10k-data-dir|--a10k-outputs-dir|--outputs-dir|--methods|--method)
      if [[ $# -lt 2 ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 2
      fi
      ;;
  esac
  case "$1" in
    --data-dir|--a10k-data-dir) DATA_DIR="$2"; shift 2 ;;
    --data-dir=*|--a10k-data-dir=*) DATA_DIR="${1#*=}"; shift ;;
    --a10k-outputs-dir|--outputs-dir) OUTPUTS_DIR="$2"; shift 2 ;;
    --a10k-outputs-dir=*|--outputs-dir=*) OUTPUTS_DIR="${1#*=}"; shift ;;
    --methods|--method) METHODS_RAW="$2"; shift 2 ;;
    --methods=*|--method=*) METHODS_RAW="${1#*=}"; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    --quiet|-q) QUIET=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *)
      echo "Error: unknown option '$1'." >&2
      print_help >&2
      exit 2
      ;;
  esac
done

# Resolve the dataset root: an explicit --data-dir wins, otherwise the standard
# location beneath PHYLOGENY_DATA_DIR.
if [[ -z "$DATA_DIR" ]]; then
  if [[ -z "${PHYLOGENY_DATA_DIR:-}" ]]; then
    echo "Error: --data-dir is required (or set PHYLOGENY_DATA_DIR)." >&2
    exit 2
  fi
  DATA_DIR="${PHYLOGENY_DATA_DIR%/}/10k-astral-dataset"
fi
if [[ "$DATA_DIR" == "~/"* ]]; then
  DATA_DIR="${HOME}/${DATA_DIR:2}"
fi
if [[ ! -d "$DATA_DIR" ]]; then
  echo "Error: A10K data directory does not exist: $DATA_DIR" >&2
  exit 2
fi
DATA_DIR="$(cd "$DATA_DIR" && pwd -P)"
if [[ ! -d "${DATA_DIR}/${STELARX_A10K_REPLICATE_ROOT}" ]]; then
  echo "Error: expected ${STELARX_A10K_REPLICATE_ROOT} at ${DATA_DIR}/${STELARX_A10K_REPLICATE_ROOT}" >&2
  exit 2
fi

OUTPUTS_DIR="$(stelarx_prepare_a10k_outputs_dir "$OUTPUTS_DIR" "$DATA_DIR")" || exit 2

# Normalize the method filter to "<method>_outputs" directory names.
declare -A METHOD_FILTER=()
if [[ -n "$METHODS_RAW" ]]; then
  read -r -a method_items <<< "${METHODS_RAW//,/ }"
  for method in "${method_items[@]}"; do
    [[ -z "$method" ]] && continue
    [[ "$method" == *_outputs ]] || method="${method}_outputs"
    METHOD_FILTER["$method"]=1
  done
  if [[ ${#METHOD_FILTER[@]} -eq 0 ]]; then
    echo "Error: --methods did not name any method." >&2
    exit 2
  fi
fi

echo "STELAR-X A10K outputs mirror sync"
echo "Data directory:    $DATA_DIR"
echo "Outputs directory: $OUTPUTS_DIR"
if [[ ${#METHOD_FILTER[@]} -gt 0 ]]; then
  echo "Methods:           ${!METHOD_FILTER[*]}"
fi
[[ "$DRY_RUN" == true ]] && echo "Dry run:           yes"
echo

mirrored=0
skipped_filter=0
failed=0
declare -A seen_method=()

while IFS= read -r -d '' results_dir; do
  mapfile -t parts < <(stelarx_a10k_results_components "$DATA_DIR" "$results_dir" 2>/dev/null)
  if [[ ${#parts[@]} -ne 4 ]]; then
    echo "  Skip (unexpected layout): $results_dir" >&2
    continue
  fi
  replicate="${parts[0]}"; method_dir="${parts[1]}"; tree_type="${parts[2]}"; setting="${parts[3]}"

  if [[ ${#METHOD_FILTER[@]} -gt 0 && -z "${METHOD_FILTER[$method_dir]:-}" ]]; then
    ((skipped_filter++)) || true
    continue
  fi
  seen_method["$method_dir"]=1

  target="${OUTPUTS_DIR}/${method_dir}/${replicate}/${tree_type}/${setting}"
  if [[ "$DRY_RUN" == true ]]; then
    [[ "$QUIET" == true ]] || echo "  Would mirror: ${replicate}/${method_dir}/${tree_type}/${setting} -> ${target}"
    forbidden="$(stelarx_a10k_find_forbidden_in_mirror "$results_dir" | head -n1)"
    if [[ -n "$forbidden" ]]; then
      echo "  Error: results directory contains A10K input data and would be refused: $forbidden" >&2
      ((failed++)) || true
    else
      ((mirrored++)) || true
    fi
    continue
  fi

  if mirrored_path="$(stelarx_mirror_a10k_results "$DATA_DIR" "$OUTPUTS_DIR" "$results_dir")"; then
    ((mirrored++)) || true
    [[ "$QUIET" == true ]] || echo "  Mirrored: ${replicate}/${method_dir}/${tree_type}/${setting} -> ${mirrored_path}"
  else
    ((failed++)) || true
    echo "  Error: could not mirror $results_dir" >&2
  fi
done < <(stelarx_list_a10k_results_dirs "$DATA_DIR")

# The dataset record is what makes the mirror reproducible; write it for every
# method that has results, including on a dry run's behalf (reported, not written).
if [[ "$DRY_RUN" == false ]]; then
  for method_dir in "${!seen_method[@]}"; do
    if ! stelarx_write_a10k_dataset_record "$DATA_DIR" "${OUTPUTS_DIR}/${method_dir}"; then
      echo "  Error: could not write the dataset record for ${method_dir}" >&2
      ((failed++)) || true
    fi
  done
fi

echo
if [[ "$DRY_RUN" == true ]]; then
  echo "Summary (dry run): would mirror=$mirrored filtered-out=$skipped_filter problems=$failed"
else
  echo "Summary: mirrored=$mirrored filtered-out=$skipped_filter failed=$failed"
fi

if (( failed > 0 )); then
  exit 1
fi
