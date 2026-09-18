#!/usr/bin/env bash
# Back-fill or refresh the reproducibility mirror of simulated-run outputs.
#
# Every results directory found in the SimPhy data tree,
#   <data>/<dataset>/<replicate>/<method>_outputs/<setting>
# is copied to
#   <outputs>/<method>_outputs/<dataset>/<replicate>/<setting>
# together with the dataset's SimPhy <dataset>.command/.params files. Gene trees,
# true species trees, and SimPhy databases are never copied.
#
# New runs mirror themselves automatically (test-stelarx-simulated.sh); this tool
# exists for results produced before the mirror existed and for verifying that
# the mirror is complete before uploading it.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/phylogeny-data-dir.sh"
source "${SCRIPT_DIR}/scripts/simphy-outputs-dir.sh"

DATA_DIR=""
OUTPUTS_DIR=""
METHODS_RAW=""
DRY_RUN=false
QUIET=false

print_help() {
  cat <<EOF
sync-simulated-outputs.sh

Copies every <dataset>/<replicate>/<method>_outputs/<setting> results directory
from the SimPhy data tree into the outputs mirror, alongside each dataset's
SimPhy .command/.params files. Existing mirror leaves are replaced so they equal
the current results exactly. Simulated gene trees, true species trees, and
SimPhy databases are never copied.

Options:
  --simphy-data-dir PATH     Data tree to read (default: \$PHYLOGENY_DATA_DIR/simphy/data)
  --simphy-outputs-dir PATH  Mirror root to write (default: \$PHYLOGENY_DATA_DIR/outputs/simphy,
                             i.e. ".../simphy/data" mirrors into ".../outputs/simphy")
  --methods LIST             Only mirror these methods, comma/space separated
                             (e.g. "stelarx" or "stelarx_outputs"; default: all)
  --dry-run                  Show what would be mirrored without writing
  --quiet, -q                Print only the summary and problems
  --help, -h                 Show this message

Examples:
  ./sync-simulated-outputs.sh --dry-run
  ./sync-simulated-outputs.sh
  ./sync-simulated-outputs.sh --methods stelarx
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --simphy-data-dir|--data-dir|--simphy-outputs-dir|--outputs-dir|--methods)
      if [[ $# -lt 2 ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 2
      fi
      ;;
  esac
  case "$1" in
    --simphy-data-dir|--data-dir) DATA_DIR="$2"; shift 2 ;;
    --simphy-data-dir=*|--data-dir=*) DATA_DIR="${1#*=}"; shift ;;
    --simphy-outputs-dir|--outputs-dir) OUTPUTS_DIR="$2"; shift 2 ;;
    --simphy-outputs-dir=*|--outputs-dir=*) OUTPUTS_DIR="${1#*=}"; shift ;;
    --methods) METHODS_RAW="$2"; shift 2 ;;
    --methods=*) METHODS_RAW="${1#*=}"; shift ;;
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

DATA_DIR="$(stelarx_prepare_simphy_data_dir "$DATA_DIR")" || exit 2
OUTPUTS_DIR="$(stelarx_prepare_simphy_outputs_dir "$OUTPUTS_DIR" "$DATA_DIR")" || exit 2

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

echo "STELAR-X simulated outputs mirror sync"
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
datasets_without_command=()
declare -A seen_dataset=()

while IFS= read -r -d '' results_dir; do
  mapfile -t parts < <(stelarx_simphy_results_components "$DATA_DIR" "$results_dir" 2>/dev/null)
  if [[ ${#parts[@]} -ne 4 ]]; then
    echo "  Skip (unexpected layout): $results_dir" >&2
    continue
  fi
  dataset="${parts[0]}"; replicate="${parts[1]}"; method_dir="${parts[2]}"; setting="${parts[3]}"

  if [[ ${#METHOD_FILTER[@]} -gt 0 && -z "${METHOD_FILTER[$method_dir]:-}" ]]; then
    ((skipped_filter++)) || true
    continue
  fi

  target="${OUTPUTS_DIR}/${method_dir}/${dataset}/${replicate}/${setting}"
  if [[ "$DRY_RUN" == true ]]; then
    [[ "$QUIET" == true ]] || echo "  Would mirror: ${dataset}/${replicate}/${method_dir}/${setting} -> ${target}"
    forbidden="$(stelarx_simphy_find_forbidden_in_mirror "$results_dir" | head -n1)"
    if [[ -n "$forbidden" ]]; then
      echo "  Error: results directory contains simulated input data and would be refused: $forbidden" >&2
      ((failed++)) || true
    else
      ((mirrored++)) || true
    fi
    if [[ -z "${seen_dataset[${method_dir}/${dataset}]:-}" ]]; then
      seen_dataset["${method_dir}/${dataset}"]=1
      command_source="${DATA_DIR}/${dataset}/${dataset}.command"
      if [[ ! -f "$command_source" ]]; then
        if [[ "$dataset" == *_incomplete && -f "${DATA_DIR}/${dataset%_incomplete}/${dataset%_incomplete}.command" ]]; then
          :
        else
          datasets_without_command+=("${method_dir}/${dataset}")
        fi
      fi
    fi
    continue
  fi

  if mirrored_path="$(stelarx_mirror_simulated_results "$DATA_DIR" "$OUTPUTS_DIR" "$results_dir")"; then
    ((mirrored++)) || true
    [[ "$QUIET" == true ]] || echo "  Mirrored: ${dataset}/${replicate}/${method_dir}/${setting} -> ${mirrored_path}"
  else
    ((failed++)) || true
    echo "  Error: could not mirror $results_dir" >&2
  fi

  if [[ -z "${seen_dataset[${method_dir}/${dataset}]:-}" ]]; then
    seen_dataset["${method_dir}/${dataset}"]=1
    if [[ ! -f "${OUTPUTS_DIR}/${method_dir}/${dataset}/${dataset}.command" ]] &&
       ! { [[ "$dataset" == *_incomplete ]] && [[ -f "${OUTPUTS_DIR}/${method_dir}/${dataset}/${dataset%_incomplete}.command" ]]; }; then
      datasets_without_command+=("${method_dir}/${dataset}")
    fi
  fi
done < <(stelarx_list_simulated_results_dirs "$DATA_DIR")

echo
if [[ "$DRY_RUN" == true ]]; then
  echo "Summary (dry run): would mirror=$mirrored filtered-out=$skipped_filter problems=$failed"
else
  echo "Summary: mirrored=$mirrored filtered-out=$skipped_filter failed=$failed"
fi
if [[ ${#datasets_without_command[@]} -gt 0 ]]; then
  echo "Warning: ${#datasets_without_command[@]} dataset(s) have no SimPhy .command file in the data tree:" >&2
  for entry in "${datasets_without_command[@]}"; do
    echo "  $entry" >&2
  done
  echo "  Their mirrors lack the simulation command needed for full reproducibility." >&2
fi

if (( failed > 0 )); then
  exit 1
fi
