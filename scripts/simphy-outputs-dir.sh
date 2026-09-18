#!/usr/bin/env bash

# Reproducibility mirror for simulated-run outputs.
#
# Every inference run on a SimPhy replicate keeps writing its results inside the
# data tree exactly as before:
#
#   <data>/<dataset>/<replicate>/<method>_outputs/<setting>/...
#
# The functions below additionally maintain a compact, shareable copy that
# contains only the small artifacts (inferred trees, CSVs, run markers/logs) plus
# the SimPhy command files, and never the gene trees or true species trees:
#
#   <outputs>/<method>_outputs/<dataset>/<dataset>.command
#   <outputs>/<method>_outputs/<dataset>/<dataset>.params
#   <outputs>/<method>_outputs/<dataset>/<replicate>/<setting>/...
#
# <outputs> defaults to $PHYLOGENY_DATA_DIR/outputs/simphy for the standard
# $PHYLOGENY_DATA_DIR/simphy/data tree (".../simphy/data" -> ".../outputs/simphy").
# Any other data directory named "data" mirrors into its "outputs" sibling, and
# any other explicit data directory into "<data-dir>_outputs", unless an outputs
# directory is given explicitly.
#
# The A10K dataset has its own mirror with the same guarantees; see
# scripts/a10k-outputs-dir.sh. Primitives shared by both live in
# scripts/outputs-mirror-common.sh.

# shellcheck source=scripts/outputs-mirror-common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/outputs-mirror-common.sh"

# Files that describe how a dataset was simulated. Copied per dataset.
STELARX_SIMPHY_COMMAND_SUFFIXES=(".command" ".params")

# Files that must never appear in the outputs mirror. These are the large
# simulated inputs and SimPhy databases that stay in the data tree only.
STELARX_SIMPHY_FORBIDDEN_MIRROR_NAMES=(
  "all_gt.tre" "s_tree.trees" "l_trees.trees" "g_trees*.trees"
  "*.db" "*.db-journal" "*.zip" "stat-sim.csv"
)

# Print the default outputs mirror root for a resolved SimPhy data directory.
stelarx_default_simphy_outputs_dir() {
  local data_dir="${1%/}"
  if [[ -z "$data_dir" ]]; then
    stelarx__outputs_error "a SimPhy data directory is required to derive the outputs directory."
    return 2
  fi
  local parent_dir
  parent_dir="$(dirname -- "$data_dir")"
  if [[ "$(basename -- "$data_dir")" == "data" && "$(basename -- "$parent_dir")" == "simphy" ]]; then
    # Standard layout: <root>/simphy/data -> <root>/outputs/simphy
    printf '%s/outputs/simphy\n' "$(dirname -- "$parent_dir")"
  elif [[ "$(basename -- "$data_dir")" == "data" ]]; then
    printf '%s/outputs\n' "$parent_dir"
  else
    printf '%s_outputs\n' "$data_dir"
  fi
}

# Resolve, validate, and create the outputs mirror root.
#   stelarx_prepare_simphy_outputs_dir REQUESTED_OUTPUTS_DIR RESOLVED_DATA_DIR
# An empty REQUESTED_OUTPUTS_DIR selects the default described above.
stelarx_prepare_simphy_outputs_dir() {
  local requested_dir="${1:-}" data_dir="${2:-}"

  if [[ -z "$data_dir" ]]; then
    stelarx__outputs_error "internal: resolved SimPhy data directory is required."
    return 2
  fi
  if [[ -z "$requested_dir" ]]; then
    requested_dir="$(stelarx_default_simphy_outputs_dir "$data_dir")" || return 2
  fi
  stelarx__resolve_outputs_root "$requested_dir" "$data_dir" "SimPhy"
}

# True when a file name matches one of the forbidden mirror patterns.
stelarx_simphy_name_is_forbidden_in_mirror() {
  stelarx__name_matches_any "$1" "${STELARX_SIMPHY_FORBIDDEN_MIRROR_NAMES[@]}"
}

# Print every forbidden file found beneath a directory (one per line).
stelarx_simphy_find_forbidden_in_mirror() {
  stelarx__find_names_under "$1" "${STELARX_SIMPHY_FORBIDDEN_MIRROR_NAMES[@]}"
}

# Copy the SimPhy command/params files of a dataset into its mirror directory.
#   stelarx_mirror_simphy_dataset_commands DATA_DIR MIRROR_DATASET_DIR DATASET_NAME
# For "<name>_incomplete" datasets the base dataset's files are copied as well,
# because the incomplete variant is derived from that simulation.
# Returns 0 when at least one .command file was mirrored, 1 otherwise.
stelarx_mirror_simphy_dataset_commands() {
  local data_dir="${1%/}" mirror_dataset_dir="${2%/}" dataset_name="$3"
  local -a names=("$dataset_name")
  local name suffix source destination copied_command=false

  if [[ "$dataset_name" == *_incomplete ]]; then
    names+=("${dataset_name%_incomplete}")
  fi

  mkdir -p -- "$mirror_dataset_dir" || return 1
  for name in "${names[@]}"; do
    for suffix in "${STELARX_SIMPHY_COMMAND_SUFFIXES[@]}"; do
      source="${data_dir}/${name}/${name}${suffix}"
      destination="${mirror_dataset_dir}/${name}${suffix}"
      [[ -f "$source" ]] || continue
      if [[ -f "$destination" ]] && cmp -s -- "$source" "$destination"; then
        [[ "$suffix" == ".command" ]] && copied_command=true
        continue
      fi
      if stelarx__copy_small_file_atomic "$source" "$destination"; then
        [[ "$suffix" == ".command" ]] && copied_command=true
      else
        stelarx__outputs_error "could not copy $source to $destination"
      fi
    done
  done

  [[ "$copied_command" == true ]]
}

# Split a results directory into its mirror components.
#   stelarx_simphy_results_components DATA_DIR RESULTS_DIR
# Prints four lines: dataset, replicate, method-outputs dir name, setting.
# Fails when RESULTS_DIR does not have the expected shape beneath DATA_DIR.
stelarx_simphy_results_components() {
  local data_dir="${1%/}" results_dir="${2%/}"
  local relative dataset replicate method_dir setting rest

  results_dir="$(realpath -m -- "$results_dir")"
  if ! stelarx__path_is_within "$results_dir" "$data_dir" || [[ "$results_dir" == "$data_dir" ]]; then
    stelarx__outputs_error "results directory is not inside the SimPhy data directory: $results_dir"
    return 1
  fi
  relative="${results_dir#"$data_dir"/}"

  dataset="${relative%%/*}"; rest="${relative#*/}"
  [[ "$rest" != "$relative" ]] || { stelarx__outputs_error "unexpected results layout: $relative"; return 1; }
  replicate="${rest%%/*}"; rest="${rest#*/}"
  [[ "$rest" != "$replicate" ]] || { stelarx__outputs_error "unexpected results layout: $relative"; return 1; }
  method_dir="${rest%%/*}"; setting="${rest#*/}"
  if [[ "$setting" == "$method_dir" || -z "$setting" || "$setting" == */* ]]; then
    stelarx__outputs_error "expected <dataset>/<replicate>/<method>_outputs/<setting>, got: $relative"
    return 1
  fi
  if [[ "$method_dir" != *_outputs ]]; then
    stelarx__outputs_error "expected a '<method>_outputs' directory, got '$method_dir' in: $relative"
    return 1
  fi
  for rest in "$dataset" "$replicate" "$method_dir" "$setting"; do
    if [[ "$rest" == . || "$rest" == .. || -z "$rest" ]]; then
      stelarx__outputs_error "unsafe path component in: $relative"
      return 1
    fi
  done

  printf '%s\n%s\n%s\n%s\n' "$dataset" "$replicate" "$method_dir" "$setting"
}

# Mirror one results directory into the outputs tree.
#   stelarx_mirror_simulated_results DATA_DIR OUTPUTS_ROOT RESULTS_DIR
# The mirror leaf is replaced atomically so it always equals the source leaf.
# Dataset command files are refreshed alongside. Prints the mirror leaf path.
stelarx_mirror_simulated_results() {
  local data_dir="${1%/}" outputs_root="${2%/}" results_dir="${3%/}"
  local -a parts=()
  local dataset replicate method_dir setting
  local mirror_dataset_dir mirror_leaf forbidden

  if [[ ! -d "$results_dir" ]]; then
    stelarx__outputs_error "results directory does not exist: $results_dir"
    return 1
  fi
  mapfile -t parts < <(stelarx_simphy_results_components "$data_dir" "$results_dir") || return 1
  [[ ${#parts[@]} -eq 4 ]] || return 1
  dataset="${parts[0]}"; replicate="${parts[1]}"; method_dir="${parts[2]}"; setting="${parts[3]}"

  forbidden="$(stelarx_simphy_find_forbidden_in_mirror "$results_dir" | head -n1)"
  if [[ -n "$forbidden" ]]; then
    stelarx__outputs_error "refusing to mirror a results directory containing simulated input data: $forbidden"
    return 1
  fi

  mirror_dataset_dir="${outputs_root}/${method_dir}/${dataset}"
  mirror_leaf="${mirror_dataset_dir}/${replicate}/${setting}"

  mkdir -p -- "${mirror_dataset_dir}/${replicate}" || return 1
  if ! stelarx_mirror_simphy_dataset_commands "$data_dir" "$mirror_dataset_dir" "$dataset"; then
    echo "Warning: no SimPhy .command file found for dataset '$dataset'; the mirror lacks its simulation command." >&2
  fi

  stelarx__mirror_leaf_atomic "$results_dir" "$mirror_leaf" || return 1

  printf '%s\n' "$mirror_leaf"
}

# Print every results directory beneath a data directory, NUL-separated:
#   <data>/<dataset>/<replicate>/<method>_outputs/<setting>
stelarx_list_simulated_results_dirs() {
  local data_dir="${1%/}"
  find "$data_dir" -mindepth 4 -maxdepth 4 -type d -path '*/*_outputs/*' \
    -not -name '.*' -print0 2>/dev/null | sort -z
}

# Resolve the outputs mirror root without a data directory (for tools that only
# read the mirror). An explicit path wins; otherwise PHYLOGENY_DATA_DIR is
# required and <PHYLOGENY_DATA_DIR>/outputs/simphy is used and created.
stelarx_prepare_simphy_outputs_dir_standalone() {
  stelarx__resolve_outputs_root_standalone "${1:-}" "outputs/simphy" "SimPhy"
}

# Validate a simulated dataset directory name as it appears in the outputs tree.
# Sets BASH_REMATCH: [1]=taxa [2]=gene trees [3]=sb [6]=spmin [7]=spmax
# [8]="_incomplete" or empty.
stelarx_simphy_dataset_name_is_valid() {
  local name="$1"
  [[ "$name" =~ ^t_([1-9][0-9]*)_g_([1-9][0-9]*)_sb_([0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?)_spmin_([1-9][0-9]*)_spmax_([1-9][0-9]*)(_incomplete)?$ ]]
}
