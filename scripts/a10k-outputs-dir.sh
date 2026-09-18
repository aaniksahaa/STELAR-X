#!/usr/bin/env bash

# Reproducibility mirror for A10K (10k-astral-dataset) run outputs.
#
# Every STELAR-X run on an A10K replicate keeps writing its results inside the
# data tree exactly as before:
#
#   <data>/10k-simphy/<replicate>/<method>_outputs/<tree-type>/<setting>/...
#
# The functions below additionally maintain a compact, shareable copy holding
# only the small artifacts (inferred trees, CSVs, run markers/logs, command
# records) plus the provenance of the dataset itself, and never the gene trees
# or true species trees:
#
#   <outputs>/<method>_outputs/a10k-dataset.command        # dataset provenance
#   <outputs>/<method>_outputs/<replicate>/inputs.tsv      # input fingerprints
#   <outputs>/<method>_outputs/<replicate>/<tree-type>/<setting>/...
#
# <outputs> defaults to the "outputs" sibling of the dataset directory:
# ".../10k-astral-dataset" -> ".../outputs/10k-astral-dataset", matching the
# SimPhy mirror's ".../simphy/data" -> ".../outputs/simphy" rule.
#
# See scripts/simphy-outputs-dir.sh for the simulated-run mirror and
# scripts/outputs-mirror-common.sh for the primitives both share.

# shellcheck source=scripts/outputs-mirror-common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/outputs-mirror-common.sh"

# The replicate tree inside an A10K dataset directory.
STELARX_A10K_REPLICATE_ROOT="10k-simphy"

# The dataset provenance record kept beside each method's mirrored results.
STELARX_A10K_DATASET_RECORD="a10k-dataset.command"

# Per-replicate fingerprint of the inputs that were NOT mirrored.
STELARX_A10K_INPUTS_RECORD="inputs.tsv"

# Gene-tree types run-a10k.sh understands; also the mirror's tree-type level.
STELARX_A10K_TREE_TYPES=("true" "estimated")

# Files that must never appear in the outputs mirror: the A10K gene trees,
# species trees, and any SimPhy artefacts shipped with the dataset.
STELARX_A10K_FORBIDDEN_MIRROR_NAMES=(
  "truegenetrees" "truegenetrees.*" "estimatedgenetrees" "estimatedgenetrees.*"
  "s_tree.trees" "all_gt.tre" "l_trees.trees" "g_trees*.trees"
  "*.db" "*.db-journal" "*.zip" "stat-sim.csv"
)

# Input files of one replicate, relative to the dataset directory. Their sizes
# are recorded in the mirror; their contents never are.
stelarx_a10k_replicate_input_paths() {
  local replicate="$1"
  printf '%s/%s/%s\n' \
    "$STELARX_A10K_REPLICATE_ROOT" "$replicate" "s_tree.trees" \
    "$STELARX_A10K_REPLICATE_ROOT" "$replicate" "truegenetrees" \
    "$STELARX_A10K_REPLICATE_ROOT" "$replicate" "estimatedgenetrees/estimatedgenetrees.tre" \
    "$STELARX_A10K_REPLICATE_ROOT" "$replicate" "estimatedgenetrees/estimatedgenetrees.rooted.tre"
}

# Print the default outputs mirror root for a resolved A10K dataset directory:
#   <parent>/<name>  ->  <parent>/outputs/<name>
stelarx_default_a10k_outputs_dir() {
  local data_dir="${1%/}"
  if [[ -z "$data_dir" ]]; then
    stelarx__outputs_error "an A10K data directory is required to derive the outputs directory."
    return 2
  fi
  local parent_dir name
  parent_dir="$(dirname -- "$data_dir")"
  name="$(basename -- "$data_dir")"
  if [[ "$(basename -- "$parent_dir")" == "outputs" ]]; then
    # The dataset already lives under an "outputs" tree; keep the mirror beside it.
    printf '%s_outputs\n' "$data_dir"
  else
    printf '%s/outputs/%s\n' "$parent_dir" "$name"
  fi
}

# Resolve, validate, and create the outputs mirror root.
#   stelarx_prepare_a10k_outputs_dir REQUESTED_OUTPUTS_DIR RESOLVED_DATA_DIR
# An empty REQUESTED_OUTPUTS_DIR selects the default described above.
stelarx_prepare_a10k_outputs_dir() {
  local requested_dir="${1:-}" data_dir="${2:-}"

  if [[ -z "$data_dir" ]]; then
    stelarx__outputs_error "internal: resolved A10K data directory is required."
    return 2
  fi
  if [[ -z "$requested_dir" ]]; then
    requested_dir="$(stelarx_default_a10k_outputs_dir "$data_dir")" || return 2
  fi
  stelarx__resolve_outputs_root "$requested_dir" "$data_dir" "A10K"
}

# Resolve the mirror root without a data directory (for tools that only read the
# mirror). An explicit path wins; otherwise
# <PHYLOGENY_DATA_DIR>/outputs/10k-astral-dataset is used and created.
stelarx_prepare_a10k_outputs_dir_standalone() {
  stelarx__resolve_outputs_root_standalone "${1:-}" "outputs/10k-astral-dataset" "A10K"
}

# True when a file name matches one of the forbidden mirror patterns.
stelarx_a10k_name_is_forbidden_in_mirror() {
  stelarx__name_matches_any "$1" "${STELARX_A10K_FORBIDDEN_MIRROR_NAMES[@]}"
}

# Print every forbidden file found beneath a directory (one per line).
stelarx_a10k_find_forbidden_in_mirror() {
  stelarx__find_names_under "$1" "${STELARX_A10K_FORBIDDEN_MIRROR_NAMES[@]}"
}

# True for a gene-tree type the A10K layout uses.
stelarx_a10k_tree_type_is_valid() {
  local candidate="$1" tree_type
  for tree_type in "${STELARX_A10K_TREE_TYPES[@]}"; do
    [[ "$candidate" == "$tree_type" ]] && return 0
  done
  return 1
}

# Split a results directory into its mirror components.
#   stelarx_a10k_results_components DATA_DIR RESULTS_DIR
# Prints four lines: replicate, method-outputs dir name, tree type, setting.
# Fails when RESULTS_DIR does not have the expected shape beneath DATA_DIR.
stelarx_a10k_results_components() {
  local data_dir="${1%/}" results_dir="${2%/}"
  local relative replicate method_dir tree_type setting rest component

  results_dir="$(realpath -m -- "$results_dir")"
  if ! stelarx__path_is_within "$results_dir" "$data_dir" || [[ "$results_dir" == "$data_dir" ]]; then
    stelarx__outputs_error "results directory is not inside the A10K data directory: $results_dir"
    return 1
  fi
  relative="${results_dir#"$data_dir"/}"

  rest="$relative"
  component="${rest%%/*}"
  if [[ "$component" != "$STELARX_A10K_REPLICATE_ROOT" || "$rest" == "$component" ]]; then
    stelarx__outputs_error "expected ${STELARX_A10K_REPLICATE_ROOT}/<replicate>/<method>_outputs/<tree-type>/<setting>, got: $relative"
    return 1
  fi
  rest="${rest#*/}"
  replicate="${rest%%/*}"; [[ "$rest" != "$replicate" ]] || { stelarx__outputs_error "unexpected results layout: $relative"; return 1; }
  rest="${rest#*/}"
  method_dir="${rest%%/*}"; [[ "$rest" != "$method_dir" ]] || { stelarx__outputs_error "unexpected results layout: $relative"; return 1; }
  rest="${rest#*/}"
  tree_type="${rest%%/*}"; setting="${rest#*/}"
  if [[ "$setting" == "$tree_type" || -z "$setting" || "$setting" == */* ]]; then
    stelarx__outputs_error "expected ${STELARX_A10K_REPLICATE_ROOT}/<replicate>/<method>_outputs/<tree-type>/<setting>, got: $relative"
    return 1
  fi
  if [[ "$method_dir" != *_outputs ]]; then
    stelarx__outputs_error "expected a '<method>_outputs' directory, got '$method_dir' in: $relative"
    return 1
  fi
  if ! stelarx_a10k_tree_type_is_valid "$tree_type"; then
    stelarx__outputs_error "expected tree type ${STELARX_A10K_TREE_TYPES[*]}, got '$tree_type' in: $relative"
    return 1
  fi
  for component in "$replicate" "$method_dir" "$tree_type" "$setting"; do
    if [[ "$component" == . || "$component" == .. || -z "$component" ]]; then
      stelarx__outputs_error "unsafe path component in: $relative"
      return 1
    fi
  done

  printf '%s\n%s\n%s\n%s\n' "$replicate" "$method_dir" "$tree_type" "$setting"
}

# Print every results directory beneath an A10K data directory, NUL-separated:
#   <data>/10k-simphy/<replicate>/<method>_outputs/<tree-type>/<setting>
stelarx_list_a10k_results_dirs() {
  local data_dir="${1%/}"
  find "${data_dir}/${STELARX_A10K_REPLICATE_ROOT}" -mindepth 4 -maxdepth 4 -type d \
    -path '*/*_outputs/*/*' -not -name '.*' -print0 2>/dev/null | sort -z
}

# Print the replicate directory names present in an A10K data directory.
stelarx_a10k_replicates() {
  local data_dir="${1%/}"
  find "${data_dir}/${STELARX_A10K_REPLICATE_ROOT}" -mindepth 1 -maxdepth 1 -type d -name 'R*' \
    -printf '%f\n' 2>/dev/null | sort -V
}

# Write the dataset provenance record into a method's mirror directory.
#   stelarx_write_a10k_dataset_record DATA_DIR MIRROR_METHOD_DIR
# The A10K dataset is an external input rather than something this repository
# simulates, so the record states where it lives, which replicates it holds, and
# which files the runs read but the mirror deliberately omits. Content is stable
# for a fixed dataset, so repeated calls rewrite nothing.
stelarx_write_a10k_dataset_record() {
  local data_dir="${1%/}" mirror_method_dir="${2%/}"
  local dataset_name replicates record
  local -a replicate_names=()

  dataset_name="$(basename -- "$data_dir")"
  mapfile -t replicate_names < <(stelarx_a10k_replicates "$data_dir")
  replicates="${replicate_names[*]}"

  record="$(cat <<RECORD
# A10K dataset record (written by scripts/a10k-outputs-dir.sh)
#
# The 10k-astral-dataset is an external input: it is downloaded, not simulated
# by this repository, so there is no simulation command to replay. This record
# pins down what the mirrored results were produced from.
#
# dataset:      ${dataset_name}
# data_dir:     ${data_dir}
# replicates:   ${#replicate_names[@]} (${replicates:-none})
#
# Inputs read by the runs and deliberately left out of the mirror:
#   <data_dir>/${STELARX_A10K_REPLICATE_ROOT}/<R>/s_tree.trees                                # true species tree, used for RF
#   <data_dir>/${STELARX_A10K_REPLICATE_ROOT}/<R>/truegenetrees                               # --tree-type true
#   <data_dir>/${STELARX_A10K_REPLICATE_ROOT}/<R>/estimatedgenetrees/estimatedgenetrees.tre   # --tree-type estimated
#   <data_dir>/${STELARX_A10K_REPLICATE_ROOT}/<R>/estimatedgenetrees/estimatedgenetrees.rooted.tre
#
# The rooted estimated gene trees are derived by run-a10k.sh on first use:
#   ./process_unrooted.sh -i <...>/estimatedgenetrees.tre -o <...>/estimatedgenetrees.rooted.tre -og 0
#
# Per-replicate sizes and modification times of these inputs are recorded in
# <R>/${STELARX_A10K_INPUTS_RECORD}. Each results directory carries the exact
# STELAR-X command that produced it in its out-stelarx.command file.
RECORD
)"

  stelarx__write_file_if_changed "${mirror_method_dir}/${STELARX_A10K_DATASET_RECORD}" "$record"
}

# Record the size and modification time of one replicate's input files.
#   stelarx_write_a10k_replicate_inputs DATA_DIR MIRROR_REPLICATE_DIR REPLICATE
# This is the light fingerprint of the data the mirror does not carry: it lets a
# consumer confirm they hold the same inputs without shipping gigabytes.
stelarx_write_a10k_replicate_inputs() {
  local data_dir="${1%/}" mirror_replicate_dir="${2%/}" replicate="$3"
  local record relative absolute bytes mtime

  record="# inputs of ${replicate} in $(basename -- "$data_dir") (contents are not mirrored)"$'\n'"path	bytes	mtime"
  while IFS= read -r relative; do
    absolute="${data_dir}/${relative}"
    if [[ -f "$absolute" ]]; then
      bytes="$(stat -c '%s' -- "$absolute" 2>/dev/null || echo NA)"
      mtime="$(date -d "@$(stat -c '%Y' -- "$absolute" 2>/dev/null || echo 0)" '+%Y-%m-%dT%H:%M:%S%z' 2>/dev/null || echo NA)"
    else
      bytes="absent"
      mtime="absent"
    fi
    record+=$'\n'"${relative}	${bytes}	${mtime}"
  done < <(stelarx_a10k_replicate_input_paths "$replicate")

  stelarx__write_file_if_changed "${mirror_replicate_dir}/${STELARX_A10K_INPUTS_RECORD}" "$record"
}

# Mirror one A10K results directory into the outputs tree.
#   stelarx_mirror_a10k_results DATA_DIR OUTPUTS_ROOT RESULTS_DIR
# The mirror leaf is replaced atomically so it always equals the source leaf; the
# dataset record and the replicate's input fingerprints are refreshed alongside.
# Prints the mirror leaf path.
stelarx_mirror_a10k_results() {
  local data_dir="${1%/}" outputs_root="${2%/}" results_dir="${3%/}"
  local -a parts=()
  local replicate method_dir tree_type setting
  local mirror_method_dir mirror_replicate_dir mirror_leaf forbidden

  if [[ ! -d "$results_dir" ]]; then
    stelarx__outputs_error "results directory does not exist: $results_dir"
    return 1
  fi
  mapfile -t parts < <(stelarx_a10k_results_components "$data_dir" "$results_dir") || return 1
  [[ ${#parts[@]} -eq 4 ]] || return 1
  replicate="${parts[0]}"; method_dir="${parts[1]}"; tree_type="${parts[2]}"; setting="${parts[3]}"

  forbidden="$(stelarx_a10k_find_forbidden_in_mirror "$results_dir" | head -n1)"
  if [[ -n "$forbidden" ]]; then
    stelarx__outputs_error "refusing to mirror a results directory containing A10K input data: $forbidden"
    return 1
  fi

  mirror_method_dir="${outputs_root}/${method_dir}"
  mirror_replicate_dir="${mirror_method_dir}/${replicate}"
  mirror_leaf="${mirror_replicate_dir}/${tree_type}/${setting}"

  mkdir -p -- "${mirror_leaf%/*}" || return 1
  if ! stelarx_write_a10k_dataset_record "$data_dir" "$mirror_method_dir"; then
    echo "Warning: could not write the A10K dataset record in ${mirror_method_dir}." >&2
  fi
  if ! stelarx_write_a10k_replicate_inputs "$data_dir" "$mirror_replicate_dir" "$replicate"; then
    echo "Warning: could not write the A10K input fingerprints for ${replicate}." >&2
  fi

  stelarx__mirror_leaf_atomic "$results_dir" "$mirror_leaf" || return 1

  printf '%s\n' "$mirror_leaf"
}
