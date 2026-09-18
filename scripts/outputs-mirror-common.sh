#!/usr/bin/env bash

# Shared primitives for the reproducibility output mirrors.
#
# Two mirrors are built on top of this file and behave identically:
#   * scripts/simphy-outputs-dir.sh  - SimPhy simulated runs
#   * scripts/a10k-outputs-dir.sh    - the 10k-astral-dataset (A10K) runs
#
# Both keep the data tree untouched and maintain a separate, small, shareable
# copy of every results directory (inferred trees, CSVs, run markers/logs,
# command records) plus the provenance of the dataset itself. Simulated or
# downloaded inputs - gene trees, species trees, databases - are never copied.

[[ -n "${STELARX_OUTPUTS_MIRROR_COMMON_SOURCED:-}" ]] && return 0
STELARX_OUTPUTS_MIRROR_COMMON_SOURCED=1

stelarx__outputs_error() {
  echo "Error: $*" >&2
}

stelarx__path_is_within() {
  # stelarx__path_is_within CHILD PARENT -> true when CHILD == PARENT or is below it.
  local child="${1%/}" parent="${2%/}"
  [[ "$child" == "$parent" || "$child" == "$parent"/* ]]
}

# Copy one small file atomically (temp file + rename), preserving timestamps.
stelarx__copy_small_file_atomic() {
  local source="$1" destination="$2" temp
  temp="$(mktemp -- "${destination}.tmp.XXXXXX")" || return 1
  if cp -p -- "$source" "$temp" && mv -f -- "$temp" "$destination"; then
    return 0
  fi
  rm -f -- "$temp"
  return 1
}

# Write CONTENT to PATH only when it differs, so re-running a mirror never
# rewrites an unchanged provenance record. Returns 0 on success.
stelarx__write_file_if_changed() {
  local destination="$1" content="$2" temp
  if [[ -f "$destination" ]] && [[ "$(cat -- "$destination" 2>/dev/null)" == "$content" ]]; then
    return 0
  fi
  mkdir -p -- "$(dirname -- "$destination")" || return 1
  temp="$(mktemp -- "${destination}.tmp.XXXXXX")" || return 1
  if printf '%s\n' "$content" > "$temp" && mv -f -- "$temp" "$destination"; then
    return 0
  fi
  rm -f -- "$temp"
  return 1
}

# True when a file name matches one of the given glob patterns.
stelarx__name_matches_any() {
  local name="$1" pattern
  shift
  for pattern in "$@"; do
    # shellcheck disable=SC2053
    [[ "$name" == $pattern ]] && return 0
  done
  return 1
}

# Print every file beneath ROOT whose name matches one of the given patterns.
#   stelarx__find_names_under ROOT PATTERN [PATTERN...]
stelarx__find_names_under() {
  local root="$1"
  shift
  local -a find_args=()
  local pattern first=true
  for pattern in "$@"; do
    if [[ "$first" == true ]]; then
      first=false
    else
      find_args+=(-o)
    fi
    find_args+=(-name "$pattern")
  done
  find "$root" -type f \( "${find_args[@]}" \) -print 2>/dev/null
}

# Resolve, validate, and create a mirror root that must stay disjoint from the
# data directory.
#   stelarx__resolve_outputs_root REQUESTED_DIR RESOLVED_DATA_DIR LABEL
stelarx__resolve_outputs_root() {
  local requested_dir="${1:-}" data_dir="${2:-}" label="${3:-Outputs}"
  local resolved_dir

  if [[ -z "$requested_dir" ]]; then
    stelarx__outputs_error "internal: an outputs directory is required."
    return 2
  fi
  if [[ -z "$data_dir" ]]; then
    stelarx__outputs_error "internal: resolved ${label} data directory is required."
    return 2
  fi
  if [[ "$requested_dir" == "~/"* ]]; then
    requested_dir="${HOME}/${requested_dir:2}"
  fi
  if [[ -e "$requested_dir" && ! -d "$requested_dir" ]]; then
    stelarx__outputs_error "${label} outputs path exists but is not a directory: $requested_dir"
    return 2
  fi

  resolved_dir="$(realpath -m -- "$requested_dir")" || {
    stelarx__outputs_error "could not resolve ${label} outputs directory: $requested_dir"
    return 2
  }
  if stelarx__path_is_within "$resolved_dir" "$data_dir"; then
    stelarx__outputs_error "${label} outputs directory must not be the data directory or inside it: $resolved_dir"
    return 2
  fi
  if stelarx__path_is_within "$data_dir" "$resolved_dir"; then
    stelarx__outputs_error "${label} outputs directory must not contain the data directory: $resolved_dir"
    return 2
  fi
  if ! mkdir -p -- "$resolved_dir"; then
    stelarx__outputs_error "could not create ${label} outputs directory: $resolved_dir"
    return 2
  fi
  if ! resolved_dir="$(cd "$resolved_dir" && pwd -P)"; then
    stelarx__outputs_error "could not enter ${label} outputs directory: $requested_dir"
    return 2
  fi
  printf '%s\n' "$resolved_dir"
}

# Resolve a mirror root without a data directory, for tools that only read the
# mirror. An explicit path wins; otherwise PHYLOGENY_DATA_DIR is required.
#   stelarx__resolve_outputs_root_standalone REQUESTED_DIR DEFAULT_SUFFIX LABEL
stelarx__resolve_outputs_root_standalone() {
  local requested_dir="${1:-}" default_suffix="${2:-}" label="${3:-Outputs}"
  local resolved_dir

  if [[ -z "$requested_dir" ]]; then
    if [[ -z "${PHYLOGENY_DATA_DIR:-}" ]]; then
      stelarx__outputs_error "PHYLOGENY_DATA_DIR is not set. Set it or pass an explicit ${label} outputs directory."
      return 2
    fi
    requested_dir="${PHYLOGENY_DATA_DIR%/}/${default_suffix#/}"
  fi
  if [[ "$requested_dir" == "~/"* ]]; then
    requested_dir="${HOME}/${requested_dir:2}"
  fi
  if [[ -e "$requested_dir" && ! -d "$requested_dir" ]]; then
    stelarx__outputs_error "${label} outputs path exists but is not a directory: $requested_dir"
    return 2
  fi
  if ! mkdir -p -- "$requested_dir"; then
    stelarx__outputs_error "could not create ${label} outputs directory: $requested_dir"
    return 2
  fi
  if ! resolved_dir="$(cd "$requested_dir" && pwd -P)"; then
    stelarx__outputs_error "could not resolve ${label} outputs directory: $requested_dir"
    return 2
  fi
  printf '%s\n' "$resolved_dir"
}

# Replace MIRROR_LEAF with an exact copy of RESULTS_DIR, atomically: the copy is
# built in a temporary sibling and swapped in, so stale files from earlier runs
# never survive and a reader never sees a half-written leaf.
#   stelarx__mirror_leaf_atomic RESULTS_DIR MIRROR_LEAF
stelarx__mirror_leaf_atomic() {
  local results_dir="${1%/}" mirror_leaf="${2%/}"
  local mirror_parent setting temp_new temp_old

  mirror_parent="$(dirname -- "$mirror_leaf")"
  setting="$(basename -- "$mirror_leaf")"
  mkdir -p -- "$mirror_parent" || return 1

  temp_new="$(mktemp -d -- "${mirror_parent}/.${setting}.mirror.XXXXXX")" || return 1
  if ! cp -a -- "${results_dir}/." "${temp_new}/"; then
    rm -rf -- "$temp_new"
    stelarx__outputs_error "could not copy results into the outputs mirror: $results_dir"
    return 1
  fi

  temp_old=""
  if [[ -e "$mirror_leaf" ]]; then
    temp_old="$(mktemp -d -u -- "${mirror_parent}/.${setting}.previous.XXXXXX")"
    if ! mv -- "$mirror_leaf" "$temp_old"; then
      rm -rf -- "$temp_new"
      stelarx__outputs_error "could not replace the existing mirror: $mirror_leaf"
      return 1
    fi
  fi
  if ! mv -- "$temp_new" "$mirror_leaf"; then
    [[ -n "$temp_old" && -d "$temp_old" ]] && mv -- "$temp_old" "$mirror_leaf" 2>/dev/null
    rm -rf -- "$temp_new"
    stelarx__outputs_error "could not install the outputs mirror: $mirror_leaf"
    return 1
  fi
  [[ -n "$temp_old" ]] && rm -rf -- "$temp_old"
  return 0
}
