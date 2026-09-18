#!/usr/bin/env bash
# Upload the reproducibility mirror of A10K run outputs to Hugging Face.
#
# Local layout (maintained by run-a10k.sh / sync-a10k-outputs.sh):
#   <outputs>/<method>_outputs/a10k-dataset.command
#   <outputs>/<method>_outputs/<replicate>/inputs.tsv
#   <outputs>/<method>_outputs/<replicate>/<tree-type>/<setting>/<trees, CSVs, markers>
#
# Remote layout (same shape, one folder upload per <method>_outputs/<replicate>):
#   <remote-dir>/<method>_outputs/<replicate>/...
#
# A10K gene trees and species trees are never part of the mirror; every
# replicate directory is checked for them before upload.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/hf-python.sh"
source "${SCRIPT_DIR}/scripts/a10k-outputs-dir.sh"

OUTPUTS_DIR=""
DATA_DIR=""
SYNC_FIRST=false
REPO_ID="imAniksahA/blab"
REPO_TYPE="dataset"
REMOTE_DIR="ph/d/a10k/outputs"
UPLOADER="${HOME}/utils/hf-data-transfer/hf_upload.py"
PYTHON_BIN=""
METHODS_RAW="stelarx"
REPLICATES_SPEC=""
ALLOW_MISSING_COMMAND=false
DRY_RUN=false
ASSUME_YES=false

print_help() {
  cat <<EOF
upload-a10k-outputs.sh

Discovers <method>_outputs/<replicate> directories in the A10K outputs mirror
and uploads each one as a folder to the Hugging Face repository under
<remote-dir>/<method>_outputs/<replicate>, together with each method's
${STELARX_A10K_DATASET_RECORD} dataset record. The complete plan is shown before
one confirmation. Only inferred trees, CSVs, run markers/logs, command records
and input fingerprints are uploaded; a replicate directory containing A10K gene
trees or species trees is refused.

Options:
  --outputs-dir PATH       Outputs mirror to upload
                            (default: \$PHYLOGENY_DATA_DIR/outputs/10k-astral-dataset)
  --a10k-outputs-dir PATH  Alias for --outputs-dir
  --sync                   Run ./sync-a10k-outputs.sh first so the mirror
                           reflects every result in the data tree
  --data-dir PATH          A10K dataset root used by --sync
                            (default: \$PHYLOGENY_DATA_DIR/10k-astral-dataset)
  --method LIST, --methods LIST
                           Only upload these methods, comma/space separated
                            (e.g. "stelarx,aster"; default: ${METHODS_RAW}).
                           Use "all" to upload every method in the mirror.
  --replicates SPEC        Only upload these replicates, e.g. "1-10" or "R1,R3"
                            (default: all replicates in the mirror)
  --allow-missing-command  Upload methods whose ${STELARX_A10K_DATASET_RECORD}
                           is absent (they are otherwise blocked as not
                           reproducible; ./sync-a10k-outputs.sh restores it)
  --repo-id ID             Hugging Face repository (default: ${REPO_ID})
  --repo-type TYPE         dataset, model, or space (default: ${REPO_TYPE})
  --remote-dir PATH        Destination directory inside the repository
                            (default: ${REMOTE_DIR})
  --uploader PATH          Path to hf_upload.py (default: ${UPLOADER})
  --python COMMAND         Python interpreter with huggingface_hub (default: the
                           first of python3, python, conda base python that has it)
  --dry-run                Validate and print commands without uploading
  --yes, -y                Do not ask for confirmation
  --help, -h               Show this message

Re-running is cheap: files already present on the Hub are skipped by the
uploader, so this doubles as an incremental sync of the outputs mirror.

Examples:
  ./upload-a10k-outputs.sh --dry-run
  ./upload-a10k-outputs.sh --sync
  ./upload-a10k-outputs.sh --replicates 1-10
  ./upload-a10k-outputs.sh --method all
EOF
}

expand_home() {
  local path="$1"
  if [[ "$path" == "~/"* ]]; then
    printf '%s/%s\n' "$HOME" "${path:2}"
  else
    printf '%s\n' "$path"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outputs-dir|--a10k-outputs-dir|--data-dir|--a10k-data-dir|--method|--methods|\
    --replicates|--repo-id|--repo-type|--remote-dir|--uploader|--python)
      if [[ $# -lt 2 ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 2
      fi
      ;;
  esac

  case "$1" in
    --outputs-dir|--a10k-outputs-dir) OUTPUTS_DIR="$2"; shift 2 ;;
    --outputs-dir=*|--a10k-outputs-dir=*) OUTPUTS_DIR="${1#*=}"; shift ;;
    --data-dir|--a10k-data-dir) DATA_DIR="$2"; shift 2 ;;
    --data-dir=*|--a10k-data-dir=*) DATA_DIR="${1#*=}"; shift ;;
    --sync) SYNC_FIRST=true; shift ;;
    --method|--methods) METHODS_RAW="$2"; shift 2 ;;
    --method=*|--methods=*) METHODS_RAW="${1#*=}"; shift ;;
    --replicates) REPLICATES_SPEC="$2"; shift 2 ;;
    --replicates=*) REPLICATES_SPEC="${1#*=}"; shift ;;
    --allow-missing-command) ALLOW_MISSING_COMMAND=true; shift ;;
    --repo-id) REPO_ID="$2"; shift 2 ;;
    --repo-id=*) REPO_ID="${1#*=}"; shift ;;
    --repo-type) REPO_TYPE="$2"; shift 2 ;;
    --repo-type=*) REPO_TYPE="${1#*=}"; shift ;;
    --remote-dir) REMOTE_DIR="$2"; shift 2 ;;
    --remote-dir=*) REMOTE_DIR="${1#*=}"; shift ;;
    --uploader) UPLOADER="$2"; shift 2 ;;
    --uploader=*) UPLOADER="${1#*=}"; shift ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --python=*) PYTHON_BIN="${1#*=}"; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    --yes|-y) ASSUME_YES=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *)
      echo "Error: unknown option '$1'." >&2
      print_help >&2
      exit 2
      ;;
  esac
done

OUTPUTS_DIR="$(expand_home "$OUTPUTS_DIR")"
if [[ "$SYNC_FIRST" == true ]]; then
  if [[ -z "$DATA_DIR" ]]; then
    if [[ -z "${PHYLOGENY_DATA_DIR:-}" ]]; then
      echo "Error: --sync needs --data-dir (or PHYLOGENY_DATA_DIR)." >&2
      exit 2
    fi
    DATA_DIR="${PHYLOGENY_DATA_DIR%/}/10k-astral-dataset"
  fi
  DATA_DIR="$(expand_home "$DATA_DIR")"
  if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: A10K data directory does not exist: $DATA_DIR" >&2
    exit 2
  fi
  DATA_DIR="$(cd "$DATA_DIR" && pwd -P)"
  OUTPUTS_DIR="$(stelarx_prepare_a10k_outputs_dir "$OUTPUTS_DIR" "$DATA_DIR")" || exit 2
else
  OUTPUTS_DIR="$(stelarx_prepare_a10k_outputs_dir_standalone "$OUTPUTS_DIR")" || exit 2
fi
UPLOADER="$(realpath -m "$(expand_home "$UPLOADER")")"
REMOTE_DIR="${REMOTE_DIR%/}"

if [[ ! "$REPO_ID" =~ ^[^/[:space:]]+/[^/[:space:]]+$ ]]; then
  echo "Error: invalid --repo-id '$REPO_ID'; expected owner/repository." >&2
  exit 2
fi
case "$REPO_TYPE" in
  dataset|model|space) ;;
  *) echo "Error: --repo-type must be dataset, model, or space." >&2; exit 2 ;;
esac
if [[ -z "$REMOTE_DIR" || "$REMOTE_DIR" == /* || "$REMOTE_DIR" =~ (^|/)[.][.](/|$) ]]; then
  echo "Error: unsafe --remote-dir '$REMOTE_DIR'." >&2
  exit 2
fi
if [[ ! -f "$UPLOADER" ]]; then
  echo "Error: uploader was not found: $UPLOADER" >&2
  exit 2
fi
PYTHON_BIN="$(stelarx_find_hf_python "$PYTHON_BIN")" || exit 2

# The mirror is uploaded folder by folder, which needs the folder-aware
# hf_upload.py (its --help lists --include/--exclude). Older copies of the
# helper only accept single files and fail on every replicate directory.
if [[ "$UPLOADER" == *.py ]] &&
   ! "$PYTHON_BIN" "$UPLOADER" --help 2>/dev/null | grep -q -- '--include'; then
  echo "Error: $UPLOADER does not support folder uploads (no --include/--exclude in its --help)." >&2
  echo "  This is an outdated hf_upload.py. Replace it with the current folder-aware version" >&2
  echo "  (the one whose --help mentions '--local-path may be a folder'), then rerun." >&2
  exit 2
fi

# "all" (or an empty value) uploads every method found in the mirror.
declare -A METHOD_FILTER=()
if [[ "${METHODS_RAW,,}" == "all" ]]; then
  METHODS_RAW=""
fi
if [[ -n "$METHODS_RAW" ]]; then
  read -r -a method_items <<< "${METHODS_RAW//,/ }"
  for method in "${method_items[@]}"; do
    [[ -z "$method" ]] && continue
    [[ "$method" == *_outputs ]] || method="${method}_outputs"
    METHOD_FILTER["$method"]=1
  done
  if [[ ${#METHOD_FILTER[@]} -eq 0 ]]; then
    echo "Error: --method did not name any method." >&2
    exit 2
  fi
fi

# --replicates accepts "1-10", "R1,R3" or "1,3"; empty means every replicate.
declare -A REPLICATE_FILTER=()
if [[ -n "$REPLICATES_SPEC" ]]; then
  if [[ "$REPLICATES_SPEC" =~ ^[0-9]+-[0-9]+$ ]]; then
    for i in $(seq "${REPLICATES_SPEC%-*}" "${REPLICATES_SPEC#*-}"); do
      REPLICATE_FILTER["R${i}"]=1
    done
  else
    read -r -a replicate_items <<< "${REPLICATES_SPEC//,/ }"
    for replicate in "${replicate_items[@]}"; do
      [[ -z "$replicate" ]] && continue
      [[ "$replicate" == R* ]] || replicate="R${replicate}"
      REPLICATE_FILTER["$replicate"]=1
    done
  fi
  if [[ ${#REPLICATE_FILTER[@]} -eq 0 ]]; then
    echo "Error: --replicates did not name any replicate." >&2
    exit 2
  fi
fi

if [[ "$SYNC_FIRST" == true ]]; then
  echo "Refreshing the outputs mirror from the data tree first..."
  SYNC_CMD=("${SCRIPT_DIR}/sync-a10k-outputs.sh" --data-dir "$DATA_DIR" --a10k-outputs-dir "$OUTPUTS_DIR" --quiet)
  [[ -n "$METHODS_RAW" ]] && SYNC_CMD+=(--methods "$METHODS_RAW")
  [[ "$DRY_RUN" == true ]] && SYNC_CMD+=(--dry-run)
  if ! "${SYNC_CMD[@]}"; then
    echo "Error: mirror sync reported problems; nothing was uploaded." >&2
    exit 1
  fi
  echo
fi

human_size() {
  du -sh -- "$1" 2>/dev/null | awk '{print $1}'
}

print_command() {
  printf '  '
  printf '%q ' "$@"
  printf '\n'
}

directory_is_nonempty() {
  [[ -d "$1" ]] && [[ -n "$(find "$1" -mindepth 1 -print -quit 2>/dev/null)" ]]
}

# Distinct tree types, settings and result files under one replicate mirror.
count_tree_types() {
  find "$1" -mindepth 1 -maxdepth 1 -type d -not -name '.*' | wc -l | tr -d ' '
}
count_settings() {
  find "$1" -mindepth 2 -maxdepth 2 -type d -not -name '.*' -printf '%f\n' 2>/dev/null | sort -u | wc -l | tr -d ' '
}
count_result_files() {
  find "$1" -mindepth 3 -type f | wc -l | tr -d ' '
}

# Files kept directly under a method's mirror directory - the dataset record and
# any merged score summary - uploaded individually, NUL-separated.
method_top_level_files() {
  find "${OUTPUTS_DIR}/${1}" -mindepth 1 -maxdepth 1 -type f -not -name '.*' -print0 2>/dev/null | sort -z
}

echo "STELAR-X A10K outputs uploader"
echo "Outputs directory: $OUTPUTS_DIR"
echo "Repository:        $REPO_ID ($REPO_TYPE)"
echo "Remote path:       $REMOTE_DIR/"
echo "Python:            $PYTHON_BIN"
echo "Selection:         methods: $([[ ${#METHOD_FILTER[@]} -gt 0 ]] && echo "${!METHOD_FILTER[*]}" || echo 'all'), replicates: $([[ ${#REPLICATE_FILTER[@]} -gt 0 ]] && echo "${REPLICATES_SPEC}" || echo 'all')"
[[ "$DRY_RUN" == true ]] && echo "Dry run:           yes"
echo

declare -a SELECTED=()   # "<method_dir>/<replicate>" entries in upload order
declare -a RECORDS=()    # "<method_dir>" entries whose dataset record is uploaded
declare -a IGNORED=()    # informational, not uploaded
blocked=0

while IFS= read -r -d '' method_path; do
  method_dir="${method_path##*/}"
  if [[ ${#METHOD_FILTER[@]} -gt 0 && -z "${METHOD_FILTER[$method_dir]:-}" ]]; then
    continue
  fi
  method_has_replicates=false
  while IFS= read -r -d '' replicate_path; do
    replicate="${replicate_path##*/}"
    if [[ ! "$replicate" =~ ^R[0-9]+$ ]]; then
      IGNORED+=("${method_dir}/${replicate}: unrecognized replicate name")
      continue
    fi
    if [[ ${#REPLICATE_FILTER[@]} -gt 0 && -z "${REPLICATE_FILTER[$replicate]:-}" ]]; then
      continue
    fi
    SELECTED+=("${method_dir}/${replicate}")
    method_has_replicates=true
  done < <(find "$method_path" -mindepth 1 -maxdepth 1 -type d -not -name '.*' -print0 | sort -zV)
  [[ "$method_has_replicates" == true ]] && RECORDS+=("$method_dir")
done < <(find "$OUTPUTS_DIR" -mindepth 1 -maxdepth 1 -type d -name '*_outputs' -print0 | sort -z)

if [[ ${#IGNORED[@]} -gt 0 ]]; then
  echo "Ignored directories:"
  for entry in "${IGNORED[@]}"; do
    echo "  $entry"
  done
  echo
fi

if [[ ${#SELECTED[@]} -eq 0 ]]; then
  echo "No <method>_outputs/<replicate> directories matched the selection under $OUTPUTS_DIR."
  echo "Run ./sync-a10k-outputs.sh (or this tool with --sync) to populate the mirror."
  exit 0
fi

printf '%-16s %-10s %-6s %-6s %-7s %-8s %-28s\n' \
  "METHOD" "REPLICATE" "TYPES" "SETS" "FILES" "SIZE" "ACTION"
printf '%-16s %-10s %-6s %-6s %-7s %-8s %-28s\n' \
  "----------------" "----------" "------" "------" "-------" "--------" "----------------------------"

declare -a UPLOADS=()
for entry in "${SELECTED[@]}"; do
  method_dir="${entry%%/*}"
  replicate="${entry#*/}"
  replicate_path="${OUTPUTS_DIR}/${entry}"
  record_file="${OUTPUTS_DIR}/${method_dir}/${STELARX_A10K_DATASET_RECORD}"

  tree_types="$(count_tree_types "$replicate_path")"
  settings="$(count_settings "$replicate_path")"
  files="$(count_result_files "$replicate_path")"
  size="$(human_size "$replicate_path")"
  action="upload"

  forbidden="$(stelarx_a10k_find_forbidden_in_mirror "$replicate_path" | head -n1)"
  if [[ -n "$forbidden" ]]; then
    action="BLOCKED: contains ${forbidden##*/}"
    ((blocked++)) || true
  elif (( files == 0 )); then
    action="BLOCKED: no result files"
    ((blocked++)) || true
  elif [[ ! -f "$record_file" ]]; then
    if [[ "$ALLOW_MISSING_COMMAND" == true ]]; then
      action="upload (NO dataset record)"
    else
      action="BLOCKED: missing ${STELARX_A10K_DATASET_RECORD}"
      ((blocked++)) || true
    fi
  fi

  printf '%-16s %-10s %-6s %-6s %-7s %-8s %-28s\n' \
    "${method_dir%_outputs}" "$replicate" "$tree_types" "$settings" "$files" "$size" "$action"

  if [[ "$action" != BLOCKED:* ]]; then
    UPLOADS+=("$entry")
  fi
done

# The dataset record and merged summaries sit directly under <method>_outputs and
# are uploaded as individual files, so they need the same forbidden-name guard
# the replicate directories get.
for method_dir in "${RECORDS[@]}"; do
  while IFS= read -r -d '' top_file; do
    top_name="${top_file##*/}"
    if stelarx_a10k_name_is_forbidden_in_mirror "$top_name"; then
      printf '%-16s %-10s %-6s %-6s %-7s %-8s %-28s\n' \
        "${method_dir%_outputs}" "-" "-" "-" "-" "-" "BLOCKED: contains ${top_name}"
      ((blocked++)) || true
    fi
  done < <(method_top_level_files "$method_dir")
done

echo
if (( blocked > 0 )); then
  echo "Error: $blocked selected entr(ies) cannot be uploaded safely." >&2
  echo "Fix the blocked entries above (./sync-a10k-outputs.sh restores the dataset record), then rerun." >&2
  echo "Nothing was uploaded." >&2
  exit 1
fi

# Collapse sorted replicate names into "R1-R4, R6" style ranges.
format_replicate_ranges() {
  local out="" start="" prev="" name n
  local -a others=()
  for name in "$@"; do
    if [[ ! "$name" =~ ^R([0-9]+)$ ]]; then
      others+=("$name"); continue
    fi
    n="${BASH_REMATCH[1]}"
    if [[ -z "$start" ]]; then start=$n; prev=$n; continue; fi
    if (( n == prev + 1 )); then prev=$n; continue; fi
    out+="${out:+, }R${start}"; (( start != prev )) && out+="-R${prev}"
    start=$n; prev=$n
  done
  if [[ -n "$start" ]]; then
    out+="${out:+, }R${start}"; (( start != prev )) && out+="-R${prev}"
  fi
  for name in "${others[@]}"; do out+="${out:+, }${name}"; done
  printf '%s' "$out"
}

echo "Cases and settings to upload (<replicates> / <tree type> / <setting>):"
declare -A case_replicates=()
total_leaves=0
for entry in "${UPLOADS[@]}"; do
  method_dir="${entry%%/*}"
  replicate="${entry#*/}"
  replicate_path="${OUTPUTS_DIR}/${entry}"
  while IFS= read -r -d '' leaf; do
    rel="${leaf#"$replicate_path"/}"
    case_replicates["${method_dir}|${rel}"]+="${replicate} "
    ((total_leaves++)) || true
  done < <(find "$replicate_path" -mindepth 2 -maxdepth 2 -type d -not -name '.*' -print0 | sort -zV)
done
while IFS= read -r key; do
  [[ -z "$key" ]] && continue
  read -r -a repl_names <<< "${case_replicates[$key]}"
  method_dir="${key%%|*}"
  rel="${key#*|}"
  echo "  [${method_dir}] $(format_replicate_ranges "${repl_names[@]}") / ${rel%%/*} / ${rel#*/}"
done < <(printf '%s\n' "${!case_replicates[@]}" | sort)

echo
echo "Plan: upload ${#UPLOADS[@]} replicate director(ies) as folders (${total_leaves} tree-type/setting result set(s))"
record_count=0
for method_dir in "${RECORDS[@]}"; do
  while IFS= read -r -d '' _record_file; do
    ((record_count++)) || true
  done < <(method_top_level_files "$method_dir")
done
echo "      plus ${record_count} dataset record / summary file(s)."
echo "Destinations:"
for method_dir in "${RECORDS[@]}"; do
  while IFS= read -r -d '' record_file; do
    record_name="${record_file##*/}"
    echo "  ${method_dir}/${record_name} -> ${REPO_ID}/${REMOTE_DIR}/${method_dir}/${record_name}"
  done < <(method_top_level_files "$method_dir")
done
for entry in "${UPLOADS[@]}"; do
  echo "  ${entry} -> ${REPO_ID}/${REMOTE_DIR}/${entry}/"
done

build_upload_command() {
  local local_path="$1" path_in_repo="$2" message="$3"
  UPLOAD_CMD=("$PYTHON_BIN" "$UPLOADER"
    --repo-id "$REPO_ID"
    --repo-type "$REPO_TYPE"
    --local-path "$local_path"
    --path-in-repo "$path_in_repo"
    --commit-message "$message")
}

# One entry per upload: "<local path>|<path in repo>|<label>".
declare -a UPLOAD_JOBS=()
for method_dir in "${RECORDS[@]}"; do
  while IFS= read -r -d '' record_file; do
    record_name="${record_file##*/}"
    UPLOAD_JOBS+=("${record_file}|${REMOTE_DIR}/${method_dir}/${record_name}|${method_dir}/${record_name}")
  done < <(method_top_level_files "$method_dir")
done
for entry in "${UPLOADS[@]}"; do
  UPLOAD_JOBS+=("${OUTPUTS_DIR}/${entry}|${REMOTE_DIR}/${entry}|${entry}")
done

if [[ "$DRY_RUN" == true ]]; then
  echo
  echo "Upload commands:"
  for job in "${UPLOAD_JOBS[@]}"; do
    build_upload_command "${job%%|*}" "$(cut -d'|' -f2 <<< "$job")" "A10K outputs: $(cut -d'|' -f3 <<< "$job")"
    print_command "${UPLOAD_CMD[@]}"
  done
  echo "Dry run complete; nothing was uploaded."
  exit 0
fi

if [[ "$ASSUME_YES" == false ]]; then
  echo
  read -r -p "Proceed with ${#UPLOAD_JOBS[@]} upload(s)? [y/N]: " confirm
  if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Cancelled; nothing was uploaded."
    exit 0
  fi
fi

echo
echo "Revalidating every replicate directory immediately before upload..."
postcheck_failed=0
for entry in "${UPLOADS[@]}"; do
  replicate_path="${OUTPUTS_DIR}/${entry}"
  forbidden="$(stelarx_a10k_find_forbidden_in_mirror "$replicate_path" | head -n1)"
  if [[ -n "$forbidden" ]]; then
    echo "  Error: A10K input data appeared in the mirror: $forbidden" >&2
    ((postcheck_failed++)) || true
  elif ! directory_is_nonempty "$replicate_path"; then
    echo "  Error: replicate directory disappeared or is empty: $replicate_path" >&2
    ((postcheck_failed++)) || true
  else
    echo "  Ready: ${entry} ($(human_size "$replicate_path"))"
  fi
done
for method_dir in "${RECORDS[@]}"; do
  while IFS= read -r -d '' top_file; do
    top_name="${top_file##*/}"
    if stelarx_a10k_name_is_forbidden_in_mirror "$top_name"; then
      echo "  Error: A10K input data appeared in the mirror: ${method_dir}/${top_name}" >&2
      ((postcheck_failed++)) || true
    fi
  done < <(method_top_level_files "$method_dir")
done
if (( postcheck_failed > 0 )); then
  echo "Error: pre-upload validation failed for $postcheck_failed director(ies); nothing was uploaded." >&2
  exit 1
fi

succeeded=0
failed=0
for i in "${!UPLOAD_JOBS[@]}"; do
  job="${UPLOAD_JOBS[$i]}"
  label="$(cut -d'|' -f3 <<< "$job")"
  echo
  echo "[$((i + 1))/${#UPLOAD_JOBS[@]}] Uploading: $label"
  build_upload_command "${job%%|*}" "$(cut -d'|' -f2 <<< "$job")" "A10K outputs: ${label}"
  if "${UPLOAD_CMD[@]}"; then
    ((succeeded++)) || true
    echo "  Done: $(cut -d'|' -f2 <<< "$job")"
  else
    ((failed++)) || true
    echo "  Error: upload failed; continuing with remaining directories." >&2
  fi
done

echo
echo "Summary: selected=${#UPLOAD_JOBS[@]} uploaded=$succeeded failed=$failed"
if (( failed > 0 )); then
  exit 1
fi
