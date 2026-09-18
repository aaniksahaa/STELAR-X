#!/usr/bin/env bash
# Upload the reproducibility mirror of simulated-run outputs to Hugging Face.
#
# Local layout (maintained by test-stelarx-simulated.sh / sync-simulated-outputs.sh):
#   <outputs>/<method>_outputs/<dataset>/<dataset>.command
#   <outputs>/<method>_outputs/<dataset>/<dataset>.params
#   <outputs>/<method>_outputs/<dataset>/<replicate>/<setting>/<trees, CSVs, markers>
#
# Remote layout (same shape, one folder upload per <method>_outputs/<dataset>):
#   <remote-dir>/<method>_outputs/<dataset>/...
#
# Simulated gene trees, true species trees, and SimPhy databases are never part
# of the mirror; every dataset directory is checked for them before upload.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/phylogeny-data-dir.sh"
source "${SCRIPT_DIR}/scripts/hf-python.sh"
source "${SCRIPT_DIR}/scripts/simphy-outputs-dir.sh"

OUTPUTS_DIR=""
DATA_DIR=""
SYNC_FIRST=false
REPO_ID="imAniksahA/blab"
REPO_TYPE="dataset"
REMOTE_DIR="ph/d/simulated/outputs"
UPLOADER="${HOME}/utils/hf-data-transfer/hf_upload.py"
PYTHON_BIN=""
METHODS_RAW="stelarx"
MIN_TAXA=1
MIN_GENE_TREES=1
INCLUDE_INCOMPLETE=true
ALLOW_MISSING_COMMAND=false
DRY_RUN=false
ASSUME_YES=false

print_help() {
  cat <<EOF
upload-bulk-simulated-outputs.sh

Discovers <method>_outputs/<dataset> directories in the simulated outputs
mirror and uploads each one as a folder to the Hugging Face repository under
<remote-dir>/<method>_outputs/<dataset>. The complete plan is shown before one
confirmation. Only inferred trees, CSVs, run markers/logs, and the SimPhy
.command/.params files are uploaded; a dataset directory containing simulated
gene trees, species trees, or SimPhy databases is refused.

Options:
  --outputs-dir PATH       Outputs mirror to upload
                            (default: \$PHYLOGENY_DATA_DIR/outputs/simphy)
  --simphy-outputs-dir PATH
                           Alias for --outputs-dir
  --sync                   Run ./sync-simulated-outputs.sh first so the mirror
                           reflects every result in the data tree
  --data-dir PATH          Data tree used by --sync
                            (default: \$PHYLOGENY_DATA_DIR/simphy/data)
  --method LIST, --methods LIST
                           Only upload these methods, comma/space separated
                            (e.g. "stelarx,aster"; default: ${METHODS_RAW}).
                           Use "all" to upload every method in the mirror.
  --min-taxa N             Minimum taxon count (default: ${MIN_TAXA})
  --min-gene-trees N       Minimum gene-tree count (default: ${MIN_GENE_TREES})
  --exclude-incomplete     Skip "<dataset>_incomplete" datasets
  --allow-missing-command  Upload datasets whose SimPhy .command file is absent
                           (they are otherwise blocked as not reproducible)
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
  ./upload-bulk-simulated-outputs.sh --dry-run
  ./upload-bulk-simulated-outputs.sh --sync
  ./upload-bulk-simulated-outputs.sh --min-taxa 1000
  ./upload-bulk-simulated-outputs.sh --method aster
  ./upload-bulk-simulated-outputs.sh --method all
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

require_positive_integer() {
  local option="$1" value="$2"
  if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: $option requires a positive integer; got '$value'." >&2
    exit 2
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outputs-dir|--simphy-outputs-dir|--data-dir|--simphy-data-dir|--method|--methods|\
    --min-taxa|--min-gene-trees|--repo-id|--repo-type|--remote-dir|--uploader|--python)
      if [[ $# -lt 2 ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 2
      fi
      ;;
  esac

  case "$1" in
    --outputs-dir|--simphy-outputs-dir) OUTPUTS_DIR="$2"; shift 2 ;;
    --outputs-dir=*|--simphy-outputs-dir=*) OUTPUTS_DIR="${1#*=}"; shift ;;
    --data-dir|--simphy-data-dir) DATA_DIR="$2"; shift 2 ;;
    --data-dir=*|--simphy-data-dir=*) DATA_DIR="${1#*=}"; shift ;;
    --sync) SYNC_FIRST=true; shift ;;
    --method|--methods) METHODS_RAW="$2"; shift 2 ;;
    --method=*|--methods=*) METHODS_RAW="${1#*=}"; shift ;;
    --min-taxa) MIN_TAXA="$2"; shift 2 ;;
    --min-taxa=*) MIN_TAXA="${1#*=}"; shift ;;
    --min-gene-trees) MIN_GENE_TREES="$2"; shift 2 ;;
    --min-gene-trees=*) MIN_GENE_TREES="${1#*=}"; shift ;;
    --exclude-incomplete) INCLUDE_INCOMPLETE=false; shift ;;
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

require_positive_integer "--min-taxa" "$MIN_TAXA"
require_positive_integer "--min-gene-trees" "$MIN_GENE_TREES"

OUTPUTS_DIR="$(expand_home "$OUTPUTS_DIR")"
if [[ "$SYNC_FIRST" == true ]]; then
  DATA_DIR="$(stelarx_prepare_simphy_data_dir "$(expand_home "$DATA_DIR")")" || exit 2
  OUTPUTS_DIR="$(stelarx_prepare_simphy_outputs_dir "$OUTPUTS_DIR" "$DATA_DIR")" || exit 2
else
  OUTPUTS_DIR="$(stelarx_prepare_simphy_outputs_dir_standalone "$OUTPUTS_DIR")" || exit 2
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

# The outputs mirror is uploaded folder by folder, which needs the folder-aware
# hf_upload.py (its --help lists --include/--exclude). Older copies of the
# helper only accept single files and fail on every dataset directory.
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

if [[ "$SYNC_FIRST" == true ]]; then
  echo "Refreshing the outputs mirror from the data tree first..."
  SYNC_CMD=("${SCRIPT_DIR}/sync-simulated-outputs.sh" --simphy-data-dir "$DATA_DIR" --simphy-outputs-dir "$OUTPUTS_DIR" --quiet)
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

# Count replicate directories and distinct setting names under a dataset mirror.
count_replicates() {
  find "$1" -mindepth 1 -maxdepth 1 -type d -not -name '.*' | wc -l | tr -d ' '
}
count_settings() {
  find "$1" -mindepth 2 -maxdepth 2 -type d -not -name '.*' -printf '%f\n' 2>/dev/null | sort -u | wc -l | tr -d ' '
}
count_result_files() {
  find "$1" -mindepth 3 -type f | wc -l | tr -d ' '
}

echo "STELAR-X simulated outputs uploader"
echo "Outputs directory: $OUTPUTS_DIR"
echo "Repository:        $REPO_ID ($REPO_TYPE)"
echo "Remote path:       $REMOTE_DIR/"
echo "Python:            $PYTHON_BIN"
echo "Selection:         taxa >= $MIN_TAXA, gene trees >= $MIN_GENE_TREES$(
  printf ', methods: %s' "$([[ ${#METHOD_FILTER[@]} -gt 0 ]] && echo "${!METHOD_FILTER[*]}" || echo 'all')"
)$([[ "$INCLUDE_INCOMPLETE" == false ]] && printf ', excluding _incomplete')"
[[ "$DRY_RUN" == true ]] && echo "Dry run:           yes"
echo

declare -a SELECTED=()      # "<method_dir>/<dataset>" entries in upload order
declare -a IGNORED=()       # informational, not uploaded
blocked=0

while IFS= read -r -d '' method_path; do
  method_dir="${method_path##*/}"
  if [[ ${#METHOD_FILTER[@]} -gt 0 && -z "${METHOD_FILTER[$method_dir]:-}" ]]; then
    continue
  fi
  while IFS= read -r -d '' dataset_path; do
    dataset_name="${dataset_path##*/}"
    if ! stelarx_simphy_dataset_name_is_valid "$dataset_name"; then
      IGNORED+=("${method_dir}/${dataset_name}: unrecognized dataset name")
      continue
    fi
    taxa="${BASH_REMATCH[1]}"
    gene_trees="${BASH_REMATCH[2]}"
    incomplete="${BASH_REMATCH[8]:-}"
    if (( taxa < MIN_TAXA || gene_trees < MIN_GENE_TREES )); then
      continue
    fi
    if [[ -n "$incomplete" && "$INCLUDE_INCOMPLETE" == false ]]; then
      continue
    fi
    SELECTED+=("${method_dir}/${dataset_name}")
  done < <(find "$method_path" -mindepth 1 -maxdepth 1 -type d -name 't_*' -print0 | sort -zV)
done < <(find "$OUTPUTS_DIR" -mindepth 1 -maxdepth 1 -type d -name '*_outputs' -print0 | sort -z)

if [[ ${#IGNORED[@]} -gt 0 ]]; then
  echo "Ignored directories:"
  for entry in "${IGNORED[@]}"; do
    echo "  $entry"
  done
  echo
fi

if [[ ${#SELECTED[@]} -eq 0 ]]; then
  echo "No <method>_outputs/<dataset> directories matched the selection under $OUTPUTS_DIR."
  echo "Run ./sync-simulated-outputs.sh (or this tool with --sync) to populate the mirror."
  exit 0
fi

printf '%-16s %-7s %-7s %-5s %-5s %-7s %-8s %-28s %s\n' \
  "METHOD" "TAXA" "GENES" "REPL" "SETS" "FILES" "SIZE" "ACTION" "DATASET"
printf '%-16s %-7s %-7s %-5s %-5s %-7s %-8s %-28s %s\n' \
  "----------------" "-------" "-------" "-----" "-----" "-------" "--------" "----------------------------" "-------"

declare -a UPLOADS=()
for entry in "${SELECTED[@]}"; do
  method_dir="${entry%%/*}"
  dataset_name="${entry#*/}"
  dataset_path="${OUTPUTS_DIR}/${entry}"
  stelarx_simphy_dataset_name_is_valid "$dataset_name"
  taxa="${BASH_REMATCH[1]}"
  gene_trees="${BASH_REMATCH[2]}"

  replicates="$(count_replicates "$dataset_path")"
  settings="$(count_settings "$dataset_path")"
  files="$(count_result_files "$dataset_path")"
  size="$(human_size "$dataset_path")"
  action="upload"

  forbidden="$(stelarx_simphy_find_forbidden_in_mirror "$dataset_path" | head -n1)"
  command_file="${dataset_path}/${dataset_name}.command"
  base_command_file="${dataset_path}/${dataset_name%_incomplete}.command"

  if [[ -n "$forbidden" ]]; then
    action="BLOCKED: contains ${forbidden##*/}"
    ((blocked++)) || true
  elif (( files == 0 )); then
    action="BLOCKED: no result files"
    ((blocked++)) || true
  elif [[ ! -f "$command_file" && ! -f "$base_command_file" ]]; then
    if [[ "$ALLOW_MISSING_COMMAND" == true ]]; then
      action="upload (NO .command)"
    else
      action="BLOCKED: missing .command"
      ((blocked++)) || true
    fi
  elif [[ "$dataset_name" == *_incomplete && ! -f "$command_file" ]]; then
    # Base SimPhy command is present; the derivation command is optional but
    # recorded by newer sim_incomplete.sh runs.
    action="upload (base .command only)"
  fi

  printf '%-16s %-7s %-7s %-5s %-5s %-7s %-8s %-28s %s\n' \
    "${method_dir%_outputs}" "$taxa" "$gene_trees" "$replicates" "$settings" "$files" "$size" "$action" "$dataset_name"

  if [[ "$action" != BLOCKED:* ]]; then
    UPLOADS+=("$entry")
  fi
done

echo
if (( blocked > 0 )); then
  echo "Error: $blocked selected dataset director(ies) cannot be uploaded safely." >&2
  echo "Fix the blocked entries above (./sync-simulated-outputs.sh restores missing .command files), then rerun." >&2
  echo "Nothing was uploaded." >&2
  exit 1
fi

# Collapse sorted replicate names into "R1-R4, R6" style ranges; names that
# are not R<n> are listed verbatim.
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

echo "Cases and settings to upload (<dataset> / <replicates> / <setting>):"
total_leaves=0
for entry in "${UPLOADS[@]}"; do
  dataset_path="${OUTPUTS_DIR}/${entry}"
  dataset_name="${entry#*/}"
  declare -A setting_replicates=()
  while IFS= read -r -d '' leaf; do
    rel="${leaf#"$dataset_path"/}"
    setting="${rel#*/}"
    setting_replicates["$setting"]+="${rel%%/*} "
    ((total_leaves++)) || true
  done < <(find "$dataset_path" -mindepth 2 -maxdepth 2 -type d -not -name '.*' -print0 | sort -zV)
  while IFS= read -r setting; do
    [[ -z "$setting" ]] && continue
    read -r -a repl_names <<< "${setting_replicates[$setting]}"
    echo "  [${entry%%/*}] ${dataset_name} / $(format_replicate_ranges "${repl_names[@]}") / ${setting}"
  done < <(printf '%s\n' "${!setting_replicates[@]}" | sort)
  unset setting_replicates
done
echo
echo "Plan: upload ${#UPLOADS[@]} dataset director(ies) as folders (${total_leaves} replicate/setting result set(s))."
echo "Destinations:"
for entry in "${UPLOADS[@]}"; do
  echo "  ${entry} -> ${REPO_ID}/${REMOTE_DIR}/${entry}/"
done

build_upload_command() {
  local entry="$1"
  UPLOAD_CMD=("$PYTHON_BIN" "$UPLOADER"
    --repo-id "$REPO_ID"
    --repo-type "$REPO_TYPE"
    --local-path "${OUTPUTS_DIR}/${entry}"
    --path-in-repo "${REMOTE_DIR}/${entry}"
    --commit-message "Simulated outputs: ${entry}")
}

if [[ "$DRY_RUN" == true ]]; then
  echo
  echo "Upload commands:"
  for entry in "${UPLOADS[@]}"; do
    build_upload_command "$entry"
    print_command "${UPLOAD_CMD[@]}"
  done
  echo "Dry run complete; nothing was uploaded."
  exit 0
fi

if [[ "$ASSUME_YES" == false ]]; then
  echo
  read -r -p "Proceed with ${#UPLOADS[@]} folder upload(s)? [y/N]: " confirm
  if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Cancelled; nothing was uploaded."
    exit 0
  fi
fi

echo
echo "Revalidating every dataset directory immediately before upload..."
postcheck_failed=0
for entry in "${UPLOADS[@]}"; do
  dataset_path="${OUTPUTS_DIR}/${entry}"
  forbidden="$(stelarx_simphy_find_forbidden_in_mirror "$dataset_path" | head -n1)"
  if [[ -n "$forbidden" ]]; then
    echo "  Error: simulated input data appeared in the mirror: $forbidden" >&2
    ((postcheck_failed++)) || true
  elif ! directory_is_nonempty "$dataset_path"; then
    echo "  Error: dataset directory disappeared or is empty: $dataset_path" >&2
    ((postcheck_failed++)) || true
  else
    echo "  Ready: ${entry} ($(human_size "$dataset_path"))"
  fi
done
if (( postcheck_failed > 0 )); then
  echo "Error: pre-upload validation failed for $postcheck_failed director(ies); nothing was uploaded." >&2
  exit 1
fi

succeeded=0
failed=0
for i in "${!UPLOADS[@]}"; do
  entry="${UPLOADS[$i]}"
  echo
  echo "[$((i + 1))/${#UPLOADS[@]}] Uploading: $entry"
  build_upload_command "$entry"
  if "${UPLOAD_CMD[@]}"; then
    ((succeeded++)) || true
    echo "  Done: ${REMOTE_DIR}/${entry}/"
  else
    ((failed++)) || true
    echo "  Error: upload failed; continuing with remaining directories." >&2
  fi
done

echo
echo "Summary: selected=${#UPLOADS[@]} uploaded=$succeeded failed=$failed"
if (( failed > 0 )); then
  exit 1
fi
