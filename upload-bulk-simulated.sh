#!/usr/bin/env bash
# Package discovered SimPhy datasets and upload them to Hugging Face.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/phylogeny-data-dir.sh"
source "${SCRIPT_DIR}/scripts/hf-python.sh"

DATA_DIR=""
REPO_ID="imAniksahA/blab"
REPO_TYPE="dataset"
REMOTE_DIR="ph/d/simulated/astralx-datasets/raw"
UPLOADER="${HOME}/utils/hf-data-transfer/hf_upload.py"
PYTHON_BIN=""
MIN_TAXA=1000
MIN_GENE_TREES=1000
DRY_RUN=false
ASSUME_YES=false

print_help() {
  cat <<EOF
upload-bulk-simulated.sh

Discovers canonical SimPhy dataset directories and ZIP archives directly under
simphy/data. Missing, stale, or invalid ZIPs are created/rebuilt before their
datasets are uploaded. The complete plan is shown before one confirmation.

Options:
  --data-dir PATH          Directory containing datasets and ZIPs
                            (default: \$PHYLOGENY_DATA_DIR/simphy/data)
  --min-taxa N             Minimum taxon count (default: ${MIN_TAXA})
  --min-gene-trees N       Minimum gene-tree count (default: ${MIN_GENE_TREES})
  --all                    Select all positive taxa/gene-tree counts
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

Directories and ZIPs ending in "_incomplete" are intentionally excluded.
Existing current ZIPs are reused. Missing, stale, and invalid ZIPs with a
matching non-empty source directory are safely built before upload.

Examples:
  ./upload-bulk-simulated.sh --dry-run
  ./upload-bulk-simulated.sh
  ./upload-bulk-simulated.sh --all --remote-dir ph/d/simulated/stelarx-datasets/raw
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
    --data-dir|--min-taxa|--min-gene-trees|--repo-id|--repo-type|\
    --remote-dir|--uploader|--python)
      if [[ $# -lt 2 ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 2
      fi
      ;;
  esac

  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --data-dir=*) DATA_DIR="${1#*=}"; shift ;;
    --min-taxa) MIN_TAXA="$2"; shift 2 ;;
    --min-taxa=*) MIN_TAXA="${1#*=}"; shift ;;
    --min-gene-trees) MIN_GENE_TREES="$2"; shift 2 ;;
    --min-gene-trees=*) MIN_GENE_TREES="${1#*=}"; shift ;;
    --all) MIN_TAXA=1; MIN_GENE_TREES=1; shift ;;
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

DATA_DIR="$(expand_home "$DATA_DIR")"
DATA_DIR="$(stelarx_prepare_simphy_data_dir "$DATA_DIR")"
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
if ! command -v unzip >/dev/null 2>&1; then
  echo "Error: unzip is required but was not found." >&2
  exit 2
fi
if ! command -v zip >/dev/null 2>&1; then
  echo "Error: zip is required but was not found." >&2
  exit 2
fi
if [[ ! -f "$UPLOADER" ]]; then
  echo "Error: uploader was not found: $UPLOADER" >&2
  exit 2
fi
PYTHON_BIN="$(stelarx_find_hf_python "$PYTHON_BIN")" || exit 2

dataset_name_is_valid() {
  local name="$1"
  [[ "$name" =~ ^t_([1-9][0-9]*)_g_([1-9][0-9]*)_sb_([0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?)_spmin_([1-9][0-9]*)_spmax_([1-9][0-9]*)$ ]]
}

archive_has_expected_root() {
  local archive="$1" dataset_name="$2"
  local entry normalized found=false

  if ! unzip -tq "$archive" >/dev/null 2>&1; then
    return 1
  fi

  while IFS= read -r entry; do
    normalized="${entry#./}"
    [[ -z "$normalized" ]] && continue
    if [[ "$normalized" == /* || "$normalized" =~ (^|/)[.][.](/|$) ]]; then
      return 1
    fi
    case "$normalized" in
      "$dataset_name"|"$dataset_name/"*) found=true ;;
      *) return 1 ;;
    esac
  done < <(unzip -Z1 "$archive")

  [[ "$found" == true ]]
}

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

CURRENT_TEMP_DIR=""
cleanup() {
  if [[ -n "$CURRENT_TEMP_DIR" && -d "$CURRENT_TEMP_DIR" ]]; then
    rm -rf -- "$CURRENT_TEMP_DIR"
  fi
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

archive_needs_rebuild() {
  local dataset_name="$1"
  local dataset_path="${DATA_DIR}/${dataset_name}"
  local archive_path="${DATA_DIR}/${dataset_name}.zip"

  [[ ! -f "$archive_path" ]] && return 0
  archive_has_expected_root "$archive_path" "$dataset_name" || return 0
  if directory_is_nonempty "$dataset_path" &&
     [[ -n "$(find "$dataset_path" -type f -newer "$archive_path" -print -quit 2>/dev/null)" ]]; then
    return 0
  fi
  return 1
}

build_archive() {
  local dataset_name="$1"
  local dataset_path="${DATA_DIR}/${dataset_name}"
  local archive_path="${DATA_DIR}/${dataset_name}.zip"
  local temp_archive

  if ! directory_is_nonempty "$dataset_path"; then
    echo "  Error: source directory is missing or empty: $dataset_path" >&2
    return 1
  fi

  CURRENT_TEMP_DIR="$(mktemp -d "${DATA_DIR}/.${dataset_name}.archive.XXXXXX")" || {
    echo "  Error: could not create a temporary archive directory." >&2
    return 1
  }
  temp_archive="${CURRENT_TEMP_DIR}/${dataset_name}.zip"

  echo "  Compressing..."
  if ! (cd "$DATA_DIR" && zip -rq "$temp_archive" "$dataset_name"); then
    echo "  Error: zip failed; any existing archive was left untouched." >&2
    rm -rf -- "$CURRENT_TEMP_DIR"
    CURRENT_TEMP_DIR=""
    return 1
  fi

  echo "  Testing archive integrity and layout..."
  if ! archive_has_expected_root "$temp_archive" "$dataset_name"; then
    echo "  Error: archive validation failed; any existing archive was left untouched." >&2
    rm -rf -- "$CURRENT_TEMP_DIR"
    CURRENT_TEMP_DIR=""
    return 1
  fi

  if ! mv -f -- "$temp_archive" "$archive_path"; then
    echo "  Error: could not install the validated archive." >&2
    rm -rf -- "$CURRENT_TEMP_DIR"
    CURRENT_TEMP_DIR=""
    return 1
  fi
  rmdir "$CURRENT_TEMP_DIR"
  CURRENT_TEMP_DIR=""
  echo "  Ready: $archive_path ($(human_size "$archive_path"))"
}

declare -A CANDIDATES=()
declare -a ARCHIVES=()
needs_archive=0
blocked=0

# Dataset directories are the primary source. Valid orphan ZIPs are also
# supported so exported archives can be uploaded after their sources move.
while IFS= read -r -d '' dataset_path; do
  dataset_name="${dataset_path##*/}"
  if ! dataset_name_is_valid "$dataset_name"; then
    continue
  fi
  taxa="${BASH_REMATCH[1]}"
  gene_trees="${BASH_REMATCH[2]}"
  if (( taxa >= MIN_TAXA && gene_trees >= MIN_GENE_TREES )); then
    CANDIDATES["$dataset_name"]=1
  fi
done < <(find "$DATA_DIR" -maxdepth 1 -mindepth 1 -type d -name 't_*' -print0 | sort -zV)

while IFS= read -r -d '' archive_path; do
  archive_name="${archive_path##*/}"
  dataset_name="${archive_name%.zip}"
  if ! dataset_name_is_valid "$dataset_name"; then
    continue
  fi
  taxa="${BASH_REMATCH[1]}"
  gene_trees="${BASH_REMATCH[2]}"
  if (( taxa >= MIN_TAXA && gene_trees >= MIN_GENE_TREES )); then
    CANDIDATES["$dataset_name"]=1
  fi
done < <(find "$DATA_DIR" -maxdepth 1 -mindepth 1 -type f -name 't_*_g_*.zip' -print0 | sort -zV)

echo "STELAR-X simulated dataset uploader"
echo "Data directory: $DATA_DIR"
echo "Repository:     $REPO_ID ($REPO_TYPE)"
echo "Remote path:    $REMOTE_DIR/"
echo "Python:         $PYTHON_BIN"
echo "Selection:      taxa >= $MIN_TAXA, gene trees >= $MIN_GENE_TREES"
[[ "$DRY_RUN" == true ]] && echo "Dry run:        yes"
echo

if [[ ${#CANDIDATES[@]} -eq 0 ]]; then
  echo "No complete dataset directories or canonical ZIPs matched the selection."
  exit 0
fi

printf '%-7s %-7s %-11s %-10s %-27s %s\n' \
  "TAXA" "GENES" "SOURCE" "ZIP" "ACTION" "DATASET"
printf '%-7s %-7s %-11s %-10s %-27s %s\n' \
  "-------" "-------" "-----------" "----------" "---------------------------" "-------"

while IFS= read -r dataset_name; do
  [[ -z "$dataset_name" ]] && continue
  dataset_name_is_valid "$dataset_name"

  taxa="${BASH_REMATCH[1]}"
  gene_trees="${BASH_REMATCH[2]}"
  dataset_path="${DATA_DIR}/${dataset_name}"
  archive_path="${DATA_DIR}/${dataset_name}.zip"
  source_size="-"
  zip_size="missing"
  source_available=false

  if directory_is_nonempty "$dataset_path"; then
    source_available=true
    source_size="$(human_size "$dataset_path")"
  elif [[ -d "$dataset_path" ]]; then
    source_size="empty"
  fi

  if [[ -f "$archive_path" ]]; then
    zip_size="$(human_size "$archive_path")"
  fi

  action="reuse ZIP, then upload"
  if [[ ! -f "$archive_path" ]]; then
    if [[ "$source_available" == true ]]; then
      action="CREATE ZIP, then upload"
      ((needs_archive++)) || true
    else
      action="BLOCKED: no usable source"
      ((blocked++)) || true
    fi
  elif ! archive_has_expected_root "$archive_path" "$dataset_name"; then
    if [[ "$source_available" == true ]]; then
      action="REBUILD invalid ZIP, upload"
      ((needs_archive++)) || true
    else
      action="BLOCKED: invalid ZIP"
      ((blocked++)) || true
    fi
  elif [[ "$source_available" == true ]] &&
       [[ -n "$(find "$dataset_path" -type f -newer "$archive_path" -print -quit 2>/dev/null)" ]]; then
    action="REBUILD stale ZIP, upload"
    ((needs_archive++)) || true
  fi

  printf '%-7s %-7s %-11s %-10s %-27s %s\n' \
    "$taxa" "$gene_trees" "$source_size" "$zip_size" "$action" "$dataset_name"

  if [[ "$action" != BLOCKED:* ]]; then
    ARCHIVES+=("$archive_path")
  fi
done < <(printf '%s\n' "${!CANDIDATES[@]}" | sort -V)

echo
if (( blocked > 0 )); then
  echo "Error: $blocked selected dataset(s) cannot be safely prepared." >&2
  echo "Fix the blocked source/archive entries shown above, then rerun." >&2
  exit 1
fi

echo "Plan: create/rebuild $needs_archive ZIP(s), reuse $((${#ARCHIVES[@]} - needs_archive)) ZIP(s), upload ${#ARCHIVES[@]} archive(s)."
echo "Destinations:"
for i in "${!ARCHIVES[@]}"; do
  archive_name="${ARCHIVES[$i]##*/}"
  echo "  ${archive_name} -> ${REPO_ID}/${REMOTE_DIR}/${archive_name}"
done

if [[ "$DRY_RUN" == true ]]; then
  echo
  if (( needs_archive > 0 )); then
    echo "ZIP preparation is performed internally and atomically before upload."
  fi
  echo "Upload commands:"
  for archive_path in "${ARCHIVES[@]}"; do
    archive_name="${archive_path##*/}"
    print_command "$PYTHON_BIN" "$UPLOADER" \
      --repo-id "$REPO_ID" \
      --repo-type "$REPO_TYPE" \
      --local-path "$archive_path" \
      --path-in-repo "${REMOTE_DIR}/${archive_name}"
  done
  echo "Dry run complete; no ZIPs were created and nothing was uploaded."
  exit 0
fi

if [[ "$ASSUME_YES" == false ]]; then
  echo
  read -r -p "Proceed with $needs_archive ZIP operation(s) and ${#ARCHIVES[@]} upload(s)? [y/N]: " confirm
  if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Cancelled; no ZIPs were changed and nothing was uploaded."
    exit 0
  fi
fi

if (( needs_archive > 0 )); then
  echo
  echo "Preparing missing, stale, or invalid ZIP archives..."
  prepared=0
  prepare_failed=0
  for archive_path in "${ARCHIVES[@]}"; do
    archive_name="${archive_path##*/}"
    dataset_name="${archive_name%.zip}"
    if ! archive_needs_rebuild "$dataset_name"; then
      continue
    fi

    ((prepared++)) || true
    echo
    echo "[$prepared/$needs_archive] Creating or rebuilding: $archive_name"
    if ! build_archive "$dataset_name"; then
      ((prepare_failed++)) || true
    fi
  done

  if (( prepare_failed > 0 )); then
    echo "Error: $prepare_failed archive operation(s) failed; uploads were not started." >&2
    exit 1
  fi
fi

echo
echo "Revalidating every archive immediately before upload..."
postcheck_failed=0
for i in "${!ARCHIVES[@]}"; do
  archive_path="${ARCHIVES[$i]}"
  archive_name="${archive_path##*/}"
  dataset_name="${archive_name%.zip}"
  dataset_path="${DATA_DIR}/${dataset_name}"

  if [[ ! -f "$archive_path" ]] || ! archive_has_expected_root "$archive_path" "$dataset_name"; then
    echo "  Error: missing or invalid archive: $archive_path" >&2
    ((postcheck_failed++)) || true
  elif directory_is_nonempty "$dataset_path" &&
       [[ -n "$(find "$dataset_path" -type f -newer "$archive_path" -print -quit 2>/dev/null)" ]]; then
    echo "  Error: source changed after archive creation: $dataset_path" >&2
    ((postcheck_failed++)) || true
  else
    echo "  Ready: $archive_name ($(human_size "$archive_path"))"
  fi
done

if (( postcheck_failed > 0 )); then
  echo "Error: pre-upload validation failed for $postcheck_failed archive(s); nothing was uploaded." >&2
  exit 1
fi

succeeded=0
failed=0
for i in "${!ARCHIVES[@]}"; do
  archive_path="${ARCHIVES[$i]}"
  archive_name="${archive_path##*/}"
  echo
  echo "[$((i + 1))/${#ARCHIVES[@]}] Uploading: $archive_name"
  if "$PYTHON_BIN" "$UPLOADER" \
      --repo-id "$REPO_ID" \
      --repo-type "$REPO_TYPE" \
      --local-path "$archive_path" \
      --path-in-repo "${REMOTE_DIR}/${archive_name}"; then
    ((succeeded++)) || true
    echo "  Done: ${REMOTE_DIR}/${archive_name}"
  else
    ((failed++)) || true
    echo "  Error: upload failed; continuing with remaining archives." >&2
  fi
done

echo
echo "Summary: selected=${#ARCHIVES[@]} uploaded=$succeeded failed=$failed"
if (( failed > 0 )); then
  exit 1
fi
