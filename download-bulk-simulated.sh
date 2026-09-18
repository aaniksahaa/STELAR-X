#!/usr/bin/env bash
# Download simulated STELAR-X datasets from Hugging Face.
#
# Edit the parameter lists below or override them with the matching command-line
# options. Every Cartesian-product combination maps to one dataset archive.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/phylogeny-data-dir.sh"

# --------------------------------
# Parameter lists (edit as needed)
# --------------------------------
T_LIST=(175000 200000)
G_LIST=(1000)
SB_LIST=(0.000001)
SPMIN_LIST=(100000)
SPMAX_LIST=(200000)

REPO_ID="imAniksahA/blab"
REPO_TYPE="dataset"
REMOTE_DIR="ph/d/simulated/astralx-datasets/raw"
LOCAL_DATA_DIR=""
DOWNLOAD_SCRIPT="${HOME}/utils/hf-data-transfer/hf_download.sh"

DOWNLOAD_ONLY=false
KEEP_ZIP=false
DRY_RUN=false

print_help() {
  cat <<EOF
download-bulk-simulated.sh

Downloads every dataset formed by the configured parameter-list Cartesian
product. By default, each ZIP is validated, extracted, and then deleted.

Options:
  --taxa-list LIST          Taxon counts, comma- or space-separated
  --gene-trees-list LIST    Gene-tree counts, comma- or space-separated
  --sb-list LIST            Substitution/birthrate values
  --spmin-list LIST         Minimum population sizes
  --spmax-list LIST         Maximum population sizes
  --repo-id ID              Hugging Face repository (default: ${REPO_ID})
  --repo-type TYPE          dataset, model, or space (default: ${REPO_TYPE})
  --local-dir PATH          Local simphy data directory
                            (default: \$PHYLOGENY_DATA_DIR/simphy/data)
  --download-script PATH    Path to hf_download.sh
                            (default: ${DOWNLOAD_SCRIPT})
  --download-only           Download and validate ZIPs without extracting them
  --no-unzip                Alias for --download-only
  --keep-zip                Keep each ZIP after successful extraction
  --dry-run                 Print generated paths and commands without changes
  --help, -h                Show this message

Existing non-empty dataset directories are skipped. Existing ZIPs are reused
after validation. Failed downloads and extractions are retained/reported, and
processing continues with the remaining combinations.

Examples:
  ./download-bulk-simulated.sh
  ./download-bulk-simulated.sh --download-only --taxa-list "175000 200000"
  ./download-bulk-simulated.sh --keep-zip --local-dir ./simphy/data
EOF
}

parse_list() {
  local raw="${1//,/ }"
  local -n destination="$2"
  read -r -a destination <<< "$raw"
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
    --taxa-list|--gene-trees-list|--sb-list|--spmin-list|--spmax-list|\
    --repo-id|--repo-type|--local-dir|--download-script)
      if [[ $# -lt 2 ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 2
      fi
      ;;
  esac
  case "$1" in
    --taxa-list)       parse_list "$2" T_LIST; shift 2 ;;
    --taxa-list=*)     parse_list "${1#*=}" T_LIST; shift ;;
    --gene-trees-list)   parse_list "$2" G_LIST; shift 2 ;;
    --gene-trees-list=*) parse_list "${1#*=}" G_LIST; shift ;;
    --sb-list)         parse_list "$2" SB_LIST; shift 2 ;;
    --sb-list=*)       parse_list "${1#*=}" SB_LIST; shift ;;
    --spmin-list)      parse_list "$2" SPMIN_LIST; shift 2 ;;
    --spmin-list=*)    parse_list "${1#*=}" SPMIN_LIST; shift ;;
    --spmax-list)      parse_list "$2" SPMAX_LIST; shift 2 ;;
    --spmax-list=*)    parse_list "${1#*=}" SPMAX_LIST; shift ;;
    --repo-id)         REPO_ID="$2"; shift 2 ;;
    --repo-type)       REPO_TYPE="$2"; shift 2 ;;
    --local-dir)       LOCAL_DATA_DIR="$2"; shift 2 ;;
    --download-script) DOWNLOAD_SCRIPT="$2"; shift 2 ;;
    --download-only|--no-unzip) DOWNLOAD_ONLY=true; shift ;;
    --keep-zip)        KEEP_ZIP=true; shift ;;
    --dry-run)         DRY_RUN=true; shift ;;
    --help|-h)         print_help; exit 0 ;;
    *)
      echo "Error: unknown option '$1'." >&2
      print_help >&2
      exit 2
      ;;
  esac
done

LOCAL_DATA_DIR="$(expand_home "$LOCAL_DATA_DIR")"
DOWNLOAD_SCRIPT="$(expand_home "$DOWNLOAD_SCRIPT")"
LOCAL_DATA_DIR="$(stelarx_prepare_simphy_data_dir "$LOCAL_DATA_DIR")"
REMOTE_DIR="${REMOTE_DIR%/}"

validate_integer_list() {
  local label="$1"
  local -n values="$2"
  local value

  if [[ ${#values[@]} -eq 0 ]]; then
    echo "Error: $label cannot be empty." >&2
    return 1
  fi
  for value in "${values[@]}"; do
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
      echo "Error: invalid $label value '$value'; expected a positive integer." >&2
      return 1
    fi
  done
}

validate_decimal_list() {
  local label="$1"
  local -n values="$2"
  local value

  if [[ ${#values[@]} -eq 0 ]]; then
    echo "Error: $label cannot be empty." >&2
    return 1
  fi
  for value in "${values[@]}"; do
    if [[ ! "$value" =~ ^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$ ]]; then
      echo "Error: invalid $label value '$value'." >&2
      return 1
    fi
  done
}

validate_remote_dir() {
  local remote_dir="$1"
  if [[ -z "$remote_dir" || "$remote_dir" == /* || "$remote_dir" =~ (^|/)\.\.(/|$) ]]; then
    echo "Error: unsafe remote directory '$remote_dir'." >&2
    return 1
  fi
}

if [[ ! "$REPO_ID" =~ ^[^/[:space:]]+/[^/[:space:]]+$ ]]; then
  echo "Error: invalid --repo-id '$REPO_ID'; expected owner/repository." >&2
  exit 2
fi
case "$REPO_TYPE" in
  dataset|model|space) ;;
  *) echo "Error: --repo-type must be dataset, model, or space." >&2; exit 2 ;;
esac

validate_integer_list "taxa-list" T_LIST || exit 2
validate_integer_list "gene-trees-list" G_LIST || exit 2
validate_decimal_list "sb-list" SB_LIST || exit 2
validate_integer_list "spmin-list" SPMIN_LIST || exit 2
validate_integer_list "spmax-list" SPMAX_LIST || exit 2
validate_remote_dir "$REMOTE_DIR" || exit 2

if [[ "$DOWNLOAD_SCRIPT" == */* ]]; then
  if [[ ! -x "$DOWNLOAD_SCRIPT" ]]; then
    echo "Error: download script is not executable: $DOWNLOAD_SCRIPT" >&2
    exit 2
  fi
elif ! command -v "$DOWNLOAD_SCRIPT" >/dev/null 2>&1; then
  echo "Error: download command was not found: $DOWNLOAD_SCRIPT" >&2
  exit 2
fi

if [[ "$DRY_RUN" == false ]] && ! command -v unzip >/dev/null 2>&1; then
  echo "Error: unzip is required to validate downloaded archives." >&2
  exit 2
fi

WORK_ROOT=""
CURRENT_EXTRACT_DIR=""
cleanup() {
  if [[ -n "$CURRENT_EXTRACT_DIR" && -d "$CURRENT_EXTRACT_DIR" ]]; then
    rm -rf -- "$CURRENT_EXTRACT_DIR"
  fi
  if [[ -n "$WORK_ROOT" && -d "$WORK_ROOT" ]]; then
    rm -rf -- "$WORK_ROOT"
  fi
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

if [[ "$DRY_RUN" == false ]]; then
  WORK_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-download.XXXXXX")" || exit 2
fi

print_command() {
  printf '  '
  printf '%q ' "$@"
  printf '\n'
}

directory_is_nonempty() {
  [[ -d "$1" ]] && [[ -n "$(find "$1" -mindepth 1 -print -quit 2>/dev/null)" ]]
}

validate_archive() {
  local zip_path="$1"
  local dataset_name="$2"
  local entry normalized
  local found_dataset=false

  if ! unzip -tq "$zip_path" >/dev/null; then
    echo "  Error: ZIP validation failed: $zip_path" >&2
    return 1
  fi

  while IFS= read -r entry; do
    normalized="${entry#./}"
    if [[ -z "$normalized" ]]; then
      continue
    fi
    if [[ "$normalized" == /* || "$normalized" =~ (^|/)\.\.(/|$) ]]; then
      echo "  Error: archive contains an unsafe path: $entry" >&2
      return 1
    fi
    case "$normalized" in
      "$dataset_name"|"$dataset_name/"*) found_dataset=true ;;
      *)
        echo "  Error: archive entry is outside '$dataset_name/': $entry" >&2
        return 1
        ;;
    esac
  done < <(unzip -Z1 "$zip_path")

  if [[ "$found_dataset" == false ]]; then
    echo "  Error: archive does not contain '$dataset_name/'." >&2
    return 1
  fi
}

TOTAL=0
DOWNLOADED=0
EXTRACTED=0
VALIDATED=0
SKIPPED=0
PLANNED=0
FAILED=0

process_dataset() {
  local t="$1" g="$2" sb="$3" spmin="$4" spmax="$5"
  local dataset_name="t_${t}_g_${g}_sb_${sb}_spmin_${spmin}_spmax_${spmax}"
  local archive_name="${dataset_name}.zip"
  local remote_path zip_path dataset_dir download_dir extract_dir extracted_dataset
  local -a download_cmd

  ((TOTAL++)) || true

  if (( spmin > spmax )); then
    echo "Error: spmin=$spmin exceeds spmax=$spmax for $dataset_name." >&2
    ((FAILED++)) || true
    return
  fi

  remote_path="${REMOTE_DIR}/${archive_name}"
  zip_path="${LOCAL_DATA_DIR}/${archive_name}"
  dataset_dir="${LOCAL_DATA_DIR}/${dataset_name}"
  download_cmd=("$DOWNLOAD_SCRIPT" --repo-id "$REPO_ID" --repo-type "$REPO_TYPE"
                --path-in-repo "$remote_path" --local-path "$zip_path")

  echo
  echo "[$TOTAL] $dataset_name"
  echo "  Remote: $REPO_ID/$remote_path"
  echo "  Local:  $dataset_dir"

  if [[ "$DOWNLOAD_ONLY" == false ]] && directory_is_nonempty "$dataset_dir"; then
    echo "  Skip: dataset directory already exists and is non-empty."
    ((SKIPPED++)) || true
    return
  fi

  if [[ "$DRY_RUN" == true ]]; then
    if [[ ! -f "$zip_path" ]]; then
      echo "  Download command:"
      print_command "${download_cmd[@]}"
    else
      echo "  Reuse existing ZIP: $zip_path"
    fi
    if [[ "$DOWNLOAD_ONLY" == false ]]; then
      echo "  Would validate and extract to: $dataset_dir"
      if [[ "$KEEP_ZIP" == false ]]; then
        echo "  Would delete after successful extraction: $zip_path"
      fi
    fi
    ((PLANNED++)) || true
    return
  fi

  if [[ ! -f "$zip_path" ]]; then
    download_dir="${WORK_ROOT}/${TOTAL}"
    mkdir -p "$download_dir"
    echo "  Downloading ZIP..."
    if ! (cd "$download_dir" && "${download_cmd[@]}"); then
      echo "  Error: download failed for $remote_path" >&2
      ((FAILED++)) || true
      return
    fi
    if [[ ! -f "$zip_path" ]]; then
      echo "  Error: downloader succeeded but did not create $zip_path" >&2
      ((FAILED++)) || true
      return
    fi
    ((DOWNLOADED++)) || true
  else
    echo "  Reusing existing ZIP."
  fi

  echo "  Validating ZIP..."
  if ! validate_archive "$zip_path" "$dataset_name"; then
    echo "  ZIP retained for inspection: $zip_path" >&2
    ((FAILED++)) || true
    return
  fi
  ((VALIDATED++)) || true

  if [[ "$DOWNLOAD_ONLY" == true ]]; then
    echo "  Ready: $zip_path"
    return
  fi

  if [[ -e "$dataset_dir" ]]; then
    if [[ ! -d "$dataset_dir" ]]; then
      echo "  Error: extraction target exists but is not a directory: $dataset_dir" >&2
      ((FAILED++)) || true
      return
    fi
    if directory_is_nonempty "$dataset_dir"; then
      echo "  Skip: dataset appeared while processing; ZIP retained."
      ((SKIPPED++)) || true
      return
    fi
  fi

  extract_dir="$(mktemp -d "${LOCAL_DATA_DIR}/.${dataset_name}.extract.XXXXXX")" || {
    echo "  Error: could not create temporary extraction directory." >&2
    ((FAILED++)) || true
    return
  }
  CURRENT_EXTRACT_DIR="$extract_dir"

  echo "  Extracting ZIP..."
  if ! unzip -q "$zip_path" -d "$extract_dir"; then
    echo "  Error: extraction failed; ZIP retained: $zip_path" >&2
    rm -rf -- "$extract_dir"
    CURRENT_EXTRACT_DIR=""
    ((FAILED++)) || true
    return
  fi

  extracted_dataset="${extract_dir}/${dataset_name}"
  if ! directory_is_nonempty "$extracted_dataset"; then
    echo "  Error: extracted dataset is missing or empty; ZIP retained." >&2
    rm -rf -- "$extract_dir"
    CURRENT_EXTRACT_DIR=""
    ((FAILED++)) || true
    return
  fi

  if [[ -d "$dataset_dir" ]] && ! rmdir "$dataset_dir"; then
    echo "  Error: extraction target is no longer empty; ZIP retained." >&2
    rm -rf -- "$extract_dir"
    CURRENT_EXTRACT_DIR=""
    ((FAILED++)) || true
    return
  fi
  if ! mv "$extracted_dataset" "$dataset_dir"; then
    echo "  Error: could not install extracted dataset; ZIP retained." >&2
    rm -rf -- "$extract_dir"
    CURRENT_EXTRACT_DIR=""
    ((FAILED++)) || true
    return
  fi
  rmdir "$extract_dir"
  CURRENT_EXTRACT_DIR=""
  ((EXTRACTED++)) || true
  echo "  Extracted: $dataset_dir"

  if [[ "$KEEP_ZIP" == false ]]; then
    if rm -f -- "$zip_path"; then
      echo "  Removed ZIP after successful extraction."
    else
      echo "  Warning: could not remove ZIP: $zip_path" >&2
    fi
  else
    echo "  Kept ZIP: $zip_path"
  fi
}

echo "STELAR-X simulated dataset downloader"
echo "Repository: $REPO_ID ($REPO_TYPE)"
echo "Local data: $LOCAL_DATA_DIR"
if [[ "$DOWNLOAD_ONLY" == true ]]; then
  echo "Mode:       download and validate only"
elif [[ "$KEEP_ZIP" == true ]]; then
  echo "Mode:       download, extract, and keep ZIP"
else
  echo "Mode:       download, extract, and remove ZIP"
fi
[[ "$DRY_RUN" == true ]] && echo "Dry run:    yes"

for t in "${T_LIST[@]}"; do
  for g in "${G_LIST[@]}"; do
    for sb in "${SB_LIST[@]}"; do
      for spmin in "${SPMIN_LIST[@]}"; do
        for spmax in "${SPMAX_LIST[@]}"; do
          process_dataset "$t" "$g" "$sb" "$spmin" "$spmax"
        done
      done
    done
  done
done

echo
echo "Summary: total=$TOTAL downloaded=$DOWNLOADED validated=$VALIDATED extracted=$EXTRACTED skipped=$SKIPPED planned=$PLANNED failed=$FAILED"

if (( FAILED > 0 )); then
  exit 1
fi
