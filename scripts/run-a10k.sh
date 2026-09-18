#!/usr/bin/env bash
# run-a10k.sh
# STELAR-X runner for the 10k-simphy dataset layout.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/a10k-outputs-dir.sh"

# Propagate terminal color preference to Java subprocesses even when this
# script's output is piped through tee into the per-run log below.
[[ -t 1 || -t 2 ]] && export FORCE_COLOR=1

NTFY_CHANNEL_NAME="${NTFY_CHANNEL_NAME:-anik-phylo-stx}"

# Exact invocation of this script, appended to each run's command record.
SCRIPT_ARGV=("$0" "$@")

TREE_TYPES_RAW="estimated"
DATA_DIR=""
A10K_OUTPUTS_DIR=""
OUTPUTS_MIRROR=true
REPLICATES_SPEC=""
START_REP=""
END_REP=""
FRESH=false
STELARX_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STELARX_OPTS="--search-space S2 -vv"
STELARX_OPTS_LIST_RAW=""
TIME_MONITOR=true
GPU_MONITOR=true
NO_NOTIFY=false

source "${STELARX_ROOT}/experiment-setting-name.sh"

csv_get_field() {
  local file="$1"
  shift
  local header data
  header="$(head -n1 "$file" 2>/dev/null || true)"
  data="$(sed -n '2p' "$file" 2>/dev/null || true)"
  if [[ -z "$header" || -z "$data" ]]; then
    echo ""
    return 0
  fi
  IFS=',' read -r -a headers <<< "$header"
  IFS=',' read -r -a values <<< "$data"
  for key in "$@"; do
    for i in "${!headers[@]}"; do
      if [[ "${headers[$i]}" == "$key" ]]; then
        echo "${values[$i]:-}"
        return 0
      fi
    done
  done
  echo ""
}

# Example single setting:
# STELARX_OPTS="--search-space S2 --intersection-method I2"
#
# Example sweep over search-space presets:
# STELARX_OPTS_LIST_RAW="--search-space S1;--search-space S2;--search-space S3"
#
# This becomes search-space_S2__intersection-method_I2. Verbosity flags such as
# -v/-vv are ignored when constructing the setting name.

print_help() {
  cat <<EOF
run-a10k.sh

Required:
  --data-dir           Path to A10K dataset root containing 10k-simphy/

Optional:
  --tree-type          estimated | true, or a semicolon-separated list
                       such as "true;estimated" (default: estimated)
  --replicates         Replicates to run, e.g. "1-20" or "R1,R2"
  --start-rep, -sr     Start replicate number
  --end-rep, -er       End replicate number
  --stelarx-root       Path to STELAR-X root
  --opts, --alg-opts   Extra options for the selected algorithm
  --opts-list, --alg-opts-list
                       Semicolon-separated list of option strings to loop over
  --fresh              Force rerun even if stat-stelarx.csv exists
  --a10k-outputs-dir   Reproducibility mirror root for the small run outputs
                       (default: the "outputs" sibling of the dataset, i.e.
                        ".../10k-astral-dataset" -> ".../outputs/10k-astral-dataset")
  --no-outputs-mirror  Do not copy results into the outputs mirror
  --no-time-monitor    Disable time monitoring
  --no-gpu-monitor     Disable GPU monitoring
  --no-notify, -nn     Disable ntfy notifications

Examples:
  ./run-a10k.sh --data-dir /path/to/10k-astral-dataset --tree-type estimated --opts "--search-space S1 --intersection-method I2 -vv"
  ./run-a10k.sh --data-dir /path/to/10k-astral-dataset --tree-type "true;estimated" --opts "--search-space S1 --intersection-method I2 -vv"
  ./run-a10k.sh --data-dir /path/to/10k-astral-dataset --tree-type estimated --opts "--search-space S2 --intersection-method I2 -vv"
  ./run-a10k.sh --data-dir /path/to/10k-astral-dataset --tree-type estimated --opts-list "--search-space S1 -vv;--search-space S2 -vv;--search-space S3 -vv"
  The first example setting is search-space_S1__intersection-method_I2.
  Verbosity is ignored; other meaningful options are appended to the name.

Results are written to
  <data-dir>/10k-simphy/<R>/stelarx_outputs/<tree-type>/<setting>
as before and mirrored (with the dataset record and the replicate's input
fingerprints) to
  <outputs>/stelarx_outputs/<R>/<tree-type>/<setting>.
Gene trees and species trees are never copied into the mirror.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --tree-type) TREE_TYPES_RAW="$2"; shift 2 ;;
    --tree-type=*) TREE_TYPES_RAW="${1#*=}"; shift ;;
    --replicates) REPLICATES_SPEC="$2"; shift 2 ;;
    --start-rep|-sr) START_REP="$2"; shift 2 ;;
    --end-rep|-er) END_REP="$2"; shift 2 ;;
    --stelarx-root|--stelar-root) STELARX_ROOT="$2"; shift 2 ;;
    --opts|--alg-opts|--stelarx-opts|--stelar-opts) STELARX_OPTS="$2"; shift 2 ;;
    --opts-list|--alg-opts-list|--stelarx-opts-list|--stelar-opts-list) STELARX_OPTS_LIST_RAW="$2"; shift 2 ;;
    --fresh) FRESH=true; shift ;;
    --a10k-outputs-dir|--outputs-dir) A10K_OUTPUTS_DIR="$2"; shift 2 ;;
    --a10k-outputs-dir=*|--outputs-dir=*) A10K_OUTPUTS_DIR="${1#*=}"; shift ;;
    --no-outputs-mirror) OUTPUTS_MIRROR=false; shift ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ -z "$DATA_DIR" ]]; then
  echo "Error: --data-dir is required."
  exit 2
fi

STELARX_ROOT="$(realpath "$STELARX_ROOT")"
PYTHON_BIN="${STELARX_PYTHON:-${STELARX_ROOT}/.venv/bin/python}"
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN="python3"
DATA_DIR="$(realpath "$DATA_DIR")"
SIMPHY_DIR="${DATA_DIR%/}/10k-simphy"
if [[ ! -d "$SIMPHY_DIR" ]]; then
  echo "Error: expected 10k-simphy at $SIMPHY_DIR"
  exit 3
fi

if [[ "$OUTPUTS_MIRROR" == true ]]; then
  A10K_OUTPUTS_DIR="$(stelarx_prepare_a10k_outputs_dir "$A10K_OUTPUTS_DIR" "$DATA_DIR")" || exit 2
else
  A10K_OUTPUTS_DIR="(disabled)"
fi

# Copy one results directory into the reproducibility mirror. A mirror problem
# is reported loudly but never changes the outcome of the run itself.
mirror_results_dir() {
  local results_dir="$1" mirrored
  [[ "$OUTPUTS_MIRROR" == true ]] || return 0
  [[ -d "$results_dir" ]] || return 0
  if mirrored="$(stelarx_mirror_a10k_results "$DATA_DIR" "$A10K_OUTPUTS_DIR" "$results_dir")"; then
    echo "Mirrored outputs to: $mirrored"
  else
    echo "WARNING: outputs mirror was not updated for $results_dir" >&2
  fi
}

# Complete the command record (out-stelarx.command, written by the wrapper with
# the exact run.sh invocation) with the outer commands that produced this run.
# Reads the loop variables of the case currently being run.
append_a10k_command_context() {
  local rf_rate="$1"
  [[ -f "$COMMAND_FILE" ]] || return 0
  {
    echo "# --- A10K run context (run-a10k.sh) ---"
    echo "# dataset:      $(basename "$DATA_DIR")"
    echo "# replicate:    $REPL"
    echo "# tree type:    $TREE_TYPE"
    echo "# setting:      $SETTING_NAME"
    echo "# gene trees:   $GT_FILE"
    echo "# true tree:    $TRUE_TREE"
    echo "# rf_rate:      $rf_rate"
    printf '# invoked as: '
    printf ' %q' "${SCRIPT_ARGV[@]}"
    printf '\n'
    printf '# wrapper cmd:'
    printf ' %q' "${CMD[@]}"
    printf '\n'
  } >> "$COMMAND_FILE" 2>/dev/null || echo "Warning: could not append to command record $COMMAND_FILE" >&2
}

TREE_TYPES=()
IFS=';' read -r -a raw_tree_types <<< "$TREE_TYPES_RAW"
for tree_type in "${raw_tree_types[@]}"; do
  tree_type="${tree_type//[[:space:]]/}"
  tree_type="${tree_type,,}"
  [[ -n "$tree_type" ]] || continue
  case "$tree_type" in
    true|estimated) ;;
    *)
      echo "Error: invalid --tree-type value '$tree_type' (expected true, estimated, or a semicolon-separated list)."
      exit 2
      ;;
  esac

  duplicate=false
  for existing_tree_type in "${TREE_TYPES[@]}"; do
    if [[ "$existing_tree_type" == "$tree_type" ]]; then
      duplicate=true
      break
    fi
  done
  [[ "$duplicate" == false ]] && TREE_TYPES+=("$tree_type")
done
if [[ ${#TREE_TYPES[@]} -eq 0 ]]; then
  echo "Error: --tree-type must contain at least one of: true, estimated."
  exit 2
fi

STELARX_OPTS_LIST=()
if [[ -n "$STELARX_OPTS_LIST_RAW" ]]; then
  IFS=';' read -r -a raw_opts_list <<< "$STELARX_OPTS_LIST_RAW"
  for opts in "${raw_opts_list[@]}"; do
    opts="$(echo "$opts" | sed 's/^ *//;s/ *$//')"
    [[ -n "$opts" ]] && STELARX_OPTS_LIST+=("$opts")
  done
fi
if [[ ${#STELARX_OPTS_LIST[@]} -eq 0 ]]; then
  STELARX_OPTS_LIST+=("${STELARX_OPTS}")
fi
echo "[DEBUG] opts list (${#STELARX_OPTS_LIST[@]} items): ${STELARX_OPTS_LIST[*]}"
echo "[DEBUG] tree types (${#TREE_TYPES[@]} items): ${TREE_TYPES[*]}"
echo "[DEBUG] replicates spec: '${REPLICATES_SPEC}' | fresh: ${FRESH}"
echo "[DEBUG] outputs mirror: ${A10K_OUTPUTS_DIR}"

REPL_LIST=()
if [[ -n "$START_REP" || -n "$END_REP" ]]; then
  for i in $(seq "$START_REP" "$END_REP"); do REPL_LIST+=("R${i}"); done
elif [[ -n "$REPLICATES_SPEC" ]]; then
  if [[ "$REPLICATES_SPEC" =~ ^[0-9]+-[0-9]+$ ]]; then
    start="${REPLICATES_SPEC%-*}"
    end="${REPLICATES_SPEC#*-}"
    for i in $(seq "$start" "$end"); do REPL_LIST+=("R${i}"); done
  else
    IFS=',' read -r -a parts <<< "$REPLICATES_SPEC"
    for p in "${parts[@]}"; do
      p="${p// /}"
      [[ "$p" =~ ^R ]] || p="R${p}"
      REPL_LIST+=("$p")
    done
  fi
else
  while IFS= read -r -d '' d; do
    REPL_LIST+=("$(basename "$d")")
  done < <(find "$SIMPHY_DIR" -maxdepth 1 -type d -name 'R*' -print0 | sort -z -V)
fi

echo "[DEBUG] replicate list (${#REPL_LIST[@]} items): ${REPL_LIST[*]}"

for TREE_TYPE in "${TREE_TYPES[@]}"; do
  echo "==> Processing A10K tree type: ${TREE_TYPE}"
  for REPL in "${REPL_LIST[@]}"; do
  REPL_DIR="${SIMPHY_DIR%/}/${REPL}"
  if [[ ! -d "$REPL_DIR" ]]; then
    echo "[DEBUG] SKIP ${REPL}: directory not found: ${REPL_DIR}"
    continue
  fi

  if [[ "$TREE_TYPE" == "estimated" ]]; then
    GT_DIR="${REPL_DIR}/estimatedgenetrees"
    GT_FILE="${GT_DIR}/estimatedgenetrees.tre"
    ROOTED_GT="${GT_DIR}/estimatedgenetrees.rooted.tre"
    if [[ ! -f "$ROOTED_GT" ]]; then
      if [[ ! -x "${STELARX_ROOT%/}/process_unrooted.sh" ]]; then
        echo "Error: process_unrooted.sh not found or not executable at ${STELARX_ROOT%/}/process_unrooted.sh"
        exit 7
      fi
      echo "Rooting estimated gene trees for ${REPL} with outgroup 0..."
      "${STELARX_ROOT%/}/process_unrooted.sh" -i "$GT_FILE" -o "$ROOTED_GT" -og "0"
    fi
    GT_FILE="$ROOTED_GT"
  else
    GT_FILE="${REPL_DIR}/truegenetrees"
  fi
  TRUE_TREE="${REPL_DIR}/s_tree.trees"
  if [[ ! -f "$GT_FILE" || ! -f "$TRUE_TREE" ]]; then
    echo "[DEBUG] SKIP ${REPL}: missing files (gt=${GT_FILE} exists=$([ -f "$GT_FILE" ] && echo yes || echo no), true_tree=${TRUE_TREE} exists=$([ -f "$TRUE_TREE" ] && echo yes || echo no))"
    continue
  fi

  for STELARX_OPTS_ITEM in "${STELARX_OPTS_LIST[@]}"; do
    SETTING_NAME="$(build_setting_name_from_opts "$STELARX_OPTS_ITEM")"
    OUT_DIR="${REPL_DIR}/stelarx_outputs/${TREE_TYPE}/${SETTING_NAME}"
    OUT_FILE="${OUT_DIR}/out-stelarx.tre"
    STAT_FILE="${OUT_DIR}/stat-stelarx.csv"
    RUN_LOG="${OUT_DIR}/.stelarx_run.log"
    COMMAND_FILE="${OUT_FILE%.*}.command"

    if [[ "$FRESH" == false && -f "$STAT_FILE" ]]; then
      echo "SKIPPING: ${STAT_FILE} exists."
      # Keep the reproducibility mirror complete even for runs finished earlier.
      mirror_results_dir "$OUT_DIR"
      continue
    elif [[ "$FRESH" == true && -f "$STAT_FILE" ]]; then
      echo "[DEBUG] --fresh set, overwriting existing: ${STAT_FILE}"
    fi

    mkdir -p "$OUT_DIR"
    rm -f "$RUN_LOG"
    CMD=("${STELARX_ROOT}/run-stelarx-with-monitor.sh" -i "$GT_FILE" -o "$OUT_FILE" --stelarx-root "$STELARX_ROOT")
    if [[ "$TIME_MONITOR" == false ]]; then CMD+=(--no-time-monitor); fi
    if [[ "$GPU_MONITOR" == false ]]; then CMD+=(--no-gpu-monitor); fi
    if [[ "$NO_NOTIFY" == true ]]; then CMD+=(--no-notify); fi
    if [[ -n "$STELARX_OPTS_ITEM" ]]; then
      CMD+=(--opts "$STELARX_OPTS_ITEM")
    fi

    echo "==> Running stelarx on ${REPL} (${TREE_TYPE}, ${SETTING_NAME})"
    echo "Command: ${CMD[*]}"
    set +e
    "${CMD[@]}" 2>&1 | tee "$RUN_LOG"
    RUN_EXIT=${PIPESTATUS[0]}
    set -e

    SIDE_STATS="${OUT_FILE%.tre}_stats.csv"
    if [[ "$RUN_EXIT" -ne 0 || ! -f "$SIDE_STATS" ]]; then
      echo "Run failed for ${REPL} (${TREE_TYPE}, ${SETTING_NAME}); skipping RF/stat summary."
      # The failed run's log and command record are still worth keeping.
      append_a10k_command_context "NA"
      mirror_results_dir "$OUT_DIR"
      continue
    fi

    RUNNING_TIME="$(csv_get_field "$SIDE_STATS" "running_time_s" "running-time-s")"
    MAX_CPU_MB="$(csv_get_field "$SIDE_STATS" "max_cpu_mb" "max-cpu-mb")"
    MAX_GPU_MB="$(csv_get_field "$SIDE_STATS" "max_gpu_mb" "max-gpu-mb")"
    OPTIMAL_TRIPLET_SCORE="$(csv_get_field "$SIDE_STATS" "optimal_triplet_score" "optimal-triplet-score")"
    EXIT_CODE="$(csv_get_field "$SIDE_STATS" "exit_code" "exit-code")"
    if [[ -z "$EXIT_CODE" ]]; then
      EXIT_CODE="$RUN_EXIT"
    fi

    RF_RATE="NA"
    if [[ -f "$OUT_FILE" && -f "$TRUE_TREE" ]]; then
      rf_output=$("$PYTHON_BIN" "${STELARX_ROOT}/rf.py" "$OUT_FILE" "$TRUE_TREE" 2>&1) || true
      rf_line=$(echo "$rf_output" | grep -i "Robinson-Foulds distance" | tail -n1 || true)
      if [[ -n "$rf_line" ]]; then
        RF_RATE=$(echo "$rf_line" | grep -Eo '[0-9]+(\.[0-9]+)?' | tail -n1 || echo "NA")
      fi
    fi

    echo "alg,setting,replicate,tree_type,rf-rate,optimal-triplet-score,running-time-s,max-cpu-mb,max-gpu-mb" > "$STAT_FILE"
    echo "stelarx,${SETTING_NAME},${REPL},${TREE_TYPE},${RF_RATE},${OPTIMAL_TRIPLET_SCORE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB}" >> "$STAT_FILE"
    echo
    echo "=== A10K STELAR-X Summary ==="
    echo "Replicate:      ${REPL}"
    echo "Tree type:      ${TREE_TYPE}"
    echo "Setting:        ${SETTING_NAME}"
    echo "RF rate:        ${RF_RATE}"
    echo "Triplet score:  ${OPTIMAL_TRIPLET_SCORE}"
    echo "Running time:   ${RUNNING_TIME}s"
    echo "Max CPU RAM:    ${MAX_CPU_MB} MB"
    echo "Max GPU VRAM:   ${MAX_GPU_MB} MB"
    echo "Output tree:    ${OUT_FILE}"
    echo "Stats file:     ${STAT_FILE}"
    echo "Saved $STAT_FILE"

    append_a10k_command_context "$RF_RATE"
    mirror_results_dir "$OUT_DIR"

    if [[ "$NO_NOTIFY" == false ]] && command -v curl >/dev/null 2>&1; then
      curl -s -d "✅ STELAR-X A10K completed

Replicate: ${REPL}
Tree type: ${TREE_TYPE}
Setting: ${SETTING_NAME}

RF: ${RF_RATE}
Triplet score: ${OPTIMAL_TRIPLET_SCORE}
Time: ${RUNNING_TIME}s
CPU: ${MAX_CPU_MB} MB
GPU: ${MAX_GPU_MB} MB
Exit: ${EXIT_CODE}

Tree: $(basename "$OUT_FILE")
Stats: $(basename "$STAT_FILE")" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
    fi
    done
  done
done
