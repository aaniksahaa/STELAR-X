#!/usr/bin/env bash
# run-a10k.sh
# Run baselines on the 10k-astral dataset directory structure.
# For each replicate, run a baseline on either estimated or true gene trees,
# store outputs under <R?>/<method>_outputs/{estimated,true}, and compute RF.

set -euo pipefail

NTFY_CHANNEL_NAME="anik-test"

METHOD=""
TREE_TYPE="estimated"   # estimated | true
DATA_DIR=""
BASE_DIR=""
BASE_DIR_SET=false
STELAR_ROOT=""
STELAR_ROOT_SET=false
BASELINES_DIR=""
BASELINES_DIR_SET=false
REPLICATES_SPEC=""
FRESH=false

# Method options
ASTER_OPTS=""
ASTER_BIN=""
ASTRAL_OPTS=""
ASTRAL_XMS=""
ASTRAL_XMX=""
TREEQMC_OPTS=""
WQFM_OPTS=""
SUPERTRIPLETS_OPTS=""
TMC_OPTS=""
STELAR_OPTS=""

# Monitoring options (DEFAULT: ON)
TIME_MONITOR=true
GPU_MONITOR=true
NO_NOTIFY=false
DEBUG=0

print_help() {
  cat <<EOF
run-a10k.sh

Required:
  --method, -m         Method: stelar, aster, astral, treeqmc, wqfmtree, supertriplets, tmc
  --data-dir           Path to dataset root (e.g., /home/user/phylogeny/datasets/10k-astral-dataset)

Optional:
  --tree-type          estimated | true (default: estimated)
  --replicates         Replicates to run (e.g., "1-20" or "R1,R2,R3"; default: all R* under 10k-simphy)
  --fresh              Force rerun even if stat file exists
  --base-dir, -b       Base dir (assumes STELAR-X at BASE_DIR/STELAR-X)
  --stelar-root        Path to STELAR-X root (overrides --base-dir)
  --baselines-dir      Path to baselines directory (overrides --stelar-root)
  --no-time-monitor    Disable time-monitoring
  --no-gpu-monitor     Disable GPU-monitoring
  --no-notify, -nn     Disable ntfy.sh notifications
  --debug              Enable shell tracing

Method Options:
  --aster-opts "..."         Extra options for ASTER
  --aster-bin PATH          Path to ASTER binary (relative to ASTER root or absolute)
  --astral-opts "..."        Extra options for ASTRAL
  --astral-xms SIZE          ASTRAL Java Xms (e.g., 4g)
  --astral-xmx SIZE          ASTRAL Java Xmx (e.g., 128g)
  --astral-Xms SIZE          Same as --astral-xms
  --astral-Xmx SIZE          Same as --astral-xmx
  --treeqmc-opts "..."       Extra options for TreeQMC
  --wqfm-opts "..."          Extra options for wQFM-TREE
  --supertriplets-opts "..." Extra options for SuperTriplets wrapper
  --tmc-opts "..."           Extra options for TMC wrapper
  --stelar-opts "..."        Extra options for STELAR-X (passed to run-with-monitor.sh)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --method|-m) METHOD="$2"; shift 2 ;;
    --tree-type) TREE_TYPE="$2"; shift 2 ;;
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --replicates) REPLICATES_SPEC="$2"; shift 2 ;;
    --fresh) FRESH=true; shift ;;
    --base-dir|-b) BASE_DIR="$2"; BASE_DIR_SET=true; shift 2 ;;
    --stelar-root) STELAR_ROOT="$2"; STELAR_ROOT_SET=true; shift 2 ;;
    --baselines-dir) BASELINES_DIR="$2"; BASELINES_DIR_SET=true; shift 2 ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
    --aster-opts) ASTER_OPTS="$2"; shift 2 ;;
    --aster-opts=*) ASTER_OPTS="${1#*=}"; shift ;;
    --aster-bin) ASTER_BIN="$2"; shift 2 ;;
    --astral-opts) ASTRAL_OPTS="$2"; shift 2 ;;
    --astral-opts=*) ASTRAL_OPTS="${1#*=}"; shift ;;
    --astral-xms|--astral-Xms) ASTRAL_XMS="$2"; shift 2 ;;
    --astral-xmx|--astral-Xmx) ASTRAL_XMX="$2"; shift 2 ;;
    --treeqmc-opts) TREEQMC_OPTS="$2"; shift 2 ;;
    --treeqmc-opts=*) TREEQMC_OPTS="${1#*=}"; shift ;;
    --wqfm-opts) WQFM_OPTS="$2"; shift 2 ;;
    --wqfm-opts=*) WQFM_OPTS="${1#*=}"; shift ;;
    --supertriplets-opts) SUPERTRIPLETS_OPTS="$2"; shift 2 ;;
    --supertriplets-opts=*) SUPERTRIPLETS_OPTS="${1#*=}"; shift ;;
    --tmc-opts) TMC_OPTS="$2"; shift 2 ;;
    --tmc-opts=*) TMC_OPTS="${1#*=}"; shift ;;
    --stelar-opts) STELAR_OPTS="$2"; shift 2 ;;
    --stelar-opts=*) STELAR_OPTS="${1#*=}"; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$METHOD" || -z "$DATA_DIR" ]]; then
  echo "Error: --method and --data-dir are required."
  print_help
  exit 2
fi

case "$METHOD" in
  stelar|aster|astral|treeqmc|tree-qmc|wqfmtree|wqfm-tree|supertriplets|tmc) ;;
  *) echo "Error: unknown method '$METHOD'"; exit 3 ;;
esac
if [[ "$METHOD" == "tree-qmc" ]]; then
  METHOD="treeqmc"
fi
if [[ "$METHOD" == "wqfm-tree" ]]; then
  METHOD="wqfmtree"
fi

if [[ "$TREE_TYPE" != "estimated" && "$TREE_TYPE" != "true" ]]; then
  echo "Error: --tree-type must be 'estimated' or 'true'"
  exit 4
fi

if [[ "$STELAR_ROOT_SET" = false ]]; then
  if [[ "$BASE_DIR_SET" = true ]]; then
    STELAR_ROOT="${BASE_DIR%/}/STELAR-X"
  else
    STELAR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  fi
fi

if [[ "$BASELINES_DIR_SET" = false ]]; then
  BASELINES_DIR="${STELAR_ROOT%/}/baselines"
fi

DATA_DIR=$(realpath "$DATA_DIR")
STELAR_ROOT=$(realpath "$STELAR_ROOT")
BASELINES_DIR=$(realpath "$BASELINES_DIR")

if [[ "${DEBUG:-0}" = "1" ]]; then
  set -x
fi

SIMPHY_DIR="${DATA_DIR%/}/10k-simphy"
if [[ ! -d "$SIMPHY_DIR" ]]; then
  echo "Error: expected simphy dir at $SIMPHY_DIR"
  exit 5
fi

if [[ -n "$REPLICATES_SPEC" ]]; then
  REPL_LIST=()
  if [[ "$REPLICATES_SPEC" =~ ^[0-9]+-[0-9]+$ ]]; then
    start="${REPLICATES_SPEC%-*}"
    end="${REPLICATES_SPEC#*-}"
    for i in $(seq "$start" "$end"); do
      REPL_LIST+=("R${i}")
    done
  else
    IFS=',' read -r -a parts <<< "$REPLICATES_SPEC"
    for p in "${parts[@]}"; do
      p_trim="$(echo "$p" | sed 's/^ *//;s/ *$//')"
      if [[ "$p_trim" =~ ^R[0-9]+$ ]]; then
        REPL_LIST+=("$p_trim")
      elif [[ "$p_trim" =~ ^[0-9]+$ ]]; then
        REPL_LIST+=("R${p_trim}")
      fi
    done
  fi
else
  REPL_LIST=()
  while IFS= read -r -d '' d; do
    REPL_LIST+=("$(basename "$d")")
  done < <(find "$SIMPHY_DIR" -maxdepth 1 -type d -name 'R*' -print0 | sort -z)
fi

if [[ ${#REPL_LIST[@]} -eq 0 ]]; then
  echo "Error: no replicates found."
  exit 6
fi

for REPL in "${REPL_LIST[@]}"; do
  REPL_DIR="${SIMPHY_DIR%/}/${REPL}"
  if [[ ! -d "$REPL_DIR" ]]; then
    echo "Skipping missing replicate dir: $REPL_DIR"
    continue
  fi

  if [[ "$TREE_TYPE" == "estimated" ]]; then
    GT_DIR="${REPL_DIR%/}/estimatedgenetrees"
  else
    GT_DIR="${REPL_DIR%/}/truegenetrees"
  fi

  GT_FILE=""
  if [[ -f "${GT_DIR%/}/estimatedgenetrees.tre" ]]; then
    GT_FILE="${GT_DIR%/}/estimatedgenetrees.tre"
  elif [[ -f "${GT_DIR%/}/truegenetrees.tre" ]]; then
    GT_FILE="${GT_DIR%/}/truegenetrees.tre"
  else
    GT_FILE="$(find "$GT_DIR" -maxdepth 1 -type f -name '*.tre' | head -n1 || true)"
  fi

  TRUE_TREE="${REPL_DIR%/}/s_tree.trees"
  if [[ -z "$GT_FILE" || ! -f "$GT_FILE" ]]; then
    echo "Skipping ${REPL}: gene trees not found in $GT_DIR"
    continue
  fi
  if [[ ! -f "$TRUE_TREE" ]]; then
    echo "Skipping ${REPL}: species tree not found at $TRUE_TREE"
    continue
  fi

  # For STELAR on estimated trees, root unrooted gene trees with outgroup "0"
  ROOTED_GT=""
  if [[ "$METHOD" == "stelar" && "$TREE_TYPE" == "estimated" ]]; then
    ROOTED_GT="${GT_DIR%/}/estimatedgenetrees.rooted.tre"
    if [[ "$FRESH" = true || ! -f "$ROOTED_GT" ]]; then
      if [[ ! -x "${STELAR_ROOT%/}/process_unrooted.sh" ]]; then
        echo "Error: process_unrooted.sh not found or not executable at ${STELAR_ROOT%/}/process_unrooted.sh"
        exit 7
      fi
      echo "Rooting estimated gene trees for ${REPL} with outgroup 0..."
      "${STELAR_ROOT%/}/process_unrooted.sh" -i "$GT_FILE" -o "$ROOTED_GT" -og "0"
    fi
    GT_FILE="$ROOTED_GT"
  fi

  OUT_DIR="${REPL_DIR%/}/${METHOD}_outputs/${TREE_TYPE}"
  mkdir -p "$OUT_DIR"
  OUT_FILE="${OUT_DIR%/}/out-${METHOD}.tre"
  STAT_FILE="${OUT_DIR%/}/stat-${METHOD}.csv"

  if [[ "$FRESH" = false && -f "$STAT_FILE" ]]; then
    echo "SKIPPING: ${STAT_FILE} exists. Use --fresh to rerun."
    continue
  fi

  # Choose int128 for ASTER if not provided and taxa is large
  if [[ "$METHOD" == "aster" && -z "$ASTER_BIN" ]]; then
    ASTER_BIN="bin/astral4_int128"
  fi

  if [[ "$METHOD" == "stelar" ]]; then
    CMD=("${STELAR_ROOT%/}/run-with-monitor.sh" -i "$GT_FILE" -o "$OUT_FILE" --stelar-root "$STELAR_ROOT")
    if [[ "$TIME_MONITOR" = false ]]; then CMD+=(--no-time-monitor); fi
    if [[ "$GPU_MONITOR" = false ]]; then CMD+=(--no-gpu-monitor); fi
    if [[ "$NO_NOTIFY" = true ]]; then CMD+=(--no-notify); fi
    if [[ "$DEBUG" = 1 ]]; then CMD+=(--debug); fi
    if [[ -n "$STELAR_OPTS" ]]; then CMD+=($STELAR_OPTS); fi
  else
    CMD=("${STELAR_ROOT%/}/run-baseline-with-monitor.sh" -m "$METHOD" -i "$GT_FILE" -o "$OUT_FILE" --stelar-root "$STELAR_ROOT" --baselines-dir "$BASELINES_DIR")
    if [[ "$TIME_MONITOR" = false ]]; then CMD+=(--no-time-monitor); fi
    if [[ "$GPU_MONITOR" = false ]]; then CMD+=(--no-gpu-monitor); fi
    if [[ "$NO_NOTIFY" = true ]]; then CMD+=(--no-notify); fi
    if [[ "$DEBUG" = 1 ]]; then CMD+=(--debug); fi

    if [[ -n "$ASTER_OPTS" ]]; then CMD+=(--aster-opts "$ASTER_OPTS"); fi
    if [[ -n "$ASTER_BIN" ]]; then CMD+=(--aster-bin "$ASTER_BIN"); fi
  if [[ -n "$ASTRAL_OPTS" ]]; then CMD+=(--astral-opts "$ASTRAL_OPTS"); fi
  if [[ -n "$ASTRAL_XMS" ]]; then CMD+=(--astral-xms "$ASTRAL_XMS"); fi
  if [[ -n "$ASTRAL_XMX" ]]; then CMD+=(--astral-xmx "$ASTRAL_XMX"); fi
    if [[ -n "$TREEQMC_OPTS" ]]; then CMD+=(--treeqmc-opts "$TREEQMC_OPTS"); fi
    if [[ -n "$WQFM_OPTS" ]]; then CMD+=(--wqfm-opts "$WQFM_OPTS"); fi
    if [[ -n "$SUPERTRIPLETS_OPTS" ]]; then CMD+=(--supertriplets-opts "$SUPERTRIPLETS_OPTS"); fi
    if [[ -n "$TMC_OPTS" ]]; then CMD+=(--tmc-opts "$TMC_OPTS"); fi
  fi

  echo "==> Running ${METHOD} on ${REPL} (${TREE_TYPE})"
  echo "Command: ${CMD[*]}"
  set +e
  "${CMD[@]}"
  RUN_EXIT=$?
  set -e

  STATS_FROM_BASELINE="${OUT_FILE%.tre}_stats.csv"
  if [[ "$RUN_EXIT" -ne 0 || ! -f "$STATS_FROM_BASELINE" ]]; then
    echo "Run failed for ${REPL} (exit: ${RUN_EXIT}). Skipping RF/stat."
    continue
  fi

  RUNNING_TIME=$(awk -F, 'NR==2{print $4}' "$STATS_FROM_BASELINE")
  MAX_CPU_MB=$(awk -F, 'NR==2{print $5}' "$STATS_FROM_BASELINE")
  MAX_GPU_MB=$(awk -F, 'NR==2{print $6}' "$STATS_FROM_BASELINE")
  EXIT_CODE=$(awk -F, 'NR==2{print $7}' "$STATS_FROM_BASELINE")

  if [[ "$EXIT_CODE" != "0" ]]; then
    echo "Run failed for ${REPL} (stats exit: ${EXIT_CODE}). Skipping RF/stat."
    continue
  fi

  RF_RATE="NA"
  if [[ -f "${STELAR_ROOT%/}/rf.py" && -x "$(command -v python3)" ]]; then
    rf_output=$(python3 "${STELAR_ROOT%/}/rf.py" "$OUT_FILE" "$TRUE_TREE" 2>&1) || true
    rf_candidate=$(echo "$rf_output" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
    if [[ -n "$rf_candidate" ]]; then
      RF_RATE="$rf_candidate"
    fi
  fi

  echo "alg,replicate,tree_type,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb,exit-code" > "$STAT_FILE"
  echo "${METHOD},${REPL},${TREE_TYPE},${RF_RATE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${EXIT_CODE}" >> "$STAT_FILE"
  echo "Wrote stats to $STAT_FILE"

  if [[ "$NO_NOTIFY" != true && -n "$NTFY_CHANNEL_NAME" && -x "$(command -v curl)" ]]; then
    NTFY_MSG="✅ ${METHOD} ${REPL} (${TREE_TYPE})\n\nRF: ${RF_RATE}\nTime: ${RUNNING_TIME}s\nCPU: ${MAX_CPU_MB} MB | GPU: ${MAX_GPU_MB} MB\n\n${STAT_FILE}"
    curl -s -d "$NTFY_MSG" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
  fi
done

echo "All runs finished."
