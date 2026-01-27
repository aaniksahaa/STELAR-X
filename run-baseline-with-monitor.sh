#!/usr/bin/env bash
# run-baseline-with-monitor.sh
# Unified baseline runner with time/GPU monitoring
# Supports: aster, astral, treeqmc, wqfmtree, supertriplets, tmc

set -euo pipefail

# ===== CONFIGURATION =====
NTFY_CHANNEL_NAME="anik-test"

# Defaults
METHOD=""
INPUT_FILE=""
OUTPUT_FILE=""
BASE_DIR=""
BASE_DIR_SET=false
STELAR_ROOT=""
STELAR_ROOT_SET=false
BASELINES_DIR=""
BASELINES_DIR_SET=false

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

# Monitoring options (DEFAULT: ON)
TIME_MONITOR=true
GPU_MONITOR=true
NO_NOTIFY=false
DEBUG=0

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_help() {
  cat <<EOF
run-baseline-with-monitor.sh - Unified baseline wrapper with performance monitoring

Usage:
  \$0 -m <method> -i <input_file> -o <output_file> [options]

Required:
  -m, --method METHOD     Baseline method: aster, astral, treeqmc, wqfmtree, supertriplets, tmc
  -i, --input FILE        Path to gene trees file
  -o, --output FILE       Path to output species tree file

Path Options:
  --base-dir, -b DIR      Base directory (assumes STELAR-X is in DIR/STELAR-X)
  --stelar-root DIR       Path to STELAR-X root (overrides --base-dir)
  --baselines-dir DIR     Path to baselines directory (overrides --stelar-root)

Method Options:
  --aster-opts "..."         Extra options for ASTER
  --aster-bin PATH          Path to ASTER binary (relative to ASTER root or absolute)
  --astral-opts "..."        Extra options for ASTRAL
  --astral-xms SIZE          ASTRAL Java Xms (e.g., 4g)
  --astral-xmx SIZE          ASTRAL Java Xmx (e.g., 128g)
  --astral-Xms SIZE          Same as --astral-xms
  --astral-Xmx SIZE          Same as --astral-xmx
  --treeqmc-opts "..."       Extra options for TreeQMC
  --wqfm-opts "..."          Extra options for wQFM-TREE (if supported)
  --supertriplets-opts "..." Extra options for SuperTriplets wrapper
  --tmc-opts "..."           Extra options for TMC wrapper

Monitoring Options:
  --no-time-monitor        Disable time-monitoring
  --no-gpu-monitor         Disable GPU-monitoring
  --no-notify, -nn         Disable ntfy.sh notifications
  --debug                  Enable shell tracing

Examples:
  \$0 -m aster -i genes.tre -o out.tre --aster-opts "-t 16"
  \$0 -m astral -i genes.tre -o out.tre --astral-opts "-t 8"
  \$0 -m treeqmc -i genes.tre -o out.tre
  \$0 -m wqfmtree -i genes.tre -o out.tre
  \$0 -m supertriplets -i genes.tre -o out.tre
  \$0 -m tmc -i genes.tre -o out.tre --tmc-opts "--nd 1000"
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--method) METHOD="$2"; shift 2 ;;
    -i|--input) INPUT_FILE="$2"; shift 2 ;;
    -o|--output) OUTPUT_FILE="$2"; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; BASE_DIR_SET=true; shift 2 ;;
    --stelar-root) STELAR_ROOT="$2"; STELAR_ROOT_SET=true; shift 2 ;;
    --baselines-dir) BASELINES_DIR="$2"; BASELINES_DIR_SET=true; shift 2 ;;
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
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
    -h|--help) print_help; exit 0 ;;
    *)
      echo -e "${RED}Error: Unknown option: $1${NC}" >&2
      print_help
      exit 1
      ;;
  esac
done

if [[ -z "$METHOD" || -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
  echo -e "${RED}Error: --method, --input, and --output are required.${NC}" >&2
  print_help
  exit 2
fi

# Normalize method name
case "$METHOD" in
  aster|astral|treeqmc|tree-qmc|wqfmtree|wqfm-tree|supertriplets|tmc) ;;
  *)
    echo -e "${RED}Error: Unknown method '$METHOD'.${NC}" >&2
    exit 3
    ;;
esac

if [[ "$METHOD" == "tree-qmc" ]]; then
  METHOD="treeqmc"
fi
if [[ "$METHOD" == "wqfm-tree" ]]; then
  METHOD="wqfmtree"
fi

# Derive STELAR_ROOT
if [[ "$STELAR_ROOT_SET" = false ]]; then
  if [[ "$BASE_DIR_SET" = true ]]; then
    STELAR_ROOT="${BASE_DIR%/}/STELAR-X"
  else
    STELAR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  fi
fi

# Derive BASELINES_DIR
if [[ "$BASELINES_DIR_SET" = false ]]; then
  BASELINES_DIR="${STELAR_ROOT%/}/baselines"
fi

# Convert to absolute paths
INPUT_FILE=$(realpath "$INPUT_FILE")
OUTPUT_FILE=$(realpath "$OUTPUT_FILE")
STELAR_ROOT=$(realpath "$STELAR_ROOT")
BASELINES_DIR=$(realpath "$BASELINES_DIR")

# Debug/tracing
if [[ "${DEBUG:-0}" = "1" ]]; then
  set -x
fi

# Validate input file
if [[ ! -f "$INPUT_FILE" ]]; then
  echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}" >&2
  exit 4
fi

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Create temporary directory for logs
TEMP_DIR=$(mktemp -d)
TIME_TMP="$TEMP_DIR/baseline_time_err.log"
MON_TMP="$TEMP_DIR/baseline_gpu_mem.log"

cleanup() {
  rm -rf "$TEMP_DIR" 2>/dev/null || true
  if [[ -n "${MON_PID:-}" ]]; then
    kill "$MON_PID" 2>/dev/null || true
    wait "$MON_PID" 2>/dev/null || true
  fi
}
trap cleanup EXIT

# START timer
START_NS=$(date +%s%N)

# Detect usable 'time' command
TIME_CMD=""
if [[ "${TIME_MONITOR:-false}" = true ]]; then
  if [[ -x "/usr/bin/time" ]]; then
    TIME_CMD="/usr/bin/time"
  else
    if command -v time >/dev/null 2>&1; then
      TMP_TEST="$(mktemp)"
      sh -c "command time -v true" 2> "$TMP_TEST" >/dev/null || true
      if grep -qi "Maximum resident set size" "$TMP_TEST" 2>/dev/null; then
        TIME_CMD="$(command -v time)"
      fi
      rm -f "$TMP_TEST"
    fi
  fi

  if [[ -z "$TIME_CMD" ]]; then
    echo -e "${YELLOW}Warning: No suitable 'time -v' found. Proceeding without time-monitor.${NC}"
    TIME_MONITOR=false
  else
    echo "Using time command: $TIME_CMD"
  fi
fi

# Start GPU monitor
MON_PID=""
if [[ "$GPU_MONITOR" = true && -x "$(command -v nvidia-smi)" ]]; then
  (
    curmax=0
    while true; do
      gpu_val=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | awk 'BEGIN{m=0} {v=int($1); if(v>m) m=v} END{print m+0}')
      if [[ -n "$gpu_val" && "$gpu_val" =~ ^[0-9]+$ ]]; then
        if (( gpu_val > curmax )); then
          curmax=$gpu_val
        fi
      fi
      if [[ -f "$TEMP_DIR/.baseline_done" ]]; then
        break
      fi
      sleep 0.2
    done
    echo "$curmax" > "$MON_TMP"
  ) &
  MON_PID=$!
else
  if [[ "$GPU_MONITOR" = true ]]; then
    echo -e "${YELLOW}Warning: nvidia-smi not found. Skipping GPU monitor.${NC}"
    GPU_MONITOR=false
  fi
fi

run_with_time() {
  if [[ "${TIME_MONITOR:-false}" = true && -n "$TIME_CMD" ]]; then
    if [[ "$TIME_CMD" == "/usr/bin/time" ]]; then
      "$TIME_CMD" -v -o "$TIME_TMP" "$@"
    else
      "$TIME_CMD" -v "$@" 2> "$TIME_TMP"
    fi
  else
    "$@"
  fi
}

BASELINE_EXIT_CODE=0

case "$METHOD" in
  aster)
    ASTER_ROOT="${BASELINES_DIR%/}/ASTER"
    if [[ -n "$ASTER_BIN" ]]; then
      if [[ "$ASTER_BIN" = /* ]]; then
        ASTER_BIN_PATH="$ASTER_BIN"
      else
        ASTER_BIN_PATH="${ASTER_ROOT%/}/$ASTER_BIN"
      fi
    else
      ASTER_BIN_PATH="${ASTER_ROOT%/}/bin/astral4"
    fi

    if [[ ! -x "$ASTER_BIN_PATH" ]]; then
      echo -e "${RED}Error: ASTER binary not found or not executable at $ASTER_BIN_PATH${NC}" >&2
      exit 5
    fi

    ASTER_OPTS_ARR=()
    if [[ -n "$ASTER_OPTS" ]]; then
      # shellcheck disable=SC2206
      ASTER_OPTS_ARR=( $ASTER_OPTS )
    fi

    echo -e "${YELLOW}Running ASTER...${NC}"
    echo -e "${YELLOW}Command: $ASTER_BIN_PATH -i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\" ${ASTER_OPTS}${NC}"
    run_with_time "$ASTER_BIN_PATH" -i "$INPUT_FILE" -o "$OUTPUT_FILE" "${ASTER_OPTS_ARR[@]}"
    BASELINE_EXIT_CODE=$?
    ;;
  astral)
    ASTRAL_ROOT="${BASELINES_DIR%/}/ASTRAL"
    if [[ ! -x "${ASTRAL_ROOT%/}/run_astral.sh" ]]; then
      echo -e "${RED}Error: ASTRAL run_astral.sh not found in $ASTRAL_ROOT${NC}" >&2
      exit 6
    fi

    ASTRAL_OPTS_ARR=()
    if [[ -n "$ASTRAL_OPTS" ]]; then
      # shellcheck disable=SC2206
      ASTRAL_OPTS_ARR=( $ASTRAL_OPTS )
    fi

    echo -e "${YELLOW}Running ASTRAL...${NC}"
    ASTRAL_MEM_ARGS=()
    if [[ -n "$ASTRAL_XMS" ]]; then
      ASTRAL_MEM_ARGS+=(--xms "$ASTRAL_XMS")
    fi
    if [[ -n "$ASTRAL_XMX" ]]; then
      ASTRAL_MEM_ARGS+=(--xmx "$ASTRAL_XMX")
    fi

    echo -e "${YELLOW}Command: cd $ASTRAL_ROOT && ./run_astral.sh -i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\" ${ASTRAL_OPTS}${NC}"
    ( cd "$ASTRAL_ROOT" && run_with_time ./run_astral.sh -i "$INPUT_FILE" -o "$OUTPUT_FILE" "${ASTRAL_MEM_ARGS[@]}" "${ASTRAL_OPTS_ARR[@]}" )
    BASELINE_EXIT_CODE=$?
    ;;
  treeqmc)
    TREEQMC_ROOT="${BASELINES_DIR%/}/TREE-QMC"
    if [[ ! -x "${TREEQMC_ROOT%/}/tree-qmc" ]]; then
      TREEQMC_ROOT="${BASELINES_DIR%/}/tree-qmc"
    fi
    TREEQMC_BIN="${TREEQMC_ROOT%/}/tree-qmc"
    if [[ ! -x "$TREEQMC_BIN" ]]; then
      echo -e "${RED}Error: tree-qmc binary not found in $TREEQMC_ROOT${NC}" >&2
      exit 7
    fi

    TREEQMC_OPTS_ARR=()
    if [[ -n "$TREEQMC_OPTS" ]]; then
      # shellcheck disable=SC2206
      TREEQMC_OPTS_ARR=( $TREEQMC_OPTS )
    fi

    echo -e "${YELLOW}Running TreeQMC...${NC}"
    echo -e "${YELLOW}Command: $TREEQMC_BIN -i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\" ${TREEQMC_OPTS}${NC}"
    run_with_time "$TREEQMC_BIN" -i "$INPUT_FILE" -o "$OUTPUT_FILE" "${TREEQMC_OPTS_ARR[@]}"
    BASELINE_EXIT_CODE=$?
    ;;
  wqfmtree)
    WQFM_ROOT="${BASELINES_DIR%/}/wQFM-TREE"
    if [[ ! -x "${WQFM_ROOT%/}/run.sh" ]]; then
      echo -e "${RED}Error: wQFM-TREE run.sh not found in $WQFM_ROOT${NC}" >&2
      exit 8
    fi

    WQFM_OPTS_ARR=()
    if [[ -n "$WQFM_OPTS" ]]; then
      # shellcheck disable=SC2206
      WQFM_OPTS_ARR=( $WQFM_OPTS )
    fi

    echo -e "${YELLOW}Running wQFM-TREE...${NC}"
    echo -e "${YELLOW}Command: cd $WQFM_ROOT && ./run.sh \"$INPUT_FILE\" \"$OUTPUT_FILE\" ${WQFM_OPTS}${NC}"
    ( cd "$WQFM_ROOT" && run_with_time ./run.sh "$INPUT_FILE" "$OUTPUT_FILE" "${WQFM_OPTS_ARR[@]}" )
    BASELINE_EXIT_CODE=$?
    ;;
  supertriplets)
    SUPERTRIPLETS_ROOT="${BASELINES_DIR%/}/SuperTriplets"
    if [[ ! -x "${STELAR_ROOT%/}/run_supertriplets.sh" ]]; then
      echo -e "${RED}Error: run_supertriplets.sh not found in $STELAR_ROOT${NC}" >&2
      exit 9
    fi

    SUPERTRIPLETS_OPTS_ARR=()
    if [[ -n "$SUPERTRIPLETS_OPTS" ]]; then
      # shellcheck disable=SC2206
      SUPERTRIPLETS_OPTS_ARR=( $SUPERTRIPLETS_OPTS )
    fi

    echo -e "${YELLOW}Running SuperTriplets...${NC}"
    echo -e "${YELLOW}Command: $STELAR_ROOT/run_supertriplets.sh --script-dir $SUPERTRIPLETS_ROOT -i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\" ${SUPERTRIPLETS_OPTS}${NC}"
    run_with_time "$STELAR_ROOT/run_supertriplets.sh" --script-dir "$SUPERTRIPLETS_ROOT" -i "$INPUT_FILE" -o "$OUTPUT_FILE" "${SUPERTRIPLETS_OPTS_ARR[@]}"
    BASELINE_EXIT_CODE=$?
    ;;
  tmc)
    TMC_ROOT="${BASELINES_DIR%/}/TMC"
    if [[ ! -x "${STELAR_ROOT%/}/run_tmc.sh" ]]; then
      echo -e "${RED}Error: run_tmc.sh not found in $STELAR_ROOT${NC}" >&2
      exit 10
    fi

    TMC_OPTS_ARR=()
    if [[ -n "$TMC_OPTS" ]]; then
      # shellcheck disable=SC2206
      TMC_OPTS_ARR=( $TMC_OPTS )
    fi

    echo -e "${YELLOW}Running TMC...${NC}"
    echo -e "${YELLOW}Command: $STELAR_ROOT/run_tmc.sh --script-dir $TMC_ROOT -i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\" ${TMC_OPTS}${NC}"
    run_with_time "$STELAR_ROOT/run_tmc.sh" --script-dir "$TMC_ROOT" -i "$INPUT_FILE" -o "$OUTPUT_FILE" "${TMC_OPTS_ARR[@]}"
    BASELINE_EXIT_CODE=$?
    ;;
  *)
    echo -e "${RED}Error: Unsupported method '$METHOD'${NC}" >&2
    exit 11
    ;;
esac

# Stop GPU monitor
touch "$TEMP_DIR/.baseline_done"

END_NS=$(date +%s%N)
ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

if [[ -n "${MON_PID:-}" ]]; then
  wait "$MON_PID" 2>/dev/null || true
fi

MAX_GPU_VAL="NA"
if [[ -f "$MON_TMP" ]]; then
  MAX_GPU_VAL=$(cat "$MON_TMP" 2>/dev/null || echo "NA")
fi

if [[ "$MAX_GPU_VAL" =~ ^[0-9]+$ ]]; then
  MAX_GPU_MB=$(awk "BEGIN{printf \"%.3f\", ${MAX_GPU_VAL} * 1.024}")
else
  MAX_GPU_MB="NA"
fi

MAX_CPU_MB="NA"
if [[ -f "$TIME_TMP" && -s "$TIME_TMP" ]]; then
  if grep -qi "Maximum resident set size" "$TIME_TMP" 2>/dev/null; then
    MAX_RSS_KB=$(grep -i "Maximum resident set size" "$TIME_TMP" | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | awk '{print int($1)}' | head -n1)
  elif grep -qi "Maximum resident set size (kbytes)" "$TIME_TMP" 2>/dev/null; then
    MAX_RSS_KB=$(grep -i "Maximum resident set size (kbytes)" "$TIME_TMP" | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | awk '{print int($1)}' | head -n1)
  else
    MAX_RSS_KB=""
  fi

  if [[ -n "${MAX_RSS_KB:-}" && "${MAX_RSS_KB}" =~ ^[0-9]+$ ]]; then
    MAX_CPU_MB=$(awk "BEGIN{printf \"%.3f\", ${MAX_RSS_KB}/1024}")
  fi
fi

# Summary
echo
echo -e "${GREEN}=== Baseline Execution Summary ===${NC}"
echo "Method:           $METHOD"
echo "Status:           $(if [[ $BASELINE_EXIT_CODE -eq 0 ]]; then echo -e "${GREEN}SUCCESS${NC}"; else echo -e "${RED}FAILED (exit $BASELINE_EXIT_CODE)${NC}"; fi)"
echo "Running time:     ${RUNNING_TIME}s"
echo "Max CPU RAM:      ${MAX_CPU_MB} MB"
echo "Max GPU VRAM:     ${MAX_GPU_MB} MB"
echo "Input file:       $INPUT_FILE"
echo "Output file:      $OUTPUT_FILE"
echo "Output exists:    $(if [[ -f "$OUTPUT_FILE" ]]; then echo "Yes"; else echo "No"; fi)"

# Save stats
STATS_FILE="${OUTPUT_FILE%.tre}_stats.csv"
echo "algorithm,input_file,output_file,running_time_s,max_cpu_mb,max_gpu_mb,exit_code" > "$STATS_FILE"
echo "${METHOD},$(basename "$INPUT_FILE"),$(basename "$OUTPUT_FILE"),${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${BASELINE_EXIT_CODE}" >> "$STATS_FILE"
echo "Stats saved to:   $STATS_FILE"

# Notification
if [[ "$NO_NOTIFY" != true && -n "$NTFY_CHANNEL_NAME" ]]; then
  STATUS_EMOJI=$(if [[ $BASELINE_EXIT_CODE -eq 0 ]]; then echo "✅"; else echo "❌"; fi)
  STATUS_TEXT=$(if [[ $BASELINE_EXIT_CODE -eq 0 ]]; then echo "SUCCESS"; else echo "FAILED"; fi)

  NTFY_MSG="${STATUS_EMOJI} ${METHOD} ${STATUS_TEXT}

Input: $(basename "$INPUT_FILE")
Output: $(basename "$OUTPUT_FILE")

Time: ${RUNNING_TIME}s
CPU RAM: ${MAX_CPU_MB} MB
GPU VRAM: ${MAX_GPU_MB} MB"

  curl -s -d "$NTFY_MSG" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
  echo -e "${GREEN}Notification sent to ntfy.sh/${NTFY_CHANNEL_NAME}${NC}"
fi

if [[ $BASELINE_EXIT_CODE -eq 0 ]]; then
  echo -e "${GREEN}Program completed successfully!${NC}"
else
  echo -e "${RED}Program failed with exit code $BASELINE_EXIT_CODE${NC}"
fi

exit $BASELINE_EXIT_CODE
