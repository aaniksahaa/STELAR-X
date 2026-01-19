#!/usr/bin/env bash
# run-aster-with-monitor.sh
# ASTER wrapper that monitors time, memory, and GPU usage
# All arguments except monitoring flags are passed directly to ASTER
#
# For detailed help: ./run-aster-with-monitor.sh -h

set -euo pipefail

# ===== CONFIGURATION =====
NTFY_CHANNEL_NAME="anik-test"
ASTER_BIN="ASTER/bin/astral4"

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
run-aster-with-monitor.sh - ASTER wrapper with performance monitoring

Usage:
    $0 [monitoring-options] -- [aster-args]
    $0 -i <input> -o <output> [aster-options]
    $0 --help

Monitoring Options (parsed by this script):
    --aster-bin PATH        Path to ASTER binary (default: ASTER/bin/astral4)
    --no-time-monitor       Disable time-monitoring
    --no-gpu-monitor        Disable GPU-monitoring
    --no-notify, -nn        Disable ntfy.sh notifications
    --debug                 Enable shell tracing

All Other Arguments:
    Passed directly to ASTER. Common ASTER options include:
    -i <input>              Input gene trees file
    -o <output>             Output species tree file
    -t <threads>            Number of threads

Examples:
    # Basic usage (all args passed to ASTER)
    $0 -i genes.tre -o species.tre -t 8

    # Disable GPU monitoring
    $0 --no-gpu-monitor -i genes.tre -o species.tre

    # Use custom ASTER binary
    $0 --aster-bin /path/to/astral -i genes.tre -o out.tre

    # Disable notifications
    $0 --no-notify -i genes.tre -o species.tre -t 4
EOF
}

# Collect ASTER arguments
ASTER_ARGS=()

# Parse arguments - only extract monitoring flags, pass rest to ASTER
while [[ $# -gt 0 ]]; do
  case "$1" in
    --aster-bin)
      ASTER_BIN="$2"
      shift 2
      ;;
    --no-time-monitor)
      TIME_MONITOR=false
      shift
      ;;
    --no-gpu-monitor)
      GPU_MONITOR=false
      shift
      ;;
    --no-notify|-nn)
      NO_NOTIFY=true
      shift
      ;;
    --debug)
      DEBUG=1
      shift
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    --)
      # Everything after -- goes to ASTER
      shift
      ASTER_ARGS+=("$@")
      break
      ;;
    *)
      # Everything else goes to ASTER
      ASTER_ARGS+=("$1")
      shift
      ;;
  esac
done

# Debug/tracing
if [[ "${DEBUG:-0}" = "1" ]]; then
  set -x
fi

# Extract input/output from ASTER args for display
INPUT_FILE=""
OUTPUT_FILE=""
for i in "${!ASTER_ARGS[@]}"; do
  if [[ "${ASTER_ARGS[$i]}" == "-i" ]] && [[ $((i+1)) -lt ${#ASTER_ARGS[@]} ]]; then
    INPUT_FILE="${ASTER_ARGS[$((i+1))]}"
  fi
  if [[ "${ASTER_ARGS[$i]}" == "-o" ]] && [[ $((i+1)) -lt ${#ASTER_ARGS[@]} ]]; then
    OUTPUT_FILE="${ASTER_ARGS[$((i+1))]}"
  fi
done

echo "=== ASTER Monitor Wrapper ==="
echo "ASTER binary:     $ASTER_BIN"
echo "ASTER args:       ${ASTER_ARGS[*]:-<none>}"
echo "Input file:       ${INPUT_FILE:-<from args>}"
echo "Output file:      ${OUTPUT_FILE:-<from args>}"
echo "Time monitor:     $TIME_MONITOR"
echo "GPU monitor:      $GPU_MONITOR"
echo "Notifications:    $(if [[ "$NO_NOTIFY" = true ]]; then echo "disabled"; else echo "enabled ($NTFY_CHANNEL_NAME)"; fi)"
echo

# Validate ASTER binary
if [[ ! -x "$ASTER_BIN" ]]; then
  echo -e "${RED}Error: ASTER binary '$ASTER_BIN' not found or not executable.${NC}" >&2
  exit 1
fi

# Create temporary directory for logs
TEMP_DIR=$(mktemp -d)
TIME_TMP="$TEMP_DIR/aster_time_err.log"
MON_TMP="$TEMP_DIR/aster_gpu_mem.log"

# Cleanup function
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
      if [[ -f "$TEMP_DIR/.aster_done" ]]; then
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

# Launch ASTER
ASTER_PID=""
echo -e "${YELLOW}Running ASTER...${NC}"
echo -e "${YELLOW}Command: $ASTER_BIN ${ASTER_ARGS[*]}${NC}"
echo

if [[ "${TIME_MONITOR:-false}" = true && -n "$TIME_CMD" ]]; then
  (
    eval "$TIME_CMD -v $ASTER_BIN ${ASTER_ARGS[*]}" < /dev/null
  ) 2> "$TIME_TMP" &
  ASTER_PID=$!
else
  (
    eval "$ASTER_BIN ${ASTER_ARGS[*]}" < /dev/null
  ) &
  ASTER_PID=$!
fi

# Verify ASTER started
sleep 0.25
if ! kill -0 "$ASTER_PID" >/dev/null 2>&1; then
  echo -e "${RED}Error: ASTER process failed to start.${NC}"
  head -n 50 "$TIME_TMP" 2>/dev/null || true
  touch "$TEMP_DIR/.aster_done"
  exit 5
fi

echo "ASTER started with PID ${ASTER_PID}"

# Wait for ASTER
wait "$ASTER_PID"
ASTER_EXIT_CODE=$?

# Stop GPU monitor
touch "$TEMP_DIR/.aster_done"

END_NS=$(date +%s%N)
ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

# Get GPU memory
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

# Get CPU memory
MAX_CPU_MB="NA"
if [[ -f "$TIME_TMP" && -s "$TIME_TMP" ]]; then
  if grep -qi "Maximum resident set size" "$TIME_TMP" 2>/dev/null; then
    MAX_RSS_KB=$(grep -i "Maximum resident set size" "$TIME_TMP" | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | awk '{print int($1)}' | head -n1)
    if [[ -n "${MAX_RSS_KB:-}" && "${MAX_RSS_KB}" =~ ^[0-9]+$ ]]; then
      MAX_CPU_MB=$(awk "BEGIN{printf \"%.3f\", ${MAX_RSS_KB}/1024}")
    fi
  fi
fi

# Print summary
echo
echo -e "${GREEN}=== ASTER Execution Summary ===${NC}"
echo "Status:           $(if [[ $ASTER_EXIT_CODE -eq 0 ]]; then echo -e "${GREEN}SUCCESS${NC}"; else echo -e "${RED}FAILED (exit $ASTER_EXIT_CODE)${NC}"; fi)"
echo "Running time:     ${RUNNING_TIME}s"
echo "Max CPU RAM:      ${MAX_CPU_MB} MB"
echo "Max GPU VRAM:     ${MAX_GPU_MB} MB"
echo "Input file:       ${INPUT_FILE:-<from args>}"
echo "Output file:      ${OUTPUT_FILE:-<from args>}"
if [[ -n "$OUTPUT_FILE" && -f "$OUTPUT_FILE" ]]; then
  echo "Output exists:    Yes"
  echo "Output size:      $(wc -l < "$OUTPUT_FILE") lines"
fi

# Save stats
if [[ -n "$OUTPUT_FILE" ]]; then
  STATS_FILE="${OUTPUT_FILE%.tre}_stats.csv"
  echo "algorithm,input_file,output_file,running_time_s,max_cpu_mb,max_gpu_mb,exit_code" > "$STATS_FILE"
  echo "aster,$(basename "${INPUT_FILE:-unknown}"),$(basename "${OUTPUT_FILE:-unknown}"),${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${ASTER_EXIT_CODE}" >> "$STATS_FILE"
  echo "Stats saved to:   $STATS_FILE"
fi

# Send ntfy notification
if [[ "$NO_NOTIFY" != true && -n "$NTFY_CHANNEL_NAME" ]]; then
  STATUS_EMOJI=$(if [[ $ASTER_EXIT_CODE -eq 0 ]]; then echo "✅"; else echo "❌"; fi)
  STATUS_TEXT=$(if [[ $ASTER_EXIT_CODE -eq 0 ]]; then echo "SUCCESS"; else echo "FAILED"; fi)
  
  NTFY_MSG="${STATUS_EMOJI} ASTER ${STATUS_TEXT}

Input: $(basename "${INPUT_FILE:-unknown}")
Output: $(basename "${OUTPUT_FILE:-unknown}")

Time: ${RUNNING_TIME}s
CPU RAM: ${MAX_CPU_MB} MB
GPU VRAM: ${MAX_GPU_MB} MB"

  curl -s -d "$NTFY_MSG" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
  echo -e "${GREEN}Notification sent to ntfy.sh/${NTFY_CHANNEL_NAME}${NC}"
fi

echo
if [[ $ASTER_EXIT_CODE -eq 0 ]]; then
  echo -e "${GREEN}ASTER completed successfully!${NC}"
else
  echo -e "${RED}ASTER failed with exit code $ASTER_EXIT_CODE${NC}"
fi

exit $ASTER_EXIT_CODE
