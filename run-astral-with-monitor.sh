#!/usr/bin/env bash
# run-astral-with-monitor.sh
# Simple wrapper for ASTRAL-MP that monitors time, memory, and GPU usage
# Usage examples:
#   ./run-astral-with-monitor.sh input.tre output.tre
#   ./run-astral-with-monitor.sh input.tre output.tre --no-time-monitor --no-gpu-monitor
#   ./run-astral-with-monitor.sh input.tre output.tre --astral-opts "-t 4"
#
set -euo pipefail

# Defaults
INPUT_FILE=""
OUTPUT_FILE=""
ASTRAL_ROOT=""               # must be provided or derived from BASE_DIR
ASTRAL_ROOT_SET=false
BASE_DIR="${HOME}/phylogeny"
BASE_DIR="${HOME}/phylogeny"
ASTRAL_OPTS=""               # default ASTRAL options

# Monitoring options (DEFAULT: ON)
TIME_MONITOR=true     # when true: run `time -v` if available and capture stderr
GPU_MONITOR=true      # when true: sample nvidia-smi while astral runs
NO_NOTIFY=false       # when true: skip ntfy.sh notifications
DEBUG=0               # set DEBUG=1 to enable set -x

print_help() {
  cat <<EOF
run-astral-with-monitor.sh - ASTRAL-MP wrapper with performance monitoring

Usage: $0 <input_file> <output_file> [options]

Required:
  input_file         Path to gene trees file
  output_file        Path to output species tree file

Optional:
  --astral-root      Path to ASTRAL-MP root directory (default: \$HOME/phylogeny/ASTRAL)
  --base-dir, -b     Base directory for deriving ASTRAL root (default: $BASE_DIR)
  --astral-opts      ASTRAL options (default: "$ASTRAL_OPTS")
  --no-time-monitor  Disable time-monitoring (overrides default ON)
  --no-gpu-monitor   Disable GPU-monitoring (overrides default ON)
  --no-notify, -nn   Disable ntfy.sh notifications
  --debug            Enable shell tracing
  --help, -h         Show this message

Examples:
  $0 my_genes.tre my_output.tre
  $0 input.tre output.tre --astral-opts "-t 4"
  $0 input.tre output.tre --base-dir /home/user/phylogeny
  $0 input.tre output.tre --astral-root /path/to/ASTRAL --no-notify
  $0 input.tre output.tre --no-time-monitor --no-gpu-monitor
EOF
}

# Parse arguments
if [[ $# -lt 2 ]]; then
  echo "Error: Input and output files are required."
  print_help
  exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"
shift 2

# Parse remaining options
while [[ $# -gt 0 ]]; do
  case "$1" in
    --astral-root) ASTRAL_ROOT="$2"; ASTRAL_ROOT_SET=true; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --astral-opts) ASTRAL_OPTS="$2"; shift 2 ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Derive ASTRAL_ROOT from BASE_DIR if not explicitly set
if [[ "$ASTRAL_ROOT_SET" = false ]]; then
  ASTRAL_ROOT="${BASE_DIR%/}/ASTRAL"
fi

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Convert to absolute paths
INPUT_FILE=$(realpath "$INPUT_FILE")
OUTPUT_FILE=$(realpath "$OUTPUT_FILE")
ASTRAL_ROOT=$(realpath "$ASTRAL_ROOT")

# Debug/tracing
if [[ "${DEBUG:-0}" = "1" ]]; then
  set -x
fi

echo "=== ASTRAL-MP Monitor Wrapper ==="
echo "Input file:     $INPUT_FILE"
echo "Output file:    $OUTPUT_FILE"
echo "ASTRAL root:    $ASTRAL_ROOT"
echo "ASTRAL options: $ASTRAL_OPTS"
echo "Time monitor:   $TIME_MONITOR"
echo "GPU monitor:    $GPU_MONITOR"
echo "Notifications:  $(if [[ "$NO_NOTIFY" = true ]]; then echo "disabled"; else echo "enabled"; fi)"
echo

# Validate input file
if [[ ! -f "$INPUT_FILE" ]]; then
  echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}"
  exit 1
fi

# Validate ASTRAL root
if [[ ! -d "$ASTRAL_ROOT" ]]; then
  echo -e "${RED}Error: ASTRAL root directory '$ASTRAL_ROOT' does not exist.${NC}"
  echo "Please specify the correct path with --astral-root or --base-dir"
  exit 1
fi

if [[ ! -f "$ASTRAL_ROOT/run_astral.sh" ]]; then
  echo -e "${RED}Error: run_astral.sh not found in '$ASTRAL_ROOT'.${NC}"
  echo "Please ensure ASTRAL-MP is properly installed in the specified directory."
  exit 1
fi

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Create temporary directory for logs
TEMP_DIR=$(mktemp -d)
TIME_TMP="$TEMP_DIR/astral_time_err.log"
MON_TMP="$TEMP_DIR/astral_gpu_mem.log"

# Cleanup function
cleanup() {
  rm -rf "$TEMP_DIR" 2>/dev/null || true
  if [[ -n "${MON_PID:-}" ]]; then
    kill "$MON_PID" 2>/dev/null || true
    wait "$MON_PID" 2>/dev/null || true
  fi
}
trap cleanup EXIT

echo -e "${YELLOW}Debug Information:${NC}"
echo "Temp directory:    $TEMP_DIR"
echo "ASTRAL run_astral.sh: $ASTRAL_ROOT/run_astral.sh"
echo

# START timer
START_NS=$(date +%s%N)

# --- Detect a usable 'time' command that supports -v ---
TIME_CMD=""
if [[ "${TIME_MONITOR:-false}" = true ]]; then
  # prefer /usr/bin/time if present and executable
  if [[ -x "/usr/bin/time" ]]; then
    TIME_CMD="/usr/bin/time"
  else
    # try to find a binary 'time' that supports -v
    if command -v time >/dev/null 2>&1; then
      # test it: run 'time -v true' via sh -c and capture stderr
      TMP_TEST="$(mktemp)"
      # Use command to prefer external time when available; shell builtin 'time' won't write 'Maximum resident set size'
      sh -c "command time -v true" 2> "$TMP_TEST" >/dev/null || true
      if grep -qi "Maximum resident set size" "$TMP_TEST" 2>/dev/null; then
        TIME_CMD="$(command -v time)"
      fi
      rm -f "$TMP_TEST"
    fi
  fi

  if [[ -z "$TIME_CMD" ]]; then
    echo -e "${YELLOW}Warning: time-monitor requested (default) but no suitable 'time -v' binary found.${NC}"
    echo "  - install it on Debian/Ubuntu with: sudo apt update && sudo apt install -y time"
    echo "Proceeding without time-monitor; GPU-monitor (if enabled) will still run."
    TIME_MONITOR=false
  else
    echo "Using time command: $TIME_CMD"
  fi
fi

# Start GPU monitor if requested and available
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
      if [[ -f "$TEMP_DIR/.astral_done" ]]; then
        break
      fi
      sleep 0.2
    done
    echo "$curmax" > "$MON_TMP"
  ) &
  MON_PID=$!
else
  if [[ "$GPU_MONITOR" = true ]]; then
    echo -e "${YELLOW}Warning: GPU monitor requested but nvidia-smi not found or not executable. Skipping GPU monitor.${NC}"
    GPU_MONITOR=false
  fi
fi

# Launch ASTRAL -- prefer using TIME_CMD if available, otherwise run directly
ASTRAL_PID=""
echo -e "${YELLOW}Running ASTRAL-MP...${NC}"
echo -e "${YELLOW}Command: cd $ASTRAL_ROOT && ./run_astral.sh -i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\" $ASTRAL_OPTS${NC}"
echo

if [[ "${TIME_MONITOR:-false}" = true && -n "$TIME_CMD" ]]; then
  # Use tee to show ASTRAL output while capturing time stats
  (
    cd "$ASTRAL_ROOT"
    "$TIME_CMD" -v ./run_astral.sh -i "$INPUT_FILE" -o "$OUTPUT_FILE" $ASTRAL_OPTS 2>&1 | tee >(grep -E "(Command being timed|User time|System time|Elapsed|Maximum resident set size|Exit status)" > "$TIME_TMP")
  ) &
  ASTRAL_PID=$!
else
  (
    cd "$ASTRAL_ROOT" && ./run_astral.sh -i "$INPUT_FILE" -o "$OUTPUT_FILE" $ASTRAL_OPTS
  ) &
  ASTRAL_PID=$!
fi

# Give a short moment and verify the job didn't die immediately
sleep 0.25
if ! kill -0 "$ASTRAL_PID" >/dev/null 2>&1; then
  echo -e "${RED}Error: ASTRAL process (pid ${ASTRAL_PID}) failed to start or died immediately.${NC}"
  echo "----- Captured time/generic stderr (first 200 lines) -----"
  head -n 200 "$TIME_TMP" 2>/dev/null || true
  echo "---------------------------------------------------------"
  touch "$TEMP_DIR/.astral_done"
  exit 5
fi

echo "ASTRAL started with PID ${ASTRAL_PID} (logging to ${TIME_TMP} if time-monitor enabled)"

# Wait for astral to finish
wait "$ASTRAL_PID"
ASTRAL_EXIT_CODE=$?

# Stop GPU monitor by placing sentinel file
touch "$TEMP_DIR/.astral_done"

END_NS=$(date +%s%N)
ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

# ensure monitor stopped and read its result
if [[ -n "${MON_PID:-}" ]]; then
  wait "$MON_PID" 2>/dev/null || true
fi

MAX_GPU_VAL="NA"
if [[ -f "$MON_TMP" ]]; then
  MAX_GPU_VAL=$(cat "$MON_TMP" 2>/dev/null || echo "NA")
fi

# Convert GPU MiB -> MB (decimal) if numeric
if [[ "$MAX_GPU_VAL" =~ ^[0-9]+$ ]]; then
  MAX_GPU_MB=$(awk "BEGIN{printf \"%.3f\", ${MAX_GPU_VAL} * 1.024}")
else
  MAX_GPU_MB="NA"
fi

# If time-monitor was used, try to parse its Maximum resident set size
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

# If time-monitor wasn't used or cpu RSS not parsed, try a best-effort read of memory using ps for the astral PID (peak not available):
if [[ "$TIME_MONITOR" = false && "$MAX_CPU_MB" == "NA" ]]; then
  if ps -p "$ASTRAL_PID" >/dev/null 2>&1; then
    MAX_CPU_MB=$(ps -o rss= -p "$ASTRAL_PID" 2>/dev/null | awk '{print ($1+0)/1024}')
  else
    MAX_CPU_MB="NA"
  fi
fi

echo
echo -e "${GREEN}=== ASTRAL-MP Execution Summary ===${NC}"
echo "Status:         $(if [[ $ASTRAL_EXIT_CODE -eq 0 ]]; then echo -e "${GREEN}SUCCESS${NC}"; else echo -e "${RED}FAILED (exit code $ASTRAL_EXIT_CODE)${NC}"; fi)"
echo "Running time:   ${RUNNING_TIME}s"
echo "Max CPU RAM:    ${MAX_CPU_MB} MB"
echo "Max GPU VRAM:   ${MAX_GPU_MB} MB"
echo "Input file:     $INPUT_FILE"
echo "Output file:    $OUTPUT_FILE"
echo "Output exists:  $(if [[ -f "$OUTPUT_FILE" ]]; then echo "Yes"; else echo "No"; fi)"
if [[ -f "$OUTPUT_FILE" ]]; then
  echo "Output size:    $(wc -l < "$OUTPUT_FILE") lines"
fi

# Create a simple stats file next to the output
STATS_FILE="${OUTPUT_FILE%.tre}_stats.csv"
echo "algorithm,input_file,output_file,running_time_s,max_cpu_mb,max_gpu_mb,exit_code" > "$STATS_FILE"
echo "astral-mp,$(basename "$INPUT_FILE"),$(basename "$OUTPUT_FILE"),${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${ASTRAL_EXIT_CODE}" >> "$STATS_FILE"
echo "Stats saved to: $STATS_FILE"

# Send notification (ntfy) if enabled and curl available
if [[ "$NO_NOTIFY" = false ]] && command -v curl >/dev/null 2>&1; then
  STATUS_EMOJI=$(if [[ $ASTRAL_EXIT_CODE -eq 0 ]]; then echo "🎉"; else echo "❌"; fi)
  STATUS_TEXT=$(if [[ $ASTRAL_EXIT_CODE -eq 0 ]]; then echo "completed successfully"; else echo "failed (exit $ASTRAL_EXIT_CODE)"; fi)
  
  curl -s -d "${STATUS_EMOJI} ASTRAL-MP ${STATUS_TEXT}!

📊 Performance:
• Running time: ${RUNNING_TIME}s
• Max CPU RAM: ${MAX_CPU_MB} MB
• Max GPU VRAM: ${MAX_GPU_MB} MB

📁 Files:
• Input: $(basename "$INPUT_FILE")
• Output: $(basename "$OUTPUT_FILE")
• Stats: $(basename "$STATS_FILE")" ntfy.sh/anik-test || true
fi

echo
if [[ $ASTRAL_EXIT_CODE -eq 0 ]]; then
  echo -e "${GREEN}Program completed successfully!${NC}"
else
  echo -e "${RED}Program failed with exit code $ASTRAL_EXIT_CODE${NC}"
  if [[ -f "$TIME_TMP" ]]; then
    echo
    echo "=== Error Details (from time/stderr log) ==="
    tail -n 50 "$TIME_TMP" 2>/dev/null || true
    echo "============================================="
  fi
fi

exit $ASTRAL_EXIT_CODE
