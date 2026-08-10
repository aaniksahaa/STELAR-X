#!/usr/bin/env bash
# run-with-monitor.sh
# Simple wrapper for STELAR-X that monitors time, memory, and GPU usage
# Usage examples:
#   ./run-with-monitor.sh --input input.tre --output output.tre
#   ./run-with-monitor.sh --input input.tre --output output.tre --no-time-monitor --no-gpu-monitor
#
set -euo pipefail

# Defaults
INPUT_FILE=""
OUTPUT_FILE=""
STELAR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"  # resolve relative to script location
STELAR_OPTS=""        # optional extra args passed through to run.sh

# Monitoring options (DEFAULT: ON)
TIME_MONITOR=true     # when true: run `time -v` if available and capture stderr
GPU_MONITOR=true      # when true: sample nvidia-smi while stelar runs
NO_NOTIFY=false       # when true: skip ntfy.sh notifications
DEBUG=0               # set DEBUG=1 to enable set -x
STELAR_ARGS=()

print_help() {
  cat <<EOF
run-with-monitor.sh - STELAR-X wrapper with performance monitoring

Usage: $0 --input <input_file> --output <output_file> [options]

Required:
  --input, -i        Path to gene trees file
  --output, -o       Path to output species tree file

Optional:
  --stelar-root      Path to STELAR-X root directory (default: current directory)
  --stelar-opts      Extra STELAR options passed to run.sh (default: empty)
  -T, --threads N    Max CPU threads/cores to expose to Java
  --expansion, -e    Enable mixed bipartitions (default: OFF)
  --no-time-monitor  Disable time-monitoring (overrides default ON)
  --no-gpu-monitor   Disable GPU-monitoring (overrides default ON)
  --no-notify, -nn   Disable ntfy.sh notifications
  --debug            Enable shell tracing
  --help, -h         Show this message

Examples:
  $0 --input my_genes.tre --output my_output.tre
  $0 --input input.tre --output output.tre --cpu-parallel --expansion
  $0 --input input.tre --output output.tre --no-time-monitor --no-notify
EOF
}

# Parse arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT_FILE="$2"; shift 2 ;;
    -o|--output) OUTPUT_FILE="$2"; shift 2 ;;
    --stelar-root) STELAR_ROOT="$2"; shift 2 ;;
    --stelar-opts)
      read -r -a TMP_OPTS <<< "$2"
      STELAR_ARGS+=("${TMP_OPTS[@]}")
      shift 2
      ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
    --help|-h) print_help; exit 0 ;;
    --cpu|--cpu-parallel|--gpu|--gpu-parallel|--expansion|-e) STELAR_ARGS+=("$1"); shift ;;
    -T|--threads|-m|--mode|-s|--support|--branch-support|--lambda)
      STELAR_ARGS+=("$1" "$2"); shift 2 ;;
    -v|--verbose) STELAR_ARGS+=("$1"); shift ;;
    --) shift; STELAR_ARGS+=("$@"); break ;;
    *) POSITIONAL+=("$1"); shift ;;
  esac
done

# Positional args are no longer supported
if [[ ${#POSITIONAL[@]} -gt 0 ]]; then
  echo "Error: Positional arguments are not supported. Use -i/--input and -o/--output."
  print_help
  exit 1
fi

if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
  echo "Error: Both --input and --output are required."
  print_help
  exit 1
fi

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Convert to absolute paths
INPUT_FILE=$(realpath "$INPUT_FILE")
OUTPUT_FILE=$(realpath "$OUTPUT_FILE")
STELAR_ROOT=$(realpath "$STELAR_ROOT")

# Debug/tracing
if [[ "${DEBUG:-0}" = "1" ]]; then
  set -x
fi

echo "=== STELAR-X Monitor Wrapper ==="
echo "Input file:     $INPUT_FILE"
echo "Output file:    $OUTPUT_FILE"
echo "STELAR root:    $STELAR_ROOT"
echo "STELAR options: ${STELAR_ARGS[*]:-(defaults)}"
echo "Cross-tree recombination: enabled by --expansion"
echo "Time monitor:   $TIME_MONITOR"
echo "GPU monitor:    $GPU_MONITOR"
echo "Notifications:  $(if [[ "$NO_NOTIFY" = true ]]; then echo "disabled"; else echo "enabled"; fi)"
echo

STELAR_RUN_MODE="AUTO"
STELAR_THREADS="AUTO"
STELAR_OPTIONS_STR="${STELAR_ARGS[*]:-}"
for ((i=0; i<${#STELAR_ARGS[@]}; i++)); do
  case "${STELAR_ARGS[$i]}" in
    --cpu) STELAR_RUN_MODE="CPU_SINGLE" ;;
    --cpu-parallel) STELAR_RUN_MODE="CPU_PARALLEL" ;;
    --gpu|--gpu-parallel) STELAR_RUN_MODE="GPU_PARALLEL" ;;
    -m|--mode)
      if [[ -n "${STELAR_ARGS[$((i+1))]:-}" ]]; then
        STELAR_RUN_MODE="${STELAR_ARGS[$((i+1))]}"
      fi
      ;;
    -T|--threads)
      if [[ -n "${STELAR_ARGS[$((i+1))]:-}" ]]; then
        STELAR_THREADS="${STELAR_ARGS[$((i+1))]}"
      fi
      ;;
  esac
done
STELAR_STATS_TAG="$(printf "%s" "${STELAR_RUN_MODE,,}_t${STELAR_THREADS}" | tr -c 'a-z0-9_.-' '_')"

# Validate input file
if [[ ! -f "$INPUT_FILE" ]]; then
  echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}"
  exit 1
fi

# Validate STELAR root
if [[ ! -d "$STELAR_ROOT" ]]; then
  echo -e "${RED}Error: STELAR root directory '$STELAR_ROOT' does not exist.${NC}"
  exit 1
fi

if [[ ! -f "$STELAR_ROOT/run.sh" ]]; then
  echo -e "${RED}Error: run.sh not found in '$STELAR_ROOT'.${NC}"
  exit 1
fi

VERSION="$("$STELAR_ROOT/project-version.sh")"
JAR_PATH="$STELAR_ROOT/target/stelar-x-${VERSION}.jar"
if [[ ! -f "$JAR_PATH" ]]; then
  echo -e "${RED}Error: JAR not found. Please run ./install.sh first.${NC}"
  exit 1
fi

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Create temporary directory for logs
TEMP_DIR=$(mktemp -d)
TIME_TMP="$TEMP_DIR/stelar_time_err.log"
MON_TMP="$TEMP_DIR/stelar_gpu_mem.log"

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
echo "Temp directory: $TEMP_DIR"
echo "STELAR run.sh:  $STELAR_ROOT/run.sh"
echo "Library path:   $STELAR_ROOT/cuda"
echo "Library exists: $(if [[ -f "$STELAR_ROOT/cuda/libweight_calc.so" ]]; then echo "Yes"; else echo "No"; fi)"
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
      if [[ -f "$TEMP_DIR/.stelar_done" ]]; then
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

STELAR_PID=""
STELAR_EXIT_CODE=""
echo -e "${YELLOW}Running STELAR-X...${NC}"
echo -e "${YELLOW}Command: cd $STELAR_ROOT && ./run.sh --input \"$INPUT_FILE\" --output \"$OUTPUT_FILE\" ${STELAR_ARGS[*]}${NC}"
echo

# Capture both stdout and stderr to TIME_TMP (using tee to also show live output)
if [[ "${TIME_MONITOR:-false}" = true && -n "$TIME_CMD" ]]; then
  (
    cd "$STELAR_ROOT" && "$TIME_CMD" -v ./run.sh --input "$INPUT_FILE" --output "$OUTPUT_FILE" "${STELAR_ARGS[@]}" < /dev/null 2>&1 | tee "$TIME_TMP"
  ) &
  STELAR_PID=$!
else
  (
    cd "$STELAR_ROOT" && ./run.sh --input "$INPUT_FILE" --output "$OUTPUT_FILE" "${STELAR_ARGS[@]}" < /dev/null 2>&1 | tee "$TIME_TMP"
  ) &
  STELAR_PID=$!
fi

# Give a short moment and verify the job didn't die immediately
sleep 0.25
if ! kill -0 "$STELAR_PID" >/dev/null 2>&1; then
  touch "$TEMP_DIR/.stelar_done"
  if wait "$STELAR_PID"; then
    STELAR_EXIT_CODE=0
    echo "STELAR process (pid ${STELAR_PID}) finished before startup check completed."
  else
    STELAR_EXIT_CODE=$?
    echo -e "${RED}Error: STELAR process (pid ${STELAR_PID}) failed to start or died immediately.${NC}"
    echo "----- Captured time/generic stderr (first 200 lines) -----"
    head -n 200 "$TIME_TMP" 2>/dev/null || true
    echo "---------------------------------------------------------"
    exit 5
  fi
fi

if [[ -z "$STELAR_EXIT_CODE" ]]; then
  echo "STELAR started with PID ${STELAR_PID} (logging to ${TIME_TMP} if time-monitor enabled)"

  # Wait for stelar to finish
  if wait "$STELAR_PID"; then
    STELAR_EXIT_CODE=0
  else
    STELAR_EXIT_CODE=$?
  fi
fi

# Stop GPU monitor by placing sentinel file
touch "$TEMP_DIR/.stelar_done"

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

# Extract the optimal triplet score from STELAR output
OPTIMAL_TRIPLET_SCORE="NA"
NUM_TAXA="NA"
NUM_GENE_TREES="NA"
NORMALIZED_SCORE="NA"

if [[ -f "$TIME_TMP" ]]; then
  SCORE_LINE=$(grep "OPTIMAL_TRIPLET_SCORE:" "$TIME_TMP" 2>/dev/null | tail -n1 || true)
  if [[ -n "$SCORE_LINE" ]]; then
    OPTIMAL_TRIPLET_SCORE=$(echo "$SCORE_LINE" | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | tr -d ' ')
  fi

  # Extract taxa count and gene tree count for normalization
  # Look for "Real taxa count (dataset): 37"
  TAXA_LINE=$(grep "Real taxa count (dataset):" "$TIME_TMP" 2>/dev/null | tail -n1 || true)
  if [[ -n "$TAXA_LINE" ]]; then
     NUM_TAXA=$(echo "$TAXA_LINE" | awk -F: '{print $2}' | tr -d ' ')
  fi
  
  # Look for "Trees: 200" or similar. 
  # In the log: "Range candidates: 1085\nGene tree bipartitions: 1085\nTrees: 200"
  # Also "Inverse Index Manager Statistics:\n  Trees: 200"
  TREES_LINE=$(grep "Trees:" "$TIME_TMP" 2>/dev/null | head -n1 || true)
  if [[ -n "$TREES_LINE" ]]; then
     NUM_GENE_TREES=$(echo "$TREES_LINE" | awk -F: '{print $2}' | tr -d ' ')
  fi
fi

# Calculate normalized score if we have all components
if [[ "$OPTIMAL_TRIPLET_SCORE" != "NA" && "$NUM_TAXA" != "NA" && "$NUM_GENE_TREES" != "NA" ]]; then
  # Remove decimals for taxa/trees just in case
  n=$(echo "$NUM_TAXA" | awk '{print int($1)}')
  k=$(echo "$NUM_GENE_TREES" | awk '{print int($1)}')
  
  if [[ "$n" -ge 3 && "$k" -gt 0 ]]; then
      # nC3 = n * (n-1) * (n-2) / 6
      TOTAL_TRIPLETS=$(echo "$n * ($n - 1) * ($n - 2) / 6" | bc)
      
      # Normalized = Score / (k * nC3)
      # Using awk for floating point division
      NORMALIZED_SCORE=$(awk -v score="$OPTIMAL_TRIPLET_SCORE" -v k="$k" -v trips="$TOTAL_TRIPLETS" 'BEGIN { printf "%.6f", score / (k * trips) }')
  fi
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

# If time-monitor wasn't used or cpu RSS not parsed, try a best-effort read of memory using ps for the stelar PID (peak not available):
if [[ "$TIME_MONITOR" = false && "$MAX_CPU_MB" == "NA" ]]; then
  if ps -p "$STELAR_PID" >/dev/null 2>&1; then
    MAX_CPU_MB=$(ps -o rss= -p "$STELAR_PID" 2>/dev/null | awk '{print ($1+0)/1024}')
  else
    MAX_CPU_MB="NA"
  fi
fi

echo
echo -e "${GREEN}=== STELAR-X Execution Summary ===${NC}"
echo "Status:         $(if [[ $STELAR_EXIT_CODE -eq 0 ]]; then echo -e "${GREEN}SUCCESS${NC}"; else echo -e "${RED}FAILED (exit code $STELAR_EXIT_CODE)${NC}"; fi)"
echo "Running time:   ${RUNNING_TIME}s"
echo "Max CPU RAM:    ${MAX_CPU_MB} MB"
echo "Max GPU VRAM:   ${MAX_GPU_MB} MB"
echo
echo "Optimal Triplet Score:    ${OPTIMAL_TRIPLET_SCORE}"
echo "Normalized Triplet Score: ${NORMALIZED_SCORE}"
echo
echo "Input file:     $INPUT_FILE"
echo "Output file:    $OUTPUT_FILE"
echo "Output exists:  $(if [[ -f "$OUTPUT_FILE" ]]; then echo "Yes"; else echo "No"; fi)"
if [[ -f "$OUTPUT_FILE" ]]; then
  echo "Output size:    $(wc -l < "$OUTPUT_FILE") lines"
fi

# Create stats files next to the output. Keep the legacy path for existing scripts,
# and also write a mode-specific file so CPU/GPU runs do not overwrite each other.
STATS_FILE="${OUTPUT_FILE%.tre}_stats.csv"
MODE_STATS_FILE="${OUTPUT_FILE%.tre}_${STELAR_STATS_TAG}_stats.csv"
STATS_HEADER="algorithm,mode,threads,stelar_options,input_file,output_file,running_time_s,max_cpu_mb,max_gpu_mb,optimal_triplet_score,normalized_triplet_score,exit_code"
STATS_ROW="stelar-x,${STELAR_RUN_MODE},${STELAR_THREADS},\"${STELAR_OPTIONS_STR}\",$(basename "$INPUT_FILE"),$(basename "$OUTPUT_FILE"),${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${OPTIMAL_TRIPLET_SCORE},${NORMALIZED_SCORE},${STELAR_EXIT_CODE}"
echo "$STATS_HEADER" > "$STATS_FILE"
echo "$STATS_ROW" >> "$STATS_FILE"
echo "$STATS_HEADER" > "$MODE_STATS_FILE"
echo "$STATS_ROW" >> "$MODE_STATS_FILE"
echo "Stats saved to: $STATS_FILE"
echo "Mode stats saved to: $MODE_STATS_FILE"

# Send notification (ntfy) if enabled and curl available
if [[ "$NO_NOTIFY" = false ]] && command -v curl >/dev/null 2>&1; then
  STATUS_EMOJI=$(if [[ $STELAR_EXIT_CODE -eq 0 ]]; then echo "🎉"; else echo "❌"; fi)
  STATUS_TEXT=$(if [[ $STELAR_EXIT_CODE -eq 0 ]]; then echo "completed successfully"; else echo "failed (exit $STELAR_EXIT_CODE)"; fi)
  
  curl -s -d "${STATUS_EMOJI} STELAR-X ${STATUS_TEXT}!

📊 Performance:
• Running time: ${RUNNING_TIME}s
• Max CPU RAM: ${MAX_CPU_MB} MB
• Max GPU VRAM: ${MAX_GPU_MB} MB
• Normalized Triplet Score: ${NORMALIZED_SCORE}
• Optimal Triplet Score: ${OPTIMAL_TRIPLET_SCORE}

📁 Files:
• Input: $(basename "$INPUT_FILE")
• Output: $(basename "$OUTPUT_FILE")
• Stats: $(basename "$STATS_FILE")" ntfy.sh/anik-test || true
fi

echo
if [[ $STELAR_EXIT_CODE -eq 0 ]]; then
  echo -e "${GREEN}Program completed successfully!${NC}"
else
  echo -e "${RED}Program failed with exit code $STELAR_EXIT_CODE${NC}"
  if [[ -f "$TIME_TMP" ]]; then
    echo
    echo "=== Error Details (from time/stderr log) ==="
    tail -n 50 "$TIME_TMP" 2>/dev/null || true
    echo "============================================="
  fi
fi

exit $STELAR_EXIT_CODE
