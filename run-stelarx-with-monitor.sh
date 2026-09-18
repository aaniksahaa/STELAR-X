#!/usr/bin/env bash
# run-stelarx-with-monitor.sh
# Wrapper for STELAR-X that records wall time, CPU RAM, GPU VRAM, and summary
# stats while preserving the existing research-script structure.

set -euo pipefail

# Propagate terminal color preference to Java subprocesses even when stderr is
# piped through tee.  Evaluated here, before any pipe redirection is applied to
# this script's file descriptors, so [[ -t ]] correctly reflects whether a real
# terminal is connected.  Java's Banner.detectColor() honours FORCE_COLOR.
[[ -t 1 || -t 2 ]] && export FORCE_COLOR=1

NTFY_CHANNEL_NAME="${NTFY_CHANNEL_NAME:-anik-phylo-stx}"

INPUT_FILE=""
OUTPUT_FILE=""
STELARX_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIME_MONITOR=true
GPU_MONITOR=true
NO_NOTIFY=false
DEBUG=0
STELARX_ARGS=()
REFERENCE_SPECIES_TREE=""

print_help() {
  cat <<EOF
run-stelarx-with-monitor.sh - STELAR-X wrapper with performance monitoring

Usage: $0 --input <input_file> --output <output_file> [options]

Required:
  --input, -i           Path to gene trees file
  --output, -o          Path to output species tree file

Optional:
  --stelarx-root        Path to STELAR-X root directory (default: current directory)
  --stelar-root         Compatibility alias for --stelarx-root
  --opts "..."          Extra algorithm options passed to run.sh
  --alg-opts "..."      Alias for --opts
  --stelarx-opts "..."  Compatibility alias for --opts
  --stelar-opts "..."   Compatibility alias for --opts
  --reference-species-tree  Reference species tree for RF rate calculation
  --threads, --num-threads, -t, -T
                         CPU worker threads
  --search-space S1..S3  Search-space preset
  --intersection-method I1..I4
                         Intersection method preset
  --keep-polytomy-during-inference
                         Preserve input polytomies during inference; final scoring
                         always preserves the input topology
  --taxa-file FILE       Restrict inference or scoring to listed taxa (one per line)
  --log-file FILE        Save run messages to FILE (progress remains terminal-only)
  --no-time-monitor     Disable time monitoring
  --no-gpu-monitor      Disable GPU monitoring
  --no-notify, -nn      Disable ntfy notifications
  --debug               Enable shell tracing
  --help, -h            Show this message
EOF
}

# Exact wrapper invocation, kept for the reproducibility record written next to
# the output tree (see write_command_record below).
WRAPPER_ARGV=("$0" "$@")

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT_FILE="$2"; shift 2 ;;
    -o|--output) OUTPUT_FILE="$2"; shift 2 ;;
    --stelarx-root|--stelar-root) STELARX_ROOT="$2"; shift 2 ;;
    --opts|--alg-opts|--stelarx-opts|--stelar-opts)
      read -r -a TMP_OPTS <<< "$2"
      STELARX_ARGS+=("${TMP_OPTS[@]}")
      shift 2
      ;;
    --reference-species-tree) REFERENCE_SPECIES_TREE="$2"; shift 2 ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
    --help|-h) print_help; exit 0 ;;
    --auto|--cpu|--gpu|--gpu-strict|--rooted|--unrooted|--keep-polytomy-during-inference|--anchor-outgroup|--anchor|--no-anchor-outgroup|--no-anchor|--no-prune-search-space|--no-gpu-batch|--consensus-experimental|--stepb-fast-restriction|--stepb-quadratic-nn-balls|--stepb-random-leftover-resolution|--stepb-process-large-polytomies|--resolve-input-gene-tree-polytomies|--verify-parse|--verify-hash|--verify-clusters|--verify-partitions|--verify-dp|--verify-weights|--verify-distance-matrix|--verify-similarity-matrix|--verify-upgma|--verify-greedy-consensus|--autocomplete-incomplete-gene-trees|-v|-vv|-vvv|-q|--quiet)
      STELARX_ARGS+=("$1")
      shift
      ;;
    --search-space|--intersection-method|--im|--weight-intersection-method|--search-mode|--log-file|-t|-T|--threads|--num-threads|-m|--seeds|--anchor-taxon|--gpu-batch-size|--gpu-batches|--gpu-vram-control-factor|--gpu-vram-occupancy-factor|--gpu-treewalk-vram-cap-mb|--gpu-progress-interval|--gpu-dp-state-space-construction-output-cap|--gpu-dp-state-space-progress-time-interval|--gpu-dp-state-space-progress-max-steps|--gpu-dist-tile-size|--gpu-sim-vram-cap-mb|--completion-method|--stepb-restriction|--large-n-score-type|--large-score-type|--taxa-file|--species-list|--species-list-file)
      STELARX_ARGS+=("$1" "$2")
      shift 2
      ;;
    --xms|--Xms|--xmx|--Xmx|--no-build)
      STELARX_ARGS+=("$1")
      if [[ "$1" != "--no-build" ]]; then
        STELARX_ARGS+=("$2")
        shift 2
      else
        shift
      fi
      ;;
    --)
      shift
      STELARX_ARGS+=("$@")
      break
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

if [[ ${#POSITIONAL[@]} -gt 0 ]]; then
  echo "Error: positional arguments are not supported."
  print_help
  exit 1
fi

if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
  echo "Error: both --input and --output are required."
  exit 1
fi

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

INPUT_FILE="$(realpath "$INPUT_FILE")"
OUTPUT_FILE="$(realpath -m "$OUTPUT_FILE")"
STELARX_ROOT="$(realpath "$STELARX_ROOT")"
PYTHON_BIN="${STELARX_PYTHON:-${STELARX_ROOT}/.venv/bin/python}"
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN="python3"

if [[ "${DEBUG:-0}" == "1" ]]; then
  set -x
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo -e "${RED}Error: input file '$INPUT_FILE' does not exist.${NC}"
  exit 1
fi
if [[ ! -x "${STELARX_ROOT}/run.sh" ]]; then
  echo -e "${RED}Error: run.sh not found or not executable in '$STELARX_ROOT'.${NC}"
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"

TEMP_DIR="$(mktemp -d)"
TIME_TMP="${TEMP_DIR}/stelarx_time_err.log"
MON_TMP="${TEMP_DIR}/stelarx_gpu_mem.log"
DONE_FILE="${TEMP_DIR}/.stelarx_done"

cleanup() {
  rm -f "$DONE_FILE" 2>/dev/null || true
  if [[ -n "${MON_PID:-}" ]]; then
    kill "$MON_PID" 2>/dev/null || true
    wait "$MON_PID" 2>/dev/null || true
  fi
  rm -rf "$TEMP_DIR" 2>/dev/null || true
}
trap cleanup EXIT

TIME_CMD=""
if [[ "$TIME_MONITOR" == true ]]; then
  if [[ -x "/usr/bin/time" ]]; then
    TIME_CMD="/usr/bin/time"
  elif command -v time >/dev/null 2>&1; then
    TMP_TEST="$(mktemp)"
    sh -c "command time -v true" 2> "$TMP_TEST" >/dev/null || true
    if grep -qi "Maximum resident set size" "$TMP_TEST" 2>/dev/null; then
      TIME_CMD="$(command -v time)"
    fi
    rm -f "$TMP_TEST"
  fi
  if [[ -z "$TIME_CMD" ]]; then
    echo -e "${YELLOW}Warning: no suitable 'time -v' binary found; continuing without time monitor.${NC}"
    TIME_MONITOR=false
  fi
fi

MON_PID=""
if [[ "$GPU_MONITOR" == true ]] && command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
  (
    curmax=0
    while true; do
      gpu_val=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | awk 'BEGIN{m=0} {v=int($1); if(v>m) m=v} END{print m+0}')
      if [[ -n "$gpu_val" && "$gpu_val" =~ ^[0-9]+$ ]] && (( gpu_val > curmax )); then
        curmax=$gpu_val
      fi
      [[ -f "$DONE_FILE" ]] && break
      sleep 0.2
    done
    echo "$curmax" > "$MON_TMP"
  ) &
  MON_PID=$!
else
  GPU_MONITOR=false
fi

# Canonical weight-intersection-method (default prefix-sum) for logs/notifications.
WEIGHT_METHOD="prefix-sum"
for ((wi = 0; wi < ${#STELARX_ARGS[@]}; wi++)); do
  if [[ "${STELARX_ARGS[$wi]}" == "--weight-intersection-method" || "${STELARX_ARGS[$wi]}" == "--intersection-method" || "${STELARX_ARGS[$wi]}" == "--im" ]] && (( wi + 1 < ${#STELARX_ARGS[@]} )); then
    case "${STELARX_ARGS[$((wi + 1))],,}" in
      1|i1|smaller-side-traversal|smaller_side_traversal|smaller-side|smallerside|legacy) WEIGHT_METHOD="smaller-side-traversal" ;;
      2|i2|prefix-sum|prefix_sum|prefixsum|prefix)                                        WEIGHT_METHOD="prefix-sum" ;;
      3|i3|simple-tree-walk|tree-walk|treewalk|simple)                                    WEIGHT_METHOD="simple-tree-walk" ;;
      4|i4|bitset|bitsets|bit-set)                                                        WEIGHT_METHOD="bitset" ;;
      *)                                                                            WEIGHT_METHOD="${STELARX_ARGS[$((wi + 1))]}" ;;
    esac
  fi
done

echo "=== STELAR-X Monitor Wrapper ==="
echo "Input file:     $INPUT_FILE"
echo "Output file:    $OUTPUT_FILE"
echo "STELAR-X root:  $STELARX_ROOT"
echo "STELAR-X opts:  ${STELARX_ARGS[*]:-(defaults)}"
echo "Weight method:  $WEIGHT_METHOD"
if [[ -n "$REFERENCE_SPECIES_TREE" ]]; then
  echo "Reference tree: $REFERENCE_SPECIES_TREE"
fi
echo "Time monitor:   $TIME_MONITOR"
echo "GPU monitor:    $GPU_MONITOR"
echo "Notifications:  $(if [[ "$NO_NOTIFY" == true ]]; then echo "disabled"; else echo "enabled"; fi)"
echo

# Reproducibility record: the exact STELAR-X command (absolute paths, every
# flag) plus the code revision and the wrapper invocation, saved beside the
# output tree as <output>.command. Written before the run so it survives a
# crash; the exit code is appended afterwards.
COMMAND_FILE="${OUTPUT_FILE%.*}.command"
quote_command() {
  local out="" arg
  for arg in "$@"; do
    out+="$(printf '%q' "$arg") "
  done
  printf '%s' "${out% }"
}
write_command_record() {
  local git_commit="unavailable"
  if git -C "$STELARX_ROOT" rev-parse --short HEAD >/dev/null 2>&1; then
    git_commit="$(git -C "$STELARX_ROOT" rev-parse --short HEAD 2>/dev/null)"
    if [[ -n "$(git -C "$STELARX_ROOT" status --porcelain --untracked-files=no 2>/dev/null)" ]]; then
      git_commit+=" (with uncommitted changes)"
    fi
  fi
  {
    echo "# STELAR-X run command (written by run-stelarx-with-monitor.sh)"
    echo "# date:         $(date '+%Y-%m-%dT%H:%M:%S%z')"
    echo "# host:         $(hostname 2>/dev/null || echo unknown)"
    echo "# stelarx_root: $STELARX_ROOT"
    echo "# git_commit:   $git_commit"
    echo "# input:        $INPUT_FILE"
    echo "# output:       $OUTPUT_FILE"
    if [[ -n "$REFERENCE_SPECIES_TREE" ]]; then
      echo "# reference:    $REFERENCE_SPECIES_TREE"
    fi
    echo "# wrapper:      $(quote_command "${WRAPPER_ARGV[@]}")"
    echo "# exact STELAR-X invocation (run from stelarx_root):"
    echo "cd $(printf '%q' "$STELARX_ROOT") && $(quote_command ./run.sh --input "$INPUT_FILE" --output "$OUTPUT_FILE" "${STELARX_ARGS[@]}")"
  } > "$COMMAND_FILE" 2>/dev/null || echo -e "${YELLOW}Warning: could not write command record to $COMMAND_FILE${NC}"
}
write_command_record

START_NS=$(date +%s%N)

STELARX_PID=""
if [[ "$TIME_MONITOR" == true && -n "$TIME_CMD" ]]; then
  (
    cd "$STELARX_ROOT" && "$TIME_CMD" -v ./run.sh --input "$INPUT_FILE" --output "$OUTPUT_FILE" "${STELARX_ARGS[@]}" < /dev/null 2>&1 | tee "$TIME_TMP"
  ) &
  STELARX_PID=$!
else
  (
    cd "$STELARX_ROOT" && ./run.sh --input "$INPUT_FILE" --output "$OUTPUT_FILE" "${STELARX_ARGS[@]}" < /dev/null 2>&1 | tee "$TIME_TMP"
  ) &
  STELARX_PID=$!
fi

set +e
wait "$STELARX_PID"
STELARX_EXIT_CODE=$?
set -e
touch "$DONE_FILE"

END_NS=$(date +%s%N)
ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

if [[ -n "${MON_PID:-}" ]]; then
  wait "$MON_PID" 2>/dev/null || true
fi

MAX_GPU_VAL="NA"
if [[ -f "$MON_TMP" ]]; then
  MAX_GPU_VAL="$(cat "$MON_TMP" 2>/dev/null || echo "NA")"
fi
if [[ "$MAX_GPU_VAL" =~ ^[0-9]+$ ]]; then
  MAX_GPU_MB=$(awk "BEGIN {printf \"%.3f\", ${MAX_GPU_VAL} * 1.024}")
else
  MAX_GPU_MB="NA"
fi

MAX_CPU_MB="NA"
if [[ -f "$TIME_TMP" && -s "$TIME_TMP" ]]; then
  MAX_RSS_KB=$(grep -i "Maximum resident set size" "$TIME_TMP" 2>/dev/null | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | awk '{print int($1)}' | head -n1 || true)
  if [[ -n "${MAX_RSS_KB:-}" && "$MAX_RSS_KB" =~ ^[0-9]+$ ]]; then
    MAX_CPU_MB=$(awk "BEGIN {printf \"%.3f\", ${MAX_RSS_KB}/1024}")
  fi
fi

OPTIMAL_TRIPLET_SCORE="NA"
if [[ -f "$TIME_TMP" ]]; then
  SCORE_LINE=$(sed $'s/\033\\[[0-9;]*m//g' "$TIME_TMP" 2>/dev/null \
    | grep -E 'TRIPLET_SCORE:|[Ff]inal triplet score[[:space:]]*=|Triplet score[[:space:]]+[0-9]+' \
    | tail -n1 || true)
  if [[ -n "$SCORE_LINE" ]]; then
    SCORE_VALUE=$(sed -E -n \
      -e 's/.*TRIPLET_SCORE:[[:space:]]*([0-9]+([.][0-9]+)?).*/\1/p' \
      -e 's/.*[Ff]inal triplet score[[:space:]]*=[[:space:]]*([0-9]+([.][0-9]+)?).*/\1/p' \
      -e 's/.*Triplet score[[:space:]]+([0-9]+([.][0-9]+)?).*/\1/p' \
      <<< "$SCORE_LINE" | tail -n1)
    [[ -n "$SCORE_VALUE" ]] && OPTIMAL_TRIPLET_SCORE="$SCORE_VALUE"
  fi
fi

RF_RATE="NA"
if [[ -n "$REFERENCE_SPECIES_TREE" && -f "$OUTPUT_FILE" ]]; then
  REFERENCE_SPECIES_TREE="$(realpath "$REFERENCE_SPECIES_TREE")"
  if [[ -f "$REFERENCE_SPECIES_TREE" ]]; then
    rf_output=$(cd "$STELARX_ROOT" && "$PYTHON_BIN" rf.py "$OUTPUT_FILE" "$REFERENCE_SPECIES_TREE" 2>&1) || true
    rf_line=$(echo "$rf_output" | grep -i "Robinson-Foulds distance" | tail -n1 || true)
    if [[ -n "$rf_line" ]]; then
      RF_RATE=$(echo "$rf_line" | grep -Eo '[0-9]+(\.[0-9]+)?' | tail -n1 || echo "NA")
    fi
  else
    echo -e "${YELLOW}Warning: reference species tree not found at '$REFERENCE_SPECIES_TREE'; skipping RF.${NC}"
  fi
fi

echo
echo -e "${GREEN}=== STELAR-X Execution Summary ===${NC}"
echo "Status:         $(if [[ $STELARX_EXIT_CODE -eq 0 ]]; then echo -e "${GREEN}SUCCESS${NC}"; else echo -e "${RED}FAILED (exit code $STELARX_EXIT_CODE)${NC}"; fi)"
echo "Running time:   ${RUNNING_TIME}s"
echo "Max CPU RAM:    ${MAX_CPU_MB} MB"
echo "Max GPU VRAM:   ${MAX_GPU_MB} MB"
echo "Triplet score:  ${OPTIMAL_TRIPLET_SCORE}"
if [[ -n "$REFERENCE_SPECIES_TREE" ]]; then
  echo "RF rate:        ${RF_RATE}"
fi
echo "Output exists:  $(if [[ -f "$OUTPUT_FILE" ]]; then echo "Yes"; else echo "No"; fi)"

STATS_FILE="${OUTPUT_FILE%.*}_stats.csv"
echo "algorithm,input_file,output_file,running_time_s,max_cpu_mb,max_gpu_mb,optimal_triplet_score,rf_rate,exit_code" > "$STATS_FILE"
echo "stelar-x,$(basename "$INPUT_FILE"),$(basename "$OUTPUT_FILE"),${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${OPTIMAL_TRIPLET_SCORE},${RF_RATE},${STELARX_EXIT_CODE}" >> "$STATS_FILE"
echo "Stats saved to: $STATS_FILE"
if [[ -f "$COMMAND_FILE" ]]; then
  {
    echo "# exit_code:    ${STELARX_EXIT_CODE}"
    echo "# running_time: ${RUNNING_TIME}s"
  } >> "$COMMAND_FILE" 2>/dev/null || true
  echo "Command saved to: $COMMAND_FILE"
fi

if [[ "$NO_NOTIFY" == false ]] && command -v curl >/dev/null 2>&1; then
  STATUS_EMOJI=$(if [[ $STELARX_EXIT_CODE -eq 0 ]]; then echo "✅"; else echo "❌"; fi)
  STATUS_TEXT=$(if [[ $STELARX_EXIT_CODE -eq 0 ]]; then echo "completed"; else echo "failed (exit $STELARX_EXIT_CODE)"; fi)
  NOTIFY_BODY="${STATUS_EMOJI} STELAR-X ${STATUS_TEXT}

Weight method: ${WEIGHT_METHOD}
Running time: ${RUNNING_TIME}s
Max CPU RAM: ${MAX_CPU_MB} MB
Max GPU VRAM: ${MAX_GPU_MB} MB
Triplet score: ${OPTIMAL_TRIPLET_SCORE}"
  if [[ -n "$REFERENCE_SPECIES_TREE" ]]; then
    NOTIFY_BODY+="
RF rate: ${RF_RATE}"
  fi
  NOTIFY_BODY+="

Input: $(basename "$INPUT_FILE")
Output: $(basename "$OUTPUT_FILE")
Stats: $(basename "$STATS_FILE")"
  curl -s -d "$NOTIFY_BODY" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
fi

exit "$STELARX_EXIT_CODE"
