#!/usr/bin/env bash
# test-stelar-default-monitor.sh
# Defaults to enabling both time and GPU monitoring, but falls back gracefully if
# a suitable `time -v` binary isn't available.
# Usage examples:
#   ./test-stelar-default-monitor.sh -t 1000 -g 500
#   ./test-stelar-default-monitor.sh -t 100 -g 100 --fresh --debug
#
set -euo pipefail

# Defaults
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
BASE_DIR="${HOME}/phylogeny"
SIMPHY_DIR=""                # derived from BASE_DIR unless provided
SIMPHY_DIR_SET=false
SIMPHY_DATA_DIR=""           # custom simphy data directory (overrides SIMPHY_DIR/data)
SIMPHY_DATA_DIR_SET=false
STELAR_ROOT=""               # derived from BASE_DIR unless provided
STELAR_ROOT_SET=false

# Defaults that match run_simulator.sh
SB="0.000001"
SPMIN="500000"
SPMAX="1500000"

USE_LEGACY_LAYOUT=false
STELAR_OPTS="GPU_PARALLEL NONE"
FRESH=false

# Monitoring options (DEFAULT: ON)
TIME_MONITOR=true     # when true: run `time -v` if available and capture stderr to TIME_TMP
GPU_MONITOR=true      # when true: sample nvidia-smi while stelar runs
NO_NOTIFY=false       # when true: skip ntfy.sh notifications
DEBUG=0               # set DEBUG=1 to enable set -x

print_help() {
  cat <<EOF
test-stelar-default-monitor.sh

Required:
  --taxa_num, -t     Number of taxa (e.g. 1000)
  --gene_trees, -g   Number of gene trees (e.g. 500)

Optional:
  --replicate, -r    Replicate name (default: R1)
  --base-dir, -b     Base directory (default: ${BASE_DIR})
  --simphy-dir       Path to simphy dir (overrides --base-dir)
  --simphy-data-dir  Custom directory for simphy data storage (overrides simphy-dir/data)
  --stelar-root      Path to STELAR-X root (overrides --base-dir)
  --stelar-opts      Extra args for STELAR run (default: "$STELAR_OPTS")
  --sb               Substitution/birthrate parameter (default: ${SB})
  --spmin            Population size minimum (default: ${SPMIN})
  --spmax            Population size maximum (default: ${SPMAX})
  --use-legacy-layout  Use legacy simphy layout
  --fresh            Force rerun even if stat-stelar.csv exists
  --no-time-monitor  Disable time-monitoring (overrides default ON)
  --no-gpu-monitor   Disable GPU-monitoring (overrides default ON)
  --no-notify, -nn   Disable ntfy.sh notifications
  --debug            Enable shell tracing (set DEBUG=1)
  --help, -h         Show this message
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxa_num|-t) TAXA_NUM="$2"; shift 2 ;;
    --gene_trees|-g) GENE_TREES="$2"; shift 2 ;;
    --replicate|-r) REPLICATE="$2"; shift 2 ;;
    --simphy-dir) SIMPHY_DIR="$2"; SIMPHY_DIR_SET=true; shift 2 ;;
    --simphy-data-dir) SIMPHY_DATA_DIR="$2"; SIMPHY_DATA_DIR_SET=true; shift 2 ;;
    --stelar-root) STELAR_ROOT="$2"; STELAR_ROOT_SET=true; shift 2 ;;
    --stelar-opts) STELAR_OPTS="$2"; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --sb) SB="$2"; shift 2 ;;
    --spmin) SPMIN="$2"; shift 2 ;;
    --spmax) SPMAX="$2"; shift 2 ;;
    --use-legacy-layout) USE_LEGACY_LAYOUT=true; shift ;;
    --fresh) FRESH=true; shift ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$TAXA_NUM" || -z "$GENE_TREES" ]]; then
  echo "Error: --taxa_num and --gene_trees are required."
  print_help
  exit 2
fi

# Derive SIMPHY_DIR/STELAR_ROOT from BASE_DIR if not explicitly set
if [[ "$SIMPHY_DIR_SET" = false ]]; then
  SIMPHY_DIR="${BASE_DIR%/}/STELAR-X/simphy"
  # SIMPHY_DIR="./simphy"
fi
if [[ "$STELAR_ROOT_SET" = false ]]; then
  STELAR_ROOT="${BASE_DIR%/}/STELAR-X"
  # STELAR_ROOT="."
fi

PAIR="${TAXA_NUM}_${GENE_TREES}"

# Choose simphy run dir using SIMPHY_DATA_DIR if provided, otherwise use SIMPHY_DIR/data
if [[ "$SIMPHY_DATA_DIR_SET" = true ]]; then
  if [[ "$USE_LEGACY_LAYOUT" = true ]]; then
    SIMPHY_RUN_DIR="${SIMPHY_DATA_DIR%/}/${PAIR}/${REPLICATE}"
  else
    SIMPHY_RUN_DIR="${SIMPHY_DATA_DIR%/}/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}/${REPLICATE}"
  fi
else
  if [[ "$USE_LEGACY_LAYOUT" = true ]]; then
    SIMPHY_RUN_DIR="${SIMPHY_DIR%/}/data/${PAIR}/${REPLICATE}"
  else
    SIMPHY_RUN_DIR="${SIMPHY_DIR%/}/data/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}/${REPLICATE}"
  fi
fi

STAT_FILE="${SIMPHY_RUN_DIR%/}/stat-stelar.csv"
ALL_GT_FILE="${SIMPHY_RUN_DIR%/}/all_gt.tre"
TRUE_SPECIES_TREE="${SIMPHY_RUN_DIR%/}/s_tree.trees"
OUT_STELAR="${SIMPHY_RUN_DIR%/}/out-stelar.tre"

# Debug/tracing
if [[ "${DEBUG:-0}" = "1" ]]; then
  set -x
fi

# checkpoint: if stat file exists and --fresh not provided, skip everything
if [[ "$FRESH" = false && -f "${STAT_FILE}" ]]; then
  echo "SKIPPING: ${STAT_FILE} already exists. Use --fresh to force rerun."
  exit 0
fi

echo "Parameters:"
echo "  taxa_num:       $TAXA_NUM"
echo "  gene_trees:     $GENE_TREES"
echo "  replicate:      $REPLICATE"
if [[ "$SIMPHY_DATA_DIR_SET" = true ]]; then
  echo "  simphy data dir: $SIMPHY_DATA_DIR (custom)"
else
  echo "  simphy data dir: ${SIMPHY_DIR%/}/data (default)"
fi
echo "  simphy run dir: $SIMPHY_RUN_DIR"
echo "  out stelar:     $OUT_STELAR"
echo "  stat file:      $STAT_FILE"
echo

if [[ ! -f "$ALL_GT_FILE" ]]; then
  echo "Error: gene-tree file not found at $ALL_GT_FILE"
  exit 6
fi

echo "==> Running STELAR (output will be written to $OUT_STELAR)"

mkdir -p "${SIMPHY_RUN_DIR%/}"

# create log paths inside run dir so they're easy to inspect remotely
TIME_TMP="${SIMPHY_RUN_DIR%/}/.stelar_time_err.log"
MON_TMP="${SIMPHY_RUN_DIR%/}/.stelar_gpu_mem.log"

# Ensure old logs are removed
rm -f "$TIME_TMP" "$MON_TMP" || true

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
    echo "Warning: time-monitor requested (default) but no suitable 'time -v' binary found."
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
      gpu_val=$(nvidia-smi --query-compute-apps=used_gpu_memory --format=csv,noheader,nounits 2>/dev/null | awk 'BEGIN{m=0} {v=int($1); if(v>m) m=v} END{print m+0}')
      if [[ -n "$gpu_val" && "$gpu_val" =~ ^[0-9]+$ ]]; then
        if (( gpu_val > curmax )); then
          curmax=$gpu_val
        fi
      fi
      if [[ -f "${SIMPHY_RUN_DIR%/}/.stelar_done" ]]; then
        break
      fi
      sleep 0.2
    done
    echo "$curmax" > "$MON_TMP"
  ) &
  MON_PID=$!
else
  if [[ "$GPU_MONITOR" = true ]]; then
    echo "Warning: GPU monitor requested but nvidia-smi not found or not executable. Skipping GPU monitor."
    GPU_MONITOR=false
  fi
fi

# Launch STELAR -- prefer using TIME_CMD if available, otherwise run directly
STELAR_PID=""
if [[ "${TIME_MONITOR:-false}" = true && -n "$TIME_CMD" ]]; then
  (
    cd "$STELAR_ROOT" && "$TIME_CMD" -v ./run.sh "$ALL_GT_FILE" "$OUT_STELAR" $STELAR_OPTS < /dev/null
  ) 2> "$TIME_TMP" &
  STELAR_PID=$!
else
  (
    cd "$STELAR_ROOT" && ./run.sh "$ALL_GT_FILE" "$OUT_STELAR" $STELAR_OPTS < /dev/null
  ) &
  STELAR_PID=$!
fi

# Give a short moment and verify the job didn't die immediately
sleep 0.25
if ! kill -0 "$STELAR_PID" >/dev/null 2>&1; then
  echo "Error: STELAR process (pid ${STELAR_PID}) failed to start or died immediately."
  echo "----- Captured time/generic stderr (first 200 lines) -----"
  head -n 200 "$TIME_TMP" 2>/dev/null || true
  echo "---------------------------------------------------------"
  touch "${SIMPHY_RUN_DIR%/}/.stelar_done"
  if [[ -n "${MON_PID:-}" ]]; then
    wait "$MON_PID" 2>/dev/null || true
  fi
  exit 5
fi

echo "STELAR started with PID ${STELAR_PID} (logging to ${TIME_TMP} if time-monitor enabled)"

# Wait for stelar to finish
wait "$STELAR_PID"
STELAR_EXIT_CODE=$?

# Stop GPU monitor by placing sentinel file
touch "${SIMPHY_RUN_DIR%/}/.stelar_done"

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

# If time-monitor wasn't used or cpu RSS not parsed, try a best-effort read of memory using ps for the stelar PID (peak not available):
if [[ "$TIME_MONITOR" = false && "$MAX_CPU_MB" == "NA" ]]; then
  if ps -p "$STELAR_PID" >/dev/null 2>&1; then
    MAX_CPU_MB=$(ps -o rss= -p "$STELAR_PID" 2>/dev/null | awk '{print ($1+0)/1024}')
  else
    MAX_CPU_MB="NA"
  fi
fi

# cleanup small sentinel
rm -f "${SIMPHY_RUN_DIR%/}/.stelar_done" 2>/dev/null || true

echo "STELAR finished in ${RUNNING_TIME}s (exit code ${STELAR_EXIT_CODE})"
echo "Max CPU RAM (MB): ${MAX_CPU_MB}"
echo "Max GPU VRAM (MB): ${MAX_GPU_MB}"

# RF calculation (if rf.py exists and true species tree present)
RF_RATE="NA"
if [[ -f "$OUT_STELAR" && -f "$TRUE_SPECIES_TREE" ]]; then
  echo
  echo "==> Calculating RF rate (using rf.py)"
  rf_output=$(cd "$STELAR_ROOT" && python rf.py "$OUT_STELAR" "$TRUE_SPECIES_TREE" 2>&1) || rf_output="$rf_output"
  rf_candidate=$(echo "$rf_output" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
  if [[ -n "$rf_candidate" ]]; then
    RF_RATE="$rf_candidate"
    echo ""
    echo "STELAR RF rate: ${RF_RATE}"
    echo ""
  else
    RF_RATE="NA"
    echo "Warning: couldn't parse RF rate from rf.py output. rf.py output:"
    echo "-----"
    echo "$rf_output"
    echo "-----"
  fi
else
  echo
  echo "STELAR output or true species tree missing; skipping STELAR RF."
fi

# Write CSV (overwrite every run)
mkdir -p "$(dirname "$STAT_FILE")"
echo "alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb" > "$STAT_FILE"
CSV_ROW="stelar,${TAXA_NUM},${GENE_TREES},${REPLICATE},${SB},${SPMIN},${SPMAX},${RF_RATE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB}"
echo "$CSV_ROW" >> "$STAT_FILE"

echo "Wrote stats to $STAT_FILE"



# Notification disabled

# # Send notification (ntfy)
# if [[ "$NO_NOTIFY" = false ]] && command -v curl >/dev/null 2>&1; then
#   curl -s -d "🎉 STELAR completed for ${TAXA_NUM} taxa and ${GENE_TREES} gene trees!

# 📊 Results:
# alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb
# $CSV_ROW

# 📁 Stats saved to: $STAT_FILE" ntfy.sh/anik-test || true
# fi



# Nicely display stat-stelar.csv summary (same as before)
if [[ -f "${STAT_FILE}" ]]; then
  echo
  echo "=== STELAR run summary (from ${STAT_FILE}) ==="

  # read header and first data row
  IFS= read -r header_line < <(head -n1 "$STAT_FILE")
  IFS= read -r data_line   < <(sed -n '2p' "$STAT_FILE" || true)

  # split into arrays on comma
  IFS=, read -r -a headers <<< "$header_line"
  IFS=, read -r -a values  <<< "$data_line"

  # determine column width for labels
  maxlabel=0
  for h in "${headers[@]}"; do
    len=${#h}
    (( len > maxlabel )) && maxlabel=$len
  done
  (( maxlabel < 16 )) && maxlabel=16

  for i in "${!headers[@]}"; do
    label="${headers[$i]}"
    value="${values[$i]:-}"
    printf "  %-*s : %s\n" "$maxlabel" "$label" "$value"
  done

  echo
  echo "Raw CSV:"
  sed -n '1,2p' "$STAT_FILE" | sed -e 's/,/, /g'
  echo "============================================="
else
  echo
  echo "Warning: stat file not found at ${STAT_FILE}"
fi

echo "Done."
exit 0

