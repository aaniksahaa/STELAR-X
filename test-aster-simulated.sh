#!/usr/bin/env bash
# test-aster-simulated.sh
# Runs ASTER on simulated data with time and GPU monitoring
# Usage examples:
#   ./test-aster-simulated.sh -t 1000 -g 500
#   ./test-aster-simulated.sh -t 100 -g 100 --aster-opts="-t 32" --fresh --debug
#
set -euo pipefail

# Defaults
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
BASE_DIR=".."
SIMPHY_DIR=""                # derived from BASE_DIR unless provided
SIMPHY_DIR_SET=false
SIMPHY_DATA_DIR=""           # custom simphy data directory (overrides SIMPHY_DIR/data)
SIMPHY_DATA_DIR_SET=false
STELAR_ROOT=""               # derived from BASE_DIR unless provided (for ASTER binary location)
STELAR_ROOT_SET=false

# Defaults that match run_simulator.sh
SB="0.000001"
SPMIN="500000"
SPMAX="1500000"

USE_LEGACY_LAYOUT=false
FRESH=false

# ASTER options (passed directly to ASTER binary)
ASTER_OPTS=""

# ASTER binary path (relative to STELAR_ROOT)
ASTER_BIN="ASTER/bin/astral4"

# Monitoring options (DEFAULT: ON)
TIME_MONITOR=true     # when true: run `time -v` if available and capture stderr
GPU_MONITOR=true      # when true: sample nvidia-smi while aster runs
NO_NOTIFY=false       # when true: skip ntfy.sh notifications
DEBUG=0               # set DEBUG=1 to enable set -x

print_help() {
  cat <<EOF
test-aster-simulated.sh

Required:
  --taxa_num, -t     Number of taxa (e.g. 1000)
  --gene_trees, -g   Number of gene trees (e.g. 500)

Optional:
  --replicate, -r    Replicate name (default: R1)
  --base-dir, -b     Base directory (default: ${BASE_DIR})
  --simphy-dir       Path to simphy dir (overrides --base-dir)
  --simphy-data-dir  Custom directory for simphy data storage (overrides simphy-dir/data)
  --stelar-root      Path to STELAR-X root (overrides --base-dir)
  --aster-opts       Raw options to pass to ASTER (e.g. --aster-opts="-t 32 -u 2")
  --aster-bin        Path to ASTER binary relative to STELAR_ROOT (default: ASTER/bin/astral4)
  --sb               Substitution/birthrate parameter (default: ${SB})
  --spmin            Population size minimum (default: ${SPMIN})
  --spmax            Population size maximum (default: ${SPMAX})
  --use-legacy-layout  Use legacy simphy layout
  --fresh            Force rerun even if stat-aster.csv exists
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
    --aster-opts=*) ASTER_OPTS="${1#*=}"; shift ;;
    --aster-opts) ASTER_OPTS="$2"; shift 2 ;;
    --aster-bin) ASTER_BIN="$2"; shift 2 ;;
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
fi
if [[ "$STELAR_ROOT_SET" = false ]]; then
  STELAR_ROOT="${BASE_DIR%/}/STELAR-X"
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

STAT_FILE="${SIMPHY_RUN_DIR%/}/stat-aster.csv"
ALL_GT_FILE="${SIMPHY_RUN_DIR%/}/all_gt.tre"
TRUE_SPECIES_TREE="${SIMPHY_RUN_DIR%/}/s_tree.trees"
OUT_ASTER="${SIMPHY_RUN_DIR%/}/out-aster.tre"
LOCK_FILE="${SIMPHY_RUN_DIR%/}/.aster.lock"

# Full path to ASTER binary
ASTER_BIN_PATH="${STELAR_ROOT%/}/${ASTER_BIN}"

# Debug/tracing
if [[ "${DEBUG:-0}" = "1" ]]; then
  set -x
fi

# checkpoint: if lock file exists and --fresh not provided, skip everything
if [[ "$FRESH" = false && -f "${LOCK_FILE}" ]]; then
  echo "SKIPPING: ${LOCK_FILE} already exists. Use --fresh to force rerun."
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
echo "  out aster:      $OUT_ASTER"
echo "  stat file:      $STAT_FILE"
echo "  aster binary:   $ASTER_BIN_PATH"
echo "  aster opts:     ${ASTER_OPTS:-<none>}"
echo

if [[ ! -f "$ALL_GT_FILE" ]]; then
  echo "Error: gene-tree file not found at $ALL_GT_FILE"
  exit 6
fi

if [[ ! -x "$ASTER_BIN_PATH" ]]; then
  echo "Error: ASTER binary not found or not executable at $ASTER_BIN_PATH"
  exit 7
fi

echo "==> Running ASTER (output will be written to $OUT_ASTER)"

mkdir -p "${SIMPHY_RUN_DIR%/}"

# create log paths inside run dir so they're easy to inspect remotely
TIME_TMP="${SIMPHY_RUN_DIR%/}/.aster_time_err.log"
MON_TMP="${SIMPHY_RUN_DIR%/}/.aster_gpu_mem.log"

# Ensure old logs are removed
rm -f "$TIME_TMP" "$MON_TMP" || true

# Remove existing output file to ensure we can verify fresh output
rm -f "$OUT_ASTER" || true

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
      if [[ -f "${SIMPHY_RUN_DIR%/}/.aster_done" ]]; then
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

# Build ASTER command with input/output plus any extra opts
ASTER_CMD="\"$ASTER_BIN_PATH\" -i \"$ALL_GT_FILE\" -o \"$OUT_ASTER\""
if [[ -n "$ASTER_OPTS" ]]; then
  ASTER_CMD="$ASTER_CMD $ASTER_OPTS"
fi

# Launch ASTER -- run in foreground so output is visible
# Use time -v -o to write time stats to file while showing ASTER output on terminal
echo "----------------------------------------"
set +e
if [[ "${TIME_MONITOR:-false}" = true && -n "$TIME_CMD" ]]; then
  eval "$TIME_CMD -v -o \"$TIME_TMP\" $ASTER_CMD"
  ASTER_EXIT_CODE=$?
else
  eval "$ASTER_CMD"
  ASTER_EXIT_CODE=$?
fi
set -e
echo "----------------------------------------"

# Stop GPU monitor by placing sentinel file
touch "${SIMPHY_RUN_DIR%/}/.aster_done"

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

# If time-monitor wasn't used or cpu RSS not parsed, try a best-effort read of memory using ps for the aster PID (peak not available):
if [[ "$TIME_MONITOR" = false && "$MAX_CPU_MB" == "NA" ]]; then
  if ps -p "$ASTER_PID" >/dev/null 2>&1; then
    MAX_CPU_MB=$(ps -o rss= -p "$ASTER_PID" 2>/dev/null | awk '{print ($1+0)/1024}')
  else
    MAX_CPU_MB="NA"
  fi
fi

# cleanup small sentinel
rm -f "${SIMPHY_RUN_DIR%/}/.aster_done" 2>/dev/null || true

echo "ASTER finished in ${RUNNING_TIME}s (exit code ${ASTER_EXIT_CODE})"
echo "Max CPU RAM (MB): ${MAX_CPU_MB}"
echo "Max GPU VRAM (MB): ${MAX_GPU_MB}"

# RF calculation (if rf.py exists and true species tree present)
RF_RATE="NA"
if [[ -f "$OUT_ASTER" && -f "$TRUE_SPECIES_TREE" ]]; then
  echo
  echo "==> Calculating RF rate (using rf.py)"
  rf_output=$(cd "$STELAR_ROOT" && python rf.py "$OUT_ASTER" "$TRUE_SPECIES_TREE" 2>&1) || rf_output="$rf_output"
  rf_candidate=$(echo "$rf_output" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
  if [[ -n "$rf_candidate" ]]; then
    RF_RATE="$rf_candidate"
    echo ""
    echo "ASTER RF rate: ${RF_RATE}"
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
  echo "ASTER output or true species tree missing; skipping ASTER RF."
fi

# Only write stats and lock file if output file exists AND exit code was 0
if [[ -f "$OUT_ASTER" && "$ASTER_EXIT_CODE" -eq 0 ]]; then
  # Write CSV (overwrite every run)
  mkdir -p "$(dirname "$STAT_FILE")"
  echo "alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb" > "$STAT_FILE"
  CSV_ROW="aster,${TAXA_NUM},${GENE_TREES},${REPLICATE},${SB},${SPMIN},${SPMAX},${RF_RATE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB}"
  echo "$CSV_ROW" >> "$STAT_FILE"

  echo "Wrote stats to $STAT_FILE"

  # Create lock file to indicate successful completion
  touch "${LOCK_FILE}"
  echo -e "\033[1;32m✅ Run completed successfully. Lock file created.\033[0m"

  # Send notification (ntfy) - only on success
  if [[ "$NO_NOTIFY" = false ]] && command -v curl >/dev/null 2>&1; then
    curl -s -d "🎉 ASTER completed for ${TAXA_NUM} taxa and ${GENE_TREES} gene trees!

⚙️ Config: sb=${SB} spmin=${SPMIN} spmax=${SPMAX} rep=${REPLICATE}

📊 Results:
RF: ${RF_RATE} | Time: ${RUNNING_TIME}s
CPU: ${MAX_CPU_MB} MB | GPU: ${MAX_GPU_MB} MB

📁 ${STAT_FILE}" https://ntfy.sh/anik-test || true
  fi
else
  echo -e "\033[1;31m❌ Run failed or was interrupted (exit code: ${ASTER_EXIT_CODE}, output exists: $(test -f "$OUT_ASTER" && echo yes || echo no))\033[0m"
  echo "Skipping stat/lock file creation."
fi

# Nicely display stat-aster.csv summary
if [[ -f "${STAT_FILE}" ]]; then
  echo
  echo "=== ASTER run summary (from ${STAT_FILE}) ==="

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
