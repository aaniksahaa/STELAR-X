#!/usr/bin/env bash
# test_astral.sh
# Usage:
#   ./test_astral.sh -t 1000 -g 500
#   ./test_astral.sh --taxa_num 1000 --gene_trees 500 --replicate R2 --astral-root /path/to/ASTRAL --fresh

set -euo pipefail

# Defaults
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
BASE_DIR="${HOME}/phylogeny"
ASTRAL_DIR=""                # derived from BASE_DIR unless provided
ASTRAL_DIR_SET=false
ASTRAL_ROOT=""               # alias for ASTRAL_DIR for compatibility
ASTRAL_ROOT_SET=false

# Defaults that match run_simulator.sh (kept for directory naming parity)
SB="0.000001"
SPMIN="500000"
SPMAX="1500000"

USE_LEGACY_LAYOUT=false
ASTRAL_OPTS=""
FRESH=false

print_help() {
  cat <<EOF
test_astral.sh

Required:
  --taxa_num, -t     Number of taxa (e.g. 1000)
  --gene_trees, -g   Number of gene trees (e.g. 500)

Optional:
  --replicate, -r    Replicate name (default: R1)
  --base-dir, -b     Base directory (default: ${BASE_DIR})
  --astral-dir        Path to ASTRAL dir (overrides --base-dir)
  --astral-root       Path to ASTRAL root (alias for --astral-dir)
  --astral-opts      Extra args for ASTRAL run (default: "$ASTRAL_OPTS")
  --sb               Substitution/birthrate parameter (default: ${SB})
  --spmin            Population size minimum (default: ${SPMIN})
  --spmax            Population size maximum (default: ${SPMAX})
  --use-legacy-layout  Use legacy simphy layout
  --fresh            Force rerun even if stat-astral.csv exists
  --help, -h         Show this message
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxa_num|-t) TAXA_NUM="$2"; shift 2 ;;
    --gene_trees|-g) GENE_TREES="$2"; shift 2 ;;
    --replicate|-r) REPLICATE="$2"; shift 2 ;;
    --astral-dir) ASTRAL_DIR="$2"; ASTRAL_DIR_SET=true; shift 2 ;;
    --astral-root) ASTRAL_ROOT="$2"; ASTRAL_ROOT_SET=true; shift 2 ;;
    --astral-opts) ASTRAL_OPTS="$2"; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --sb) SB="$2"; shift 2 ;;
    --spmin) SPMIN="$2"; shift 2 ;;
    --spmax) SPMAX="$2"; shift 2 ;;
    --use-legacy-layout) USE_LEGACY_LAYOUT=true; shift ;;
    --fresh) FRESH=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$TAXA_NUM" || -z "$GENE_TREES" ]]; then
  echo "Error: --taxa_num and --gene_trees are required."
  print_help
  exit 2
fi

# Derive ASTRAL_DIR/ASTRAL_ROOT from BASE_DIR if not explicitly set
if [[ "$ASTRAL_DIR_SET" = false && "$ASTRAL_ROOT_SET" = false ]]; then
  ASTRAL_ROOT="${BASE_DIR%/}/ASTRAL"
elif [[ "$ASTRAL_DIR_SET" = true && "$ASTRAL_ROOT_SET" = false ]]; then
  ASTRAL_ROOT="$ASTRAL_DIR"
elif [[ "$ASTRAL_ROOT_SET" = true ]]; then
  # ASTRAL_ROOT already set by user
  :
fi

PAIR="${TAXA_NUM}_${GENE_TREES}"

# Construct SIMPHY_RUN_DIR early (so we can check the stat file before doing heavy work)
if [[ "$USE_LEGACY_LAYOUT" = true ]]; then
  SIMPHY_RUN_DIR="${BASE_DIR%/}/STELAR-MP/simphy/data/${PAIR}/${REPLICATE}"
else
  SIMPHY_RUN_DIR="${BASE_DIR%/}/STELAR-MP/simphy/data/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}/${REPLICATE}"
fi

STAT_FILE="${SIMPHY_RUN_DIR%/}/stat-astral.csv"
ALL_GT_FILE="${SIMPHY_RUN_DIR%/}/all_gt.tre"
TRUE_SPECIES_TREE="${SIMPHY_RUN_DIR%/}/s_tree.trees"
OUT_ASTRAL="${SIMPHY_RUN_DIR%/}/out-astral.tre"

# Checkpoint: if stat file exists and --fresh not provided, skip everything
if [[ "$FRESH" = false && -f "${STAT_FILE}" ]]; then
  echo "SKIPPING: ${STAT_FILE} already exists. Use --fresh to force rerun."
  exit 0
fi

echo "Parameters:"
echo "  taxa_num:       $TAXA_NUM"
echo "  gene_trees:     $GENE_TREES"
echo "  replicate:      $REPLICATE"
echo "  simphy run dir: $SIMPHY_RUN_DIR"
echo "  out astral:     $OUT_ASTRAL"
echo "  stat file:      $STAT_FILE"
echo

if [[ ! -f "$ALL_GT_FILE" ]]; then
  echo "Error: gene-tree file not found at $ALL_GT_FILE"
  exit 6
fi

if [[ ! -d "$ASTRAL_ROOT" ]]; then
  echo "Warning: ASTRAL root directory '$ASTRAL_ROOT' does not exist. Continuing, but run may fail."
fi

echo "==> Running ASTRAL (output will be written to $OUT_ASTRAL)"

mkdir -p "${SIMPHY_RUN_DIR%/}"

# prepare tempfiles
TIME_TMP=$(mktemp)
MON_TMP=$(mktemp)

# time + gpu-monitor wrapper
START_NS=$(date +%s%N)

# run ASTRAL under /usr/bin/time -v and capture its stderr (resource usage info)
# (
#   cd "$ASTRAL_ROOT" && /usr/bin/time -v ./run_astral.sh -i "$ALL_GT_FILE" -o "$OUT_ASTRAL" $ASTRAL_OPTS
# ) 2> "$TIME_TMP" &
# ASTRAL_WRAPPER_PID=$!

(
  cd "$ASTRAL_ROOT" && /usr/bin/time -v ./run_astral.sh -i "$ALL_GT_FILE" -o "$OUT_ASTRAL" $ASTRAL_OPTS
) 2> >(tee "$TIME_TMP" >&2) &
ASTRAL_WRAPPER_PID=$!


# start GPU monitor only if nvidia-smi exists (some clusters may have GPUs even if ASTRAL is CPU-bound)
if command -v nvidia-smi >/dev/null 2>&1; then
  (
    curmax=0
    # sample at 0.1s to catch short spikes
    while kill -0 "$ASTRAL_WRAPPER_PID" 2>/dev/null; do
      gpu_val=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | awk 'BEGIN{m=0} {v=int($1); if(v>m) m=v} END{print m}')
      if [[ -n "$gpu_val" && "$gpu_val" =~ ^[0-9]+$ ]]; then
        if (( gpu_val > curmax )); then
          curmax=$gpu_val
        fi
      fi
      sleep 0.1
    done
    echo "$curmax" > "$MON_TMP"
  ) &
  MON_PID=$!
else
  echo "NA" > "$MON_TMP"
  MON_PID=""
fi

# wait for the ASTRAL wrapper to finish and capture exit code
wait "$ASTRAL_WRAPPER_PID"
ASTRAL_EXIT_CODE=$?

END_NS=$(date +%s%N)
ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

# ensure monitor stopped and read its result
if [[ -n "${MON_PID:-}" ]]; then
  wait "$MON_PID" 2>/dev/null || true
fi

MAX_GPU_VAL=$(cat "$MON_TMP" 2>/dev/null || echo "NA")

# Normalize GPU: nvidia-smi reports MiB (integer). Convert to decimal MB with 3 decimal places (MiB * 1.024).
if [[ "$MAX_GPU_VAL" =~ ^[0-9]+$ ]]; then
  MAX_GPU_MB=$(awk "BEGIN{printf \"%.3f\", ${MAX_GPU_VAL} * 1.024}")
else
  MAX_GPU_MB="NA"
fi

# parse /usr/bin/time -v output to get Maximum resident set size (kbytes)
MAX_CPU_MB="NA"
if grep -qi "Maximum resident set size" "$TIME_TMP" 2>/dev/null; then
  MAX_RSS_KB=$(grep -i "Maximum resident set size" "$TIME_TMP" | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | awk '{print int($1)}')
elif grep -qi "Maximum resident set size (kbytes)" "$TIME_TMP" 2>/dev/null; then
  MAX_RSS_KB=$(grep -i "Maximum resident set size (kbytes)" "$TIME_TMP" | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | awk '{print int($1)}')
else
  MAX_RSS_KB=""
fi

if [[ -n "$MAX_RSS_KB" && "$MAX_RSS_KB" =~ ^[0-9]+$ ]]; then
  # convert kB -> MB (1024 kB = 1 MiB) and show 3 decimal places
  MAX_CPU_MB=$(awk "BEGIN{printf \"%.3f\", ${MAX_RSS_KB}/1024}")
fi

# clean up tempfiles
rm -f "$TIME_TMP" "$MON_TMP" 2>/dev/null || true

echo "ASTRAL finished in ${RUNNING_TIME}s (exit code ${ASTRAL_EXIT_CODE})"
echo "Max CPU RAM (MB): ${MAX_CPU_MB}"
echo "Max GPU VRAM (MB): ${MAX_GPU_MB}"

# RF calculation (if rf.py exists and true species tree present)
RF_RATE="NA"
if [[ -f "$OUT_ASTRAL" && -f "$TRUE_SPECIES_TREE" ]]; then
  echo
  echo "==> Calculating RF rate (using rf.py)"
  rf_output=$(cd "$ASTRAL_ROOT" && python rf.py "$OUT_ASTRAL" "$TRUE_SPECIES_TREE" 2>&1) || rf_output="$rf_output"
  rf_candidate=$(echo "$rf_output" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
  if [[ -n "$rf_candidate" ]]; then
    RF_RATE="$rf_candidate"
    echo ""
    echo "ASTRAL RF rate: ${RF_RATE}"
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
  echo "ASTRAL output or true species tree missing; skipping ASTRAL RF."
fi

# Write CSV (overwrite every run) — includes replicate after gene-trees, and max-cpu-mb and max-gpu-mb
echo "alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb" > "$STAT_FILE"
echo "astral,${TAXA_NUM},${GENE_TREES},${REPLICATE},${SB},${SPMIN},${SPMAX},${RF_RATE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB}" >> "$STAT_FILE"

echo "Wrote stats to $STAT_FILE"
echo "Done."
