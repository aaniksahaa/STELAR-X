#!/usr/bin/env bash
# test-stelar.sh
# Usage:
#   ./test-stelar.sh -t 1000 -g 500
#   ./test-stelar.sh --taxa_num 1000 --gene_trees 500 --replicate R2 --stelar-root /path/to/STELAR-MP

set -euo pipefail

# Defaults
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
BASE_DIR="${HOME}/phylogeny"
SIMPHY_DIR=""                # derived from BASE_DIR unless provided
SIMPHY_DIR_SET=false
STELAR_ROOT=""               # derived from BASE_DIR unless provided
STELAR_ROOT_SET=false

# Defaults that match run_simulator.sh
SB="0.000001"
SPMIN="500000"
SPMAX="1500000"

USE_LEGACY_LAYOUT=false
STELAR_OPTS="GPU_PARALLEL NONE"

print_help() {
  cat <<EOF
test-stelar.sh

Required:
  --taxa_num, -t     Number of taxa (e.g. 1000)
  --gene_trees, -g   Number of gene trees (e.g. 500)

Optional:
  --replicate, -r    Replicate name (default: R1)
  --base-dir, -b     Base directory (default: ${BASE_DIR})
  --simphy-dir       Path to simphy dir (overrides --base-dir)
  --stelar-root      Path to STELAR-MP root (overrides --base-dir)
  --stelar-opts      Extra args for STELAR run (default: "$STELAR_OPTS")
  --sb               Substitution/birthrate parameter (default: ${SB})
  --spmin            Population size minimum (default: ${SPMIN})
  --spmax            Population size maximum (default: ${SPMAX})
  --use-legacy-layout  Use legacy simphy layout
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
    --stelar-root) STELAR_ROOT="$2"; STELAR_ROOT_SET=true; shift 2 ;;
    --stelar-opts) STELAR_OPTS="$2"; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --sb) SB="$2"; shift 2 ;;
    --spmin) SPMIN="$2"; shift 2 ;;
    --spmax) SPMAX="$2"; shift 2 ;;
    --use-legacy-layout) USE_LEGACY_LAYOUT=true; shift ;;
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
  SIMPHY_DIR="${BASE_DIR%/}/STELAR-MP/simphy"
fi
if [[ "$STELAR_ROOT_SET" = false ]]; then
  STELAR_ROOT="${BASE_DIR%/}/STELAR-MP"
fi

PAIR="${TAXA_NUM}_${GENE_TREES}"

# Construct SIMPHY_RUN_DIR:
if [[ "$USE_LEGACY_LAYOUT" = true ]]; then
  SIMPHY_RUN_DIR="${SIMPHY_DIR%/}/data/${PAIR}/${REPLICATE}"
else
  SIMPHY_RUN_DIR="${SIMPHY_DIR%/}/data/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}/${REPLICATE}"
fi

ALL_GT_FILE="${SIMPHY_RUN_DIR%/}/all_gt.tre"
TRUE_SPECIES_TREE="${SIMPHY_RUN_DIR%/}/s_tree.trees"

# OUTPUTS now live in the simphy run dir
OUT_STELAR="${SIMPHY_RUN_DIR%/}/out-stelar.tre"
STAT_FILE="${SIMPHY_RUN_DIR%/}/stat-stelar.csv"

echo "Parameters:"
echo "  taxa_num:       $TAXA_NUM"
echo "  gene_trees:     $GENE_TREES"
echo "  replicate:      $REPLICATE"
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

# prepare tempfiles
TIME_TMP=$(mktemp)
MON_TMP=$(mktemp)

# time + gpu-monitor wrapper
START_NS=$(date +%s%N)

# run STELAR under /usr/bin/time -v and capture its stderr (resource usage info)
(
  cd "$STELAR_ROOT" && /usr/bin/time -v ./run.sh "$ALL_GT_FILE" "$OUT_STELAR" $STELAR_OPTS
) 2> "$TIME_TMP" &
STELAR_WRAPPER_PID=$!

# start GPU monitor only if nvidia-smi exists
if command -v nvidia-smi >/dev/null 2>&1; then
  (
    curmax=0
    # sample at 0.1s to catch short spikes
    while kill -0 "$STELAR_WRAPPER_PID" 2>/dev/null; do
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

# wait for the STELAR wrapper to finish and capture exit code
wait "$STELAR_WRAPPER_PID"
STELAR_EXIT_CODE=$?

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

# Write CSV (overwrite every run) — includes replicate after gene-trees, and max-cpu-mb and max-gpu-mb
echo "alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb" > "$STAT_FILE"
echo "stelar,${TAXA_NUM},${GENE_TREES},${REPLICATE},${SB},${SPMIN},${SPMAX},${RF_RATE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB}" >> "$STAT_FILE"

echo "Wrote stats to $STAT_FILE"
echo "Done."
