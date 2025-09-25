#!/usr/bin/env bash
# test-stelar-simple.sh
# Simple runner: runs STELAR directly (logs visible), measures elapsed time,
# computes RF (if available), writes CSV. No CPU/GPU monitoring or PID wrapping.
#
# Usage:
#   ./test-stelar-simple.sh -t 1000 -g 500
#   ./test-stelar-simple.sh --taxa_num 1000 --gene_trees 500 --replicate R2 --stelar-root /path/to/STELAR-MP --fresh

set -euo pipefail

# Defaults
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
BASE_DIR="${HOME}/phylogeny"
SIMPHY_DIR=""                # derived from BASE_DIR unless provided
SIMPHY_DIR_SET=false
SIMPHY_DATA_DIR=""           # custom simphy data directory
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

print_help() {
  cat <<EOF
test-stelar-simple.sh

Required:
  --taxa_num, -t     Number of taxa (e.g. 1000)
  --gene_trees, -g   Number of gene trees (e.g. 500)

Optional:
  --replicate, -r    Replicate name (default: R1)
  --base-dir, -b     Base directory (default: ${BASE_DIR})
  --simphy-dir       Path to simphy dir (overrides --base-dir)
  --simphy-data-dir  Custom directory for simphy data storage
  --stelar-root      Path to STELAR-MP root (overrides --base-dir)
  --stelar-opts      Extra args for STELAR run (default: "$STELAR_OPTS")
  --sb               Substitution/birthrate parameter (default: ${SB})
  --spmin            Population size minimum (default: ${SPMIN})
  --spmax            Population size maximum (default: ${SPMAX})
  --use-legacy-layout  Use legacy simphy layout
  --fresh            Force rerun even if stat-stelar.csv exists
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

# Construct SIMPHY_RUN_DIR early (so we can check the stat file before doing heavy work)
if [[ "$SIMPHY_DATA_DIR_SET" = true ]]; then
  # Use custom data directory
  if [[ "$USE_LEGACY_LAYOUT" = true ]]; then
    SIMPHY_RUN_DIR="${SIMPHY_DATA_DIR%/}/${PAIR}/${REPLICATE}"
  else
    SIMPHY_RUN_DIR="${SIMPHY_DATA_DIR%/}/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}/${REPLICATE}"
  fi
else
  # Default behavior: use data directory inside SIMPHY_DIR
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

# Checkpoint: if stat file exists and --fresh not provided, skip everything
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

# measure elapsed time with nanosecond precision
START_NS=$(date +%s%N)

# Run STELAR directly so its stdout/stderr are visible in this runner's terminal.
# We cd into STELAR_ROOT and run run.sh with arguments.
(
  cd "$STELAR_ROOT"
  ./run.sh "$ALL_GT_FILE" "$OUT_STELAR" $STELAR_OPTS
)
STELAR_EXIT_CODE=$?

END_NS=$(date +%s%N)
ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

echo "STELAR finished in ${RUNNING_TIME}s (exit code ${STELAR_EXIT_CODE})"

# RF calculation (if rf.py exists and true species tree present)
RF_RATE="NA"
if [[ -f "$OUT_STELAR" && -f "$TRUE_SPECIES_TREE" && -f "${STELAR_ROOT%/}/rf.py" ]]; then
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
  echo "STELAR output or true species tree or rf.py missing; skipping RF calculation."
fi

# Write CSV (overwrite every run) — simplified header (no cpu/gpu columns)
echo "alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s" > "$STAT_FILE"
CSV_ROW="stelar,${TAXA_NUM},${GENE_TREES},${REPLICATE},${SB},${SPMIN},${SPMAX},${RF_RATE},${RUNNING_TIME}"
echo "$CSV_ROW" >> "$STAT_FILE"

echo "Wrote stats to $STAT_FILE"

# Send notification with results
echo "Sending notification..."
curl -s -d "🎉 STELAR completed for ${TAXA_NUM} taxa and ${GENE_TREES} gene trees!

📊 Results:
alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s
$CSV_ROW

📁 Stats saved to: $STAT_FILE" ntfy.sh/anik-test

# # Cleanup should not be done, because then our csv would also be deleted
# # Cleanup: remove the simphy run directory to save disk space
# echo "Cleaning up simphy run directory: $SIMPHY_RUN_DIR"
# if [[ -d "$SIMPHY_RUN_DIR" ]]; then
#   rm -rf "$SIMPHY_RUN_DIR"
#   echo "✅ Cleanup completed - removed $SIMPHY_RUN_DIR"
# else
#   echo "⚠️  Directory $SIMPHY_RUN_DIR not found, nothing to clean up"
# fi

echo "Done."

# Exit with the same exit code as STELAR (so CI/automation can detect failure)
exit "${STELAR_EXIT_CODE}"
