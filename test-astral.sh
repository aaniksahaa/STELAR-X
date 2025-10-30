#!/usr/bin/env bash
# test_astral.sh (simplified, no CPU/GPU monitoring)
# Usage:
#   ./test_astral.sh -t 1000 -g 500
#   ./test_astral.sh --taxa_num 1000 --gene_trees 500 --replicate R2 --astral-root /path/to/ASTRAL --fresh

set -euo pipefail

# Defaults
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
BASE_DIR="${HOME}/phylogeny"

SIMPHY_DIR=""
SIMPHY_DIR_SET=false
SIMPHY_DATA_DIR=""
SIMPHY_DATA_DIR_SET=false

ASTRAL_DIR=""
ASTRAL_DIR_SET=false
ASTRAL_ROOT=""
ASTRAL_ROOT_SET=false

STELAR_ROOT="${BASE_DIR%/}/STELAR-MP"

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
  --simphy-dir       Path to simphy dir (overrides --base-dir)
  --simphy-data-dir  Custom directory for simphy data storage
  --astral-dir       Path to ASTRAL dir (overrides --base-dir)
  --astral-root      Path to ASTRAL root (alias for --astral-dir)
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
    --simphy-dir) SIMPHY_DIR="$2"; SIMPHY_DIR_SET=true; shift 2 ;;
    --simphy-data-dir) SIMPHY_DATA_DIR="$2"; SIMPHY_DATA_DIR_SET=true; shift 2 ;;
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

# Derive dirs
if [[ "$SIMPHY_DIR_SET" = false ]]; then
  SIMPHY_DIR="${BASE_DIR%/}/STELAR-MP/simphy"
fi

if [[ "$ASTRAL_DIR_SET" = false && "$ASTRAL_ROOT_SET" = false ]]; then
  ASTRAL_ROOT="${BASE_DIR%/}/ASTRAL"
elif [[ "$ASTRAL_DIR_SET" = true && "$ASTRAL_ROOT_SET" = false ]]; then
  ASTRAL_ROOT="$ASTRAL_DIR"
fi

PAIR="${TAXA_NUM}_${GENE_TREES}"

# Construct SIMPHY_RUN_DIR
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

STAT_FILE="${SIMPHY_RUN_DIR%/}/stat-astral.csv"
ALL_GT_FILE="${SIMPHY_RUN_DIR%/}/all_gt.tre"
TRUE_SPECIES_TREE="${SIMPHY_RUN_DIR%/}/s_tree.trees"
OUT_ASTRAL="${SIMPHY_RUN_DIR%/}/out-astral.tre"

# Skip if done
if [[ "$FRESH" = false && -f "${STAT_FILE}" ]]; then
  echo "SKIPPING: ${STAT_FILE} already exists. Use --fresh to force rerun."
  exit 0
fi

echo "Parameters:"
echo "  taxa_num:   $TAXA_NUM"
echo "  gene_trees: $GENE_TREES"
echo "  replicate:  $REPLICATE"
echo "  simphy run: $SIMPHY_RUN_DIR"
echo

if [[ ! -f "$ALL_GT_FILE" ]]; then
  echo "Error: gene-tree file not found at $ALL_GT_FILE"
  exit 6
fi

echo "==> Running ASTRAL (output: $OUT_ASTRAL)"

mkdir -p "${SIMPHY_RUN_DIR%/}"
START_NS=$(date +%s%N)

(
  cd "$ASTRAL_ROOT" && ./run_astral.sh -i "$ALL_GT_FILE" -o "$OUT_ASTRAL" $ASTRAL_OPTS
)
ASTRAL_EXIT_CODE=$?

END_NS=$(date +%s%N)
ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

echo "ASTRAL finished in ${RUNNING_TIME}s (exit code ${ASTRAL_EXIT_CODE})"

# RF calculation
RF_RATE="NA"
if [[ -f "$OUT_ASTRAL" && -f "$TRUE_SPECIES_TREE" ]]; then
  echo "==> Calculating RF rate (using rf.py)"
  if [[ -f "${STELAR_ROOT%/}/rf.py" ]]; then
    rf_output=$(cd "$STELAR_ROOT" && python rf.py "$OUT_ASTRAL" "$TRUE_SPECIES_TREE" 2>&1) || true
    rf_candidate=$(echo "$rf_output" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
    if [[ -n "$rf_candidate" ]]; then
      RF_RATE="$rf_candidate"
      echo "ASTRAL RF rate: ${RF_RATE}"
    else
      echo "Warning: couldn't parse RF rate from rf.py output."
      echo "$rf_output"
    fi
  else
    echo "Warning: rf.py not found; skipping RF calculation."
  fi
fi

# Write CSV (no cpu/gpu columns)
echo "alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s" > "$STAT_FILE"
CSV_ROW="astral,${TAXA_NUM},${GENE_TREES},${REPLICATE},${SB},${SPMIN},${SPMAX},${RF_RATE},${RUNNING_TIME}"
echo "$CSV_ROW" >> "$STAT_FILE"

echo "Wrote stats to $STAT_FILE"

# Notification
curl -s -d "🎉 ASTRAL completed for ${TAXA_NUM} taxa and ${GENE_TREES} gene trees!

📊 Results:
alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s
$CSV_ROW

📁 Stats saved to: $STAT_FILE" ntfy.sh/anik-test >/dev/null

echo "Done."
exit "${ASTRAL_EXIT_CODE}"
