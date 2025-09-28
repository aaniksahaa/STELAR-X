#!/usr/bin/env bash
# sim.sh (new runner for updated run_simulator.sh)
# Usage examples:
#   ./sim.sh -t 1000 -g 500
#   ./sim.sh -t 1000 -g 500 -r R1 -b "$HOME/phylogeny" --fresh

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
FRESH=false

# Defaults that match run_simulator.sh
SB="0.000001"
SPMIN="500000"
SPMAX="1500000"

print_help() {
  cat <<EOF
sim.sh - runner for run_simulator.sh

Required:
  --taxa_num, -t     Number of taxa (e.g. 1000)
  --gene_trees, -g   Number of gene trees (e.g. 500)

Optional:
  --replicate, -r    Replicate name (default: R1)
  --base-dir, -b     Base directory (default: ${BASE_DIR})
  --simphy-dir       Path to simphy dir (overrides --base-dir)
  --simphy-data-dir  Custom directory for simphy data storage
  --sb               Substitution/birthrate parameter (default: ${SB})
  --spmin            Population size minimum (default: ${SPMIN})
  --spmax            Population size maximum (default: ${SPMAX})
  --fresh            Force rerun even if stat-sim.csv exists
  --help, -h         Show this message
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxa_num|-t) TAXA_NUM="$2"; shift 2 ;;
    --gene_trees|-g) GENE_TREES="$2"; shift 2 ;;
    --replicate|-r) REPLICATE="$2"; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --simphy-dir) SIMPHY_DIR="$2"; SIMPHY_DIR_SET=true; shift 2 ;;
    --simphy-data-dir) SIMPHY_DATA_DIR="$2"; SIMPHY_DATA_DIR_SET=true; shift 2 ;;
    --sb) SB="$2"; shift 2 ;;
    --spmin) SPMIN="$2"; shift 2 ;;
    --spmax) SPMAX="$2"; shift 2 ;;
    --fresh) FRESH=true; shift 1 ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$TAXA_NUM" || -z "$GENE_TREES" ]]; then
  echo "Error: --taxa_num and --gene_trees are required."
  print_help
  exit 2
fi

# derive SIMPHY_DIR
if [[ "$SIMPHY_DIR_SET" = false ]]; then
  SIMPHY_DIR="${BASE_DIR%/}/STELAR-MP/simphy"
fi

# Construct expected output paths early
OUT_DIR="${SIMPHY_DIR%/}/data/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}"
REPL_DIR="${OUT_DIR%/}/${REPLICATE}"
CSV_FILE="${REPL_DIR%/}/stat-sim.csv"

# Checkpoint mechanism
if [[ "$FRESH" = false && -f "${CSV_FILE}" ]]; then
  echo "SKIPPING: ${CSV_FILE} already exists. Use --fresh to rerun."
  exit 0
fi

echo "Parameters:"
echo "  taxa_num:    $TAXA_NUM"
echo "  gene_trees:  $GENE_TREES"
echo "  replicate:   $REPLICATE"
echo "  base_dir:    $BASE_DIR"
echo "  simphy_dir:  $SIMPHY_DIR"
echo "  sb:          $SB"
echo "  spmin:       $SPMIN"
echo "  spmax:       $SPMAX"
echo

if [[ ! -d "$SIMPHY_DIR" ]]; then
  echo "Error: simphy directory \"$SIMPHY_DIR\" not found."
  exit 3
fi

# Path to run_simulator.sh inside simphy dir
RUN_SCRIPT="${SIMPHY_DIR%/}/run_simulator.sh"

if [[ ! -x "$RUN_SCRIPT" ]]; then
  echo "Error: $RUN_SCRIPT not found or not executable."
  exit 4
fi

# Call the main simulator script with named options (out_dir handled by run_simulator.sh)
echo "==> Running run_simulator.sh in $SIMPHY_DIR"

# Build command with optional data directory
RUN_CMD="./run_simulator.sh -t \"${TAXA_NUM}\" -g \"${GENE_TREES}\" --sb \"${SB}\" --spmin \"${SPMIN}\" --spmax \"${SPMAX}\""

# Add data directory option if specified
if [[ "$SIMPHY_DATA_DIR_SET" = true ]]; then
  # Ensure the data directory is created
  mkdir -p "$SIMPHY_DATA_DIR"
  RUN_CMD="$RUN_CMD --data_dir \"${SIMPHY_DATA_DIR}\""
  echo "Using custom simphy data directory: $SIMPHY_DATA_DIR"
fi

(
  cd "$SIMPHY_DIR"
  eval "$RUN_CMD"
)

# Construct expected out_dir the same way run_simulator.sh does
if [[ "$SIMPHY_DATA_DIR_SET" = true ]]; then
  # If custom data directory is used, construct path accordingly
  OUT_DIR="${SIMPHY_DATA_DIR%/}/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}"
else
  # Default behavior: data directory inside SIMPHY_DIR
  OUT_DIR="${SIMPHY_DIR%/}/data/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}"
fi
REPL_DIR="${OUT_DIR%/}/${REPLICATE}"
ALL_GT_FILE="${REPL_DIR%/}/all_gt.tre"

echo
if [[ -d "${REPL_DIR}" ]]; then
  echo "Simulator produced replicate directory: ${REPL_DIR}"
  if [[ -f "${ALL_GT_FILE}" ]]; then
    echo "Gene trees available at: ${ALL_GT_FILE}"
  else
    echo "Warning: expected gene-tree file not found at ${ALL_GT_FILE}"
  fi
else
  echo "Warning: expected replicate directory ${REPL_DIR} not found."
fi

# -------------------------
# New block: run analysis and write stat-sim.csv
# -------------------------
# Only proceed if replicate dir exists and all_gt.tre exists
if [[ -d "${REPL_DIR}" && -f "${ALL_GT_FILE}" ]]; then
  # Try to locate a species tree file within the replicate dir.
  SPECIES_TREE=""
  for f in "${REPL_DIR}"/s_tree* "${REPL_DIR}"/species_tree* "${REPL_DIR}"/species*; do
    if [[ -f "$f" ]]; then
      SPECIES_TREE="$f"
      break
    fi
  done

  if [[ -z "${SPECIES_TREE}" ]]; then
    echo "Warning: no species-tree file discovered in ${REPL_DIR}. Skipping analysis."
  else
    echo "Found species tree: ${SPECIES_TREE}"
    # Try to run analyze-dataset.py. We'll attempt from replicate dir first, then from SIMPHY_DIR, then rely on PATH.
    ANALYZE_OUTPUT=""
    set +e
    (
      cd "${REPL_DIR}"
      ANALYZE_OUTPUT=$(python analyze-dataset.py "${ALL_GT_FILE##*/}" "${SPECIES_TREE##*/}" 2>&1) || true
    )
    if [[ -z "${ANALYZE_OUTPUT}" ]]; then
      (
        cd "${SIMPHY_DIR}"
        ANALYZE_OUTPUT=$(python analyze-dataset.py "${ALL_GT_FILE}" "${SPECIES_TREE}" 2>&1) || true
      )
    fi
    if [[ -z "${ANALYZE_OUTPUT}" ]]; then
      ANALYZE_OUTPUT=$(python analyze-dataset.py "${ALL_GT_FILE}" "${SPECIES_TREE}" 2>&1) || true
    fi
    set -e

    if [[ -z "${ANALYZE_OUTPUT}" ]]; then
      echo "Warning: analyze-dataset.py produced no output or failed to run. Skipping CSV write."
    else
      echo "${ANALYZE_OUTPUT}"

      # Extract gt-gt and gt-st using patterns from your example output lines
      gt_gt=$(echo "${ANALYZE_OUTPUT}" | grep -F "Average pairwise discordance (gene trees):" | awk '{print $NF}' || true)
      gt_st=$(echo "${ANALYZE_OUTPUT}" | grep -F "Average discordance (gene trees vs species tree):" | awk '{print $NF}' || true)

      # If extraction failed, try slightly different phrases (be permissive)
      if [[ -z "${gt_gt}" ]]; then
        gt_gt=$(echo "${ANALYZE_OUTPUT}" | grep -E "Average pairwise.*discord|pairwise RF" | awk '{print $NF}' || true)
      fi
      if [[ -z "${gt_st}" ]]; then
        gt_st=$(echo "${ANALYZE_OUTPUT}" | grep -E "discordance .*gene trees vs species tree|gene trees vs species" | awk '{print $NF}' || true)
      fi

      echo ""
      echo "============================================================"
      echo "EXTRACTED VALUES:"
      echo "============================================================"
      echo "gt_gt: ${gt_gt:-N/A}"
      echo "gt_st: ${gt_st:-N/A}"

      # Prepare CSV path
      CSV_FILE="${REPL_DIR%/}/stat-sim.csv"

      # Write header if not present
      if [[ ! -f "${CSV_FILE}" ]]; then
        echo "num-taxa,gene-trees,replicate,sb,spmin,spmax,gt-gt,gt-st" > "${CSV_FILE}"
      fi

      # Append the row (use empty fields if values missing)
      printf "%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "${TAXA_NUM}" "${GENE_TREES}" "${REPLICATE}" "${SB}" "${SPMIN}" "${SPMAX}" \
        "${gt_gt:-}" "${gt_st:-}" >> "${CSV_FILE}"

      echo "Wrote/updated stats CSV: ${CSV_FILE}"
    fi
  fi
else
  echo "Skipping analysis: replicate directory or all_gt.tre missing."
fi

echo "Done. Output directory (if created): ${OUT_DIR}"

