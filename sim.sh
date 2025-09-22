#!/usr/bin/env bash
# sim.sh (new runner for updated run_simulator.sh)
# Usage examples:
#   ./sim.sh -t 1000 -g 500
#   ./sim.sh -t 1000 -g 500 -r R1 -b "$HOME/phylogeny" --simphy-dir /opt/simphy --sb 0.000002 --spmin 400000 --spmax 1200000

set -euo pipefail

# Defaults
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
BASE_DIR="${HOME}/phylogeny"
SIMPHY_DIR=""          # will be derived from BASE_DIR unless explicitly provided
SIMPHY_DIR_SET=false

# Defaults that match your run_simulator.sh defaults
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
  --sb               Substitution/birthrate parameter (default: ${SB})
  --spmin            Population size minimum (default: ${SPMIN})
  --spmax            Population size maximum (default: ${SPMAX})
  --help, -h         Show this message

Behavior:
  If --simphy-dir is not provided, the script uses:
    \${BASE_DIR}/STELAR-MP/simphy

Examples:
  ./sim.sh -t 1000 -g 500
  ./sim.sh -t 1000 -g 500 --replicate R2 -b /home/user/research --sb 0.000002
  ./sim.sh -t 1000 -g 500 --simphy-dir /opt/simphy
EOF
}

# parse args (supports short and long)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxa_num|-t) TAXA_NUM="$2"; shift 2 ;;
    --gene_trees|-g) GENE_TREES="$2"; shift 2 ;;
    --replicate|-r) REPLICATE="$2"; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --simphy-dir) SIMPHY_DIR="$2"; SIMPHY_DIR_SET=true; shift 2 ;;
    --sb) SB="$2"; shift 2 ;;
    --spmin) SPMIN="$2"; shift 2 ;;
    --spmax) SPMAX="$2"; shift 2 ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$TAXA_NUM" || -z "$GENE_TREES" ]]; then
  echo "Error: --taxa_num and --gene_trees are required."
  print_help
  exit 2
fi

# derive SIMPHY_DIR from BASE_DIR if not explicitly set
if [[ "$SIMPHY_DIR_SET" = false ]]; then
  SIMPHY_DIR="${BASE_DIR%/}/STELAR-MP/simphy"
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
(
  cd "$SIMPHY_DIR"
  ./run_simulator.sh -t "${TAXA_NUM}" -g "${GENE_TREES}" \
    --sb "${SB}" --spmin "${SPMIN}" --spmax "${SPMAX}"
)

# Construct expected out_dir the same way run_simulator.sh does (relative to SIMPHY_DIR)
OUT_DIR="${SIMPHY_DIR%/}/data/t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}"
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

echo "Done. Output directory (if created): ${OUT_DIR}"
