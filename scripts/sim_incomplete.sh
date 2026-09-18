#!/usr/bin/env bash
# sim_incomplete.sh
#
# Wraps sim.sh: runs the full simulation (or skips if already done), then
# generates an incomplete-gene-tree variant by randomly pruning taxa.
#
# The incomplete trees are stored alongside the complete ones, in a directory
# whose name is the standard sim directory name with "_incomplete" appended:
#
#   complete:   $PHYLOGENY_DATA_DIR/simphy/data/t_N_g_K_sb_S_spmin_A_spmax_B/Ri/all_gt.tre
#   incomplete: $PHYLOGENY_DATA_DIR/simphy/data/t_N_g_K_sb_S_spmin_A_spmax_B_incomplete/Ri/all_gt.tre
#
# The true species tree (s_tree.trees) is copied into each incomplete replicate
# directory so that test-stelarx-simulated.sh --incomplete can compute RF distance.
#
# All sim.sh options are accepted and forwarded. Incomplete-specific options:
#   --fraction F    Fraction of taxa to remove per tree (default: 0.30)
#   --seed N        Random seed for reproducibility (default: 42)
#   --min-keep N    Minimum taxa to keep per tree (default: 4)
#   --fresh-inc     Regenerate incomplete trees even if they already exist
#
# Usage examples:
#   ./sim_incomplete.sh -t 1000 -g 500
#   ./sim_incomplete.sh -t 1000 -g 500 --fraction 0.5 --seed 7
#   ./sim_incomplete.sh -t 1000 -g 500 -rs 5 -r R2 --fraction 0.3 --fresh-inc

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/phylogeny-data-dir.sh"
PYTHON_BIN="${STELARX_PYTHON:-${SCRIPT_DIR}/.venv/bin/python}"
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN="python3"

# ── Incomplete-specific defaults ──────────────────────────────────────────────
FRACTION="0.30"
SEED="42"
MIN_KEEP="4"
FRESH_INC=false

# ── sim.sh param mirrors (needed to reconstruct output dir) ───────────────────
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
REPLICATES="10"
BASE_DIR="$SCRIPT_DIR"
SIMPHY_DIR=""
SIMPHY_DIR_SET=false
SIMPHY_DATA_DIR=""
SIMPHY_DATA_DIR_SET=false
SB="0.000001"
SPMIN="500000"
SPMAX="1500000"

# Collect all args; peel out incomplete-specific ones; forward the rest to sim.sh
SIM_ARGS=()

print_help() {
  cat <<EOF
sim_incomplete.sh — Run simulation and generate incomplete gene trees.

Accepts all sim.sh options, plus:

  --fraction F     Fraction of taxa to remove per tree (default: 0.30)
  --seed N         Random seed for gen_incomplete.py (default: 42)
  --min-keep N     Minimum taxa to keep per tree (default: 4)
  --fresh-inc      Re-generate incomplete trees even if they already exist

Incomplete trees are stored in:
  <sim_output_dir>_incomplete/<replicate>/all_gt.tre

Examples:
  ./sim_incomplete.sh -t 1000 -g 500
  ./sim_incomplete.sh -t 1000 -g 500 --fraction 0.5 --seed 7 --fresh-inc
  ./sim_incomplete.sh -t 1000 -g 500 -rs 5 --fraction 0.3
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    # Incomplete-specific — consume here, do NOT forward to sim.sh
    --fraction)    FRACTION="$2";  shift 2 ;;
    --seed)        SEED="$2";      shift 2 ;;
    --min-keep)    MIN_KEEP="$2";  shift 2 ;;
    --fresh-inc)   FRESH_INC=true; shift ;;

    # sim.sh args we also need locally — mirror AND forward
    --taxa_num|-t)       TAXA_NUM="$2";  SIM_ARGS+=("$1" "$2"); shift 2 ;;
    --gene_trees|-g)     GENE_TREES="$2"; SIM_ARGS+=("$1" "$2"); shift 2 ;;
    --replicate|-r)      REPLICATE="$2"; SIM_ARGS+=("$1" "$2"); shift 2 ;;
    --replicates|-rs)    REPLICATES="$2"; SIM_ARGS+=("$1" "$2"); shift 2 ;;
    --project-root|--base-dir|-b)
                          BASE_DIR="$2"; SIM_ARGS+=("--project-root" "$2"); shift 2 ;;
    --simphy-dir)        SIMPHY_DIR="$2"; SIMPHY_DIR_SET=true; SIM_ARGS+=("$1" "$2"); shift 2 ;;
    --simphy-data-dir)   SIMPHY_DATA_DIR="$2"; SIMPHY_DATA_DIR_SET=true; shift 2 ;;
    --sb)                SB="$2";    SIM_ARGS+=("$1" "$2"); shift 2 ;;
    --spmin)             SPMIN="$2"; SIM_ARGS+=("$1" "$2"); shift 2 ;;
    --spmax)             SPMAX="$2"; SIM_ARGS+=("$1" "$2"); shift 2 ;;

    # Forward all other sim.sh args unchanged
    --fresh)             SIM_ARGS+=("$1"); shift ;;
    --help|-h)           print_help; exit 0 ;;
    *)                   SIM_ARGS+=("$1"); shift ;;
  esac
done

if [[ -z "$TAXA_NUM" || -z "$GENE_TREES" ]]; then
  echo "Error: --taxa_num (-t) and --gene_trees (-g) are required."
  print_help
  exit 1
fi

# ── Derive SIMPHY_DIR (same logic as sim.sh) ──────────────────────────────────
if [[ "$SIMPHY_DIR_SET" == false ]]; then
  SIMPHY_DIR="${BASE_DIR%/}/simphy"
fi
SIMPHY_DIR="$(realpath "$SIMPHY_DIR")"
SIMPHY_DATA_DIR="$(stelarx_prepare_simphy_data_dir "$SIMPHY_DATA_DIR")"
SIM_ARGS+=(--simphy-data-dir "$SIMPHY_DATA_DIR")

# ── Derive complete output dir (same logic as sim.sh) ────────────────────────
DATASET_NAME="t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}"
COMPLETE_DIR="${SIMPHY_DATA_DIR%/}/${DATASET_NAME}"

INCOMPLETE_DIR="${COMPLETE_DIR}_incomplete"

# ── Step 1: Run sim.sh ────────────────────────────────────────────────────────
echo "=== Step 1: Running sim.sh ==="
bash "${SCRIPT_DIR}/sim.sh" "${SIM_ARGS[@]}"
echo

# ── Step 2: Generate incomplete trees for each replicate ─────────────────────
echo "=== Step 2: Generating incomplete gene trees (fraction=${FRACTION}, seed=${SEED}) ==="
echo "  Complete dir:   ${COMPLETE_DIR}"
echo "  Incomplete dir: ${INCOMPLETE_DIR}"
echo

GEN_SCRIPT="${SCRIPT_DIR}/test/gen_incomplete.py"
if [[ ! -f "$GEN_SCRIPT" ]]; then
  echo "Error: gen_incomplete.py not found at ${GEN_SCRIPT}"
  exit 1
fi

for i in $(seq 1 "${REPLICATES}"); do
  REPL="R${i}"
  SRC_FILE="${COMPLETE_DIR}/${REPL}/all_gt.tre"
  DST_DIR="${INCOMPLETE_DIR}/${REPL}"
  DST_FILE="${DST_DIR}/all_gt.tre"
  TRUE_TREE_SRC="${COMPLETE_DIR}/${REPL}/s_tree.trees"
  TRUE_TREE_DST="${DST_DIR}/s_tree.trees"

  if [[ ! -f "$SRC_FILE" ]]; then
    echo "  [R${i}] SKIP — source not found: ${SRC_FILE}"
    continue
  fi

  if [[ "$FRESH_INC" == false && -f "$DST_FILE" ]]; then
    echo "  [R${i}] SKIP — incomplete trees already exist (use --fresh-inc to regenerate)"
    continue
  fi

  mkdir -p "$DST_DIR"

  "$PYTHON_BIN" "$GEN_SCRIPT" \
    "$SRC_FILE" "$DST_FILE" \
    --fraction "$FRACTION" \
    --seed "$SEED" \
    --min-keep "$MIN_KEEP" \
    --stats 2>&1 | sed "s/^/  [R${i}] /"

  # Copy true species tree so test-stelarx-simulated.sh --incomplete can compute RF
  if [[ -f "$TRUE_TREE_SRC" ]]; then
    cp "$TRUE_TREE_SRC" "$TRUE_TREE_DST"
  fi

  echo "  [R${i}] Done → ${DST_FILE}"
done

# Record how the incomplete variant was derived, next to SimPhy's own
# <dataset>.command in the complete directory, so the outputs mirror and any
# upload carry the full reproducibility fingerprint of this dataset.
if [[ -d "$INCOMPLETE_DIR" ]]; then
  INCOMPLETE_COMMAND_FILE="${INCOMPLETE_DIR}/${DATASET_NAME}_incomplete.command"
  INCOMPLETE_COMMAND_TMP="${INCOMPLETE_COMMAND_FILE}.tmp.$$"
  {
    printf '# Incomplete gene-tree variant derived from %s\n' "$DATASET_NAME"
    printf '# Base SimPhy command: %s\n' "${COMPLETE_DIR}/${DATASET_NAME}.command"
    printf '%q -t %q -g %q -rs %q --sb %q --spmin %q --spmax %q --fraction %q --seed %q --min-keep %q\n' \
      "${SCRIPT_DIR}/sim_incomplete.sh" "$TAXA_NUM" "$GENE_TREES" "$REPLICATES" \
      "$SB" "$SPMIN" "$SPMAX" "$FRACTION" "$SEED" "$MIN_KEEP"
    printf '# Per replicate: %q <complete>/<R>/all_gt.tre <incomplete>/<R>/all_gt.tre --fraction %q --seed %q --min-keep %q --stats\n' \
      "$GEN_SCRIPT" "$FRACTION" "$SEED" "$MIN_KEEP"
  } > "$INCOMPLETE_COMMAND_TMP"
  mv -f "$INCOMPLETE_COMMAND_TMP" "$INCOMPLETE_COMMAND_FILE"
fi

echo
echo "=== Incomplete simulation complete ==="
echo "  Dataset:    ${DATASET_NAME}"
echo "  Directory:  ${INCOMPLETE_DIR}"
echo "  Fraction:   ${FRACTION}  Seed: ${SEED}"
