#!/usr/bin/env bash
# run-bulk-simulated.sh
#
# Runs sim.sh and test-stelar.sh over all combinations of parameter lists.
#
# Usage:
#   ./run-bulk-simulated.sh
#   ./run-bulk-simulated.sh --base-dir /path/to/research
#
# Default base-dir = $HOME/phylogeny

set -euo pipefail

BASE_DIR="${HOME}/phylogeny"

print_help() {
  cat <<EOF
bulk-run.sh

Runs sim.sh and test-stelar.sh for all combinations of parameter lists.

Options:
  --base-dir, -b    Base directory (default: ${BASE_DIR})
  --help, -h        Show this message
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# -------------------------------
# parameter lists (EDIT AS NEEDED)
# -------------------------------
# T_LIST=(1000 2000 5000 10000 15000 20000 25000 30000)
# G_LIST=(1000)
# SB_LIST=(0.000001)
# SPMIN_LIST=(50000 100000)
# SPMAX_LIST=(150000 200000 250000 300000)

# T_LIST=(1000)
# G_LIST=(100 200)
# SB_LIST=(0.000001)
# SPMIN_LIST=(50000)
# SPMAX_LIST=(150000)

T_LIST=(1000)
G_LIST=(100 1000 2500 5000)
SB_LIST=(0.000001)
SPMIN_LIST=(100000)
SPMAX_LIST=(200000)

# T_LIST=(30000 40000)
# G_LIST=(1000)
# SB_LIST=(0.000001)
# SPMIN_LIST=(100000)
# SPMAX_LIST=(150000)

# Number of replicates to run
NUM_REPLICATES=5

# -------------------------------
# execution
# -------------------------------
echo "Base dir: $BASE_DIR"
echo "Starting bulk runs..."

for t in "${T_LIST[@]}"; do
  for g in "${G_LIST[@]}"; do
    for sb in "${SB_LIST[@]}"; do
      for spmin in "${SPMIN_LIST[@]}"; do
        for spmax in "${SPMAX_LIST[@]}"; do

          echo ">>> Running: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax"
          
          ./sim.sh -rs $NUM_REPLICATES --base-dir "$BASE_DIR" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" --fresh
          
          # Run replicates
          for ((i=1; i<=NUM_REPLICATES; i++)); do
            echo "  Running replicate R$i"
            ./test-astral-with-mem.sh -r "R$i" --base-dir "$BASE_DIR" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" --fresh
          done

        done
      done
    done
  done
done

echo "All runs finished."
