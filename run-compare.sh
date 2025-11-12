#!/usr/bin/env bash
# run-compare.sh
#
# Runs sim.sh and multiple test scripts over all combinations of parameter lists and replicates.
#
# Usage:
#   ./run-compare.sh --replicates 5
#   ./run-compare.sh --base-dir /path/to/research --replicates 3
#
# Default base-dir = $HOME/phylogeny

set -euo pipefail

BASE_DIR="${HOME}/phylogeny"
REPLICATES=""
NO_NOTIFY=false

print_help() {
  cat <<EOF
run-compare.sh

Runs sim.sh and multiple test scripts for all combinations of parameter lists and replicates.

Options:
  --base-dir, -b      Base directory (default: ${BASE_DIR})
  --replicates, -rs   Number of replicates (required)
  --no-notify, -nn    Disable ntfy.sh notifications in all test scripts
  --help, -h          Show this message
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --replicates|-rs) REPLICATES="$2"; shift 2 ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Validate required arguments
if [[ -z "$REPLICATES" ]]; then
  echo "Error: --replicates is required"
  print_help
  exit 1
fi

if ! [[ "$REPLICATES" =~ ^[0-9]+$ ]] || [[ "$REPLICATES" -le 0 ]]; then
  echo "Error: --replicates must be a positive integer"
  exit 1
fi

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

T_LIST=(100)
G_LIST=(1000)
SB_LIST=(0.000001)
SPMIN_LIST=(50000)
SPMAX_LIST=(1000000)

# T_LIST=(30000 40000)
# G_LIST=(1000)
# SB_LIST=(0.000001)
# SPMIN_LIST=(100000)
# SPMAX_LIST=(150000)

# -------------------------------
# execution
# -------------------------------
echo "Base dir: $BASE_DIR"
echo "Replicates: $REPLICATES"
echo "Starting bulk runs..."

for t in "${T_LIST[@]}"; do
  for g in "${G_LIST[@]}"; do
    for sb in "${SB_LIST[@]}"; do
      for spmin in "${SPMIN_LIST[@]}"; do
        for spmax in "${SPMAX_LIST[@]}"; do
          echo ">>> Running: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax"
          
          # Run simulation with specified number of replicates
          ./sim.sh -b "$BASE_DIR" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" -rs "$REPLICATES" --fresh
          
          # Loop over each replicate and run all test scripts
          for ((r=1; r<=REPLICATES; r++)); do
            replicate="R$r"
            echo "  >>> Running replicate: $replicate"
            
            # Prepare notification flag
            NOTIFY_FLAG=""
            if [[ "$NO_NOTIFY" = true ]]; then
              NOTIFY_FLAG="--no-notify"
            fi
            
            # Run all test scripts for this replicate
            ./test-stelar-with-mem.sh -b "$BASE_DIR" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" -r "$replicate" --fresh $NOTIFY_FLAG
            ./test-astral-with-mem.sh -b "$BASE_DIR" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" -r "$replicate" --fresh $NOTIFY_FLAG
            ./test-tree-qmc.sh -b "$BASE_DIR" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" -r "$replicate" --fresh $NOTIFY_FLAG
            ./test-wqfm-tree.sh -b "$BASE_DIR" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" -r "$replicate" --fresh $NOTIFY_FLAG
          done
        done
      done
    done
  done
done

echo "All runs finished."
