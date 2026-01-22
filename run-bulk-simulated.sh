#!/usr/bin/env bash
# run-bulk-simulated.sh
#
# Runs sim.sh and test-stelar-simulated.sh or test-aster-simulated.sh 
# over all combinations of parameter lists.
#
# Usage:
#   ./run-bulk-simulated.sh -m stelar
#   ./run-bulk-simulated.sh -m aster --base-dir /path/to/research
#
# Default base-dir = $HOME/phylogeny

set -euo pipefail

BASE_DIR=""
BASE_DIR_PROVIDED=false
METHOD="stelar"  # default method: stelar or aster
FRESH=false
NUM_REPLICATES=1

print_help() {
  cat <<EOF
run-bulk-simulated.sh

Runs sim.sh and test-{stelar,aster}-simulated.sh for all combinations of parameter lists.

Options:
  --method, -m      Method to use: 'stelar' or 'aster' (default: stelar)
  --base-dir, -b    Base directory (optional, passed to sub-scripts if provided)
  --num-replicates, -n  Number of replicates to run (default: 1)
  --fresh           Pass --fresh to sim.sh and test scripts (recreate outputs)
  --help, -h        Show this message

Examples:
  ./run-bulk-simulated.sh -m stelar
  ./run-bulk-simulated.sh -m aster --base-dir /path/to/research
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --method|-m) METHOD="$2"; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; BASE_DIR_PROVIDED=true; shift 2 ;;
    --num-replicates|-n) NUM_REPLICATES="$2"; shift 2 ;;
    --fresh) FRESH=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Validate method
if [[ "$METHOD" != "stelar" && "$METHOD" != "aster" ]]; then
  echo "Error: --method must be 'stelar' or 'aster', got: $METHOD"
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

T_LIST=(10)
G_LIST=(100 200 1000 2500 5000)
SB_LIST=(0.000001)
SPMIN_LIST=(100000)
SPMAX_LIST=(200000)

# T_LIST=(30000 40000)
# G_LIST=(1000)
# SB_LIST=(0.000001)
# SPMIN_LIST=(100000)
# SPMAX_LIST=(150000)

# Number of replicates to run
# NUM_REPLICATES=5  # Now set via --num-replicates flag (default: 1)

# -------------------------------
# execution
# -------------------------------

# Build base-dir argument if provided
if $BASE_DIR_PROVIDED; then
  BASE_DIR_ARG="--base-dir $BASE_DIR"
  echo "Base dir: $BASE_DIR"
else
  BASE_DIR_ARG=""
  echo "Base dir: (not specified, scripts will use their defaults)"
fi

# Build fresh argument if provided
if $FRESH; then
  FRESH_ARG="--fresh"
  echo "Fresh:    yes"
else
  FRESH_ARG=""
  echo "Fresh:    no"
fi
echo "Method:   $METHOD"
echo "Replicates: $NUM_REPLICATES"
echo "Starting bulk runs..."

for t in "${T_LIST[@]}"; do
  for g in "${G_LIST[@]}"; do
    for sb in "${SB_LIST[@]}"; do
      for spmin in "${SPMIN_LIST[@]}"; do
        for spmax in "${SPMAX_LIST[@]}"; do

          echo ">>> Running: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax (method=$METHOD)"
          
          ./sim.sh -rs $NUM_REPLICATES $BASE_DIR_ARG -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" $FRESH_ARG
          
          # Run replicates
          for ((i=1; i<=NUM_REPLICATES; i++)); do
            echo "  Running replicate R$i with $METHOD"
            
            if [[ "$METHOD" == "stelar" ]]; then
              ./test-stelar-simulated.sh -r "R$i" $BASE_DIR_ARG -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" $FRESH_ARG
            else
              ./test-aster-simulated.sh -r "R$i" $BASE_DIR_ARG -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" --aster-opts="-t 32" $FRESH_ARG
            fi
          done

        done
      done
    done
  done
done

echo "All runs finished."
