#!/usr/bin/env bash
# run-bulk-simulated.sh
#
# Runs sim.sh and test-stelar-simulated.sh or test-baseline-simulated.sh
# over all combinations of parameter lists.
#
# Usage:
#   ./run-bulk-simulated.sh -m stelar
#   ./run-bulk-simulated.sh -m aster --base-dir /path/to/research
#   ./run-bulk-simulated.sh -m astral --base-dir /path/to/research
#
# Default base-dir = $HOME/phylogeny

set -euo pipefail

BASE_DIR=""
BASE_DIR_PROVIDED=false
METHOD="stelar"  # default method
FRESH=false
NUM_REPLICATES=1

# Method-specific options (passed through)
ASTER_OPTS=""
ASTER_BIN=""
ASTRAL_OPTS=""
ASTRAL_XMS=""
ASTRAL_XMX=""
TREEQMC_OPTS=""
WQFM_OPTS=""
SUPERTRIPLETS_OPTS=""
TMC_OPTS=""

print_help() {
  cat <<EOF
run-bulk-simulated.sh

Runs sim.sh and test-stelar-simulated.sh or test-baseline-simulated.sh for all combinations of parameter lists.

Options:
  --method, -m      Method to use: stelar, aster, astral, treeqmc, wqfmtree, supertriplets, stp-nni, tmc (default: stelar)
  --base-dir, -b    Base directory (optional, passed to sub-scripts if provided)
  --num-replicates, -n  Number of replicates to run (default: 1)
  --fresh           Pass --fresh to sim.sh and test scripts (recreate outputs)
  --aster-opts      Extra options for ASTER (for method=aster)
  --aster-bin       Path to ASTER binary (for method=aster)
  --astral-opts     Extra options for ASTRAL (for method=astral)
  --astral-xms      ASTRAL Java Xms (for method=astral)
  --astral-xmx      ASTRAL Java Xmx (for method=astral)
  --astral-Xms      Same as --astral-xms
  --astral-Xmx      Same as --astral-xmx
  --treeqmc-opts    Extra options for TreeQMC (for method=treeqmc)
  --wqfm-opts       Extra options for wQFM-TREE (for method=wqfmtree)
  --supertriplets-opts  Extra options for SuperTriplets (for method=supertriplets)
  --tmc-opts        Extra options for TMC (for method=tmc)
  --help, -h        Show this message

Examples:
  ./run-bulk-simulated.sh -m stelar
  ./run-bulk-simulated.sh -m aster --base-dir /path/to/research
  ./run-bulk-simulated.sh -m astral --base-dir /path/to/research
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --method|-m) METHOD="$2"; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; BASE_DIR_PROVIDED=true; shift 2 ;;
    --num-replicates|-n) NUM_REPLICATES="$2"; shift 2 ;;
    --aster-opts) ASTER_OPTS="$2"; shift 2 ;;
    --aster-opts=*) ASTER_OPTS="${1#*=}"; shift ;;
    --aster-bin) ASTER_BIN="$2"; shift 2 ;;
    --astral-opts) ASTRAL_OPTS="$2"; shift 2 ;;
    --astral-opts=*) ASTRAL_OPTS="${1#*=}"; shift ;;
    --astral-xms|--astral-Xms) ASTRAL_XMS="$2"; shift 2 ;;
    --astral-xmx|--astral-Xmx) ASTRAL_XMX="$2"; shift 2 ;;
    --treeqmc-opts) TREEQMC_OPTS="$2"; shift 2 ;;
    --treeqmc-opts=*) TREEQMC_OPTS="${1#*=}"; shift ;;
    --wqfm-opts) WQFM_OPTS="$2"; shift 2 ;;
    --wqfm-opts=*) WQFM_OPTS="${1#*=}"; shift ;;
    --supertriplets-opts) SUPERTRIPLETS_OPTS="$2"; shift 2 ;;
    --supertriplets-opts=*) SUPERTRIPLETS_OPTS="${1#*=}"; shift ;;
    --tmc-opts) TMC_OPTS="$2"; shift 2 ;;
    --tmc-opts=*) TMC_OPTS="${1#*=}"; shift ;;
    --fresh) FRESH=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Validate method
case "$METHOD" in
  stelar|aster|astral|treeqmc|tree-qmc|wqfmtree|wqfm-tree|supertriplets|stp-nni|tmc) ;;
  *)
  echo "Error: --method must be one of: stelar, aster, astral, treeqmc, wqfmtree, supertriplets, stp-nni, tmc. Got: $METHOD"
  exit 1
  ;;
esac

if [[ "$METHOD" == "tree-qmc" ]]; then
  METHOD="treeqmc"
fi
if [[ "$METHOD" == "wqfm-tree" ]]; then
  METHOD="wqfmtree"
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
G_LIST=(100 200 1000 2500 5000)
SB_LIST=(0.000001)
SPMIN_LIST=(100000)
SPMAX_LIST=(200000)

T_LIST=(7500)
G_LIST=(100 200 1000 2500 5000)
SB_LIST=(0.000001)
SPMIN_LIST=(100000)
SPMAX_LIST=(200000)

# T_LIST=(10)
# G_LIST=(10)
# SB_LIST=(0.000001)
# SPMIN_LIST=(100000)
# SPMAX_LIST=(200000)

T_LIST=(100 200 500)
G_LIST=(1000)
SB_LIST=(0.000001)
SPMIN_LIST=(50000)
SPMAX_LIST=(1000000)


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
          
          ./sim.sh -rs $NUM_REPLICATES $BASE_DIR_ARG -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" --fresh
          
          # Run replicates
          for ((i=1; i<=NUM_REPLICATES; i++)); do
            echo "  Running replicate R$i with $METHOD"
            
            if [[ "$METHOD" == "stelar" ]]; then
              ./test-stelar-simulated.sh -r "R$i" $BASE_DIR_ARG -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" $FRESH_ARG
            else
              METHOD_ARGS=(--method "$METHOD")
              if [[ -n "$ASTER_OPTS" ]]; then
                METHOD_ARGS+=(--aster-opts "$ASTER_OPTS")
              fi
              if [[ -n "$ASTER_BIN" ]]; then
                METHOD_ARGS+=(--aster-bin "$ASTER_BIN")
              fi
              if [[ -n "$ASTRAL_OPTS" ]]; then
                METHOD_ARGS+=(--astral-opts "$ASTRAL_OPTS")
              fi
              if [[ -n "$ASTRAL_XMS" ]]; then
                METHOD_ARGS+=(--astral-xms "$ASTRAL_XMS")
              fi
              if [[ -n "$ASTRAL_XMX" ]]; then
                METHOD_ARGS+=(--astral-xmx "$ASTRAL_XMX")
              fi
              if [[ -n "$TREEQMC_OPTS" ]]; then
                METHOD_ARGS+=(--treeqmc-opts "$TREEQMC_OPTS")
              fi
              if [[ -n "$WQFM_OPTS" ]]; then
                METHOD_ARGS+=(--wqfm-opts "$WQFM_OPTS")
              fi
              if [[ -n "$SUPERTRIPLETS_OPTS" ]]; then
                METHOD_ARGS+=(--supertriplets-opts "$SUPERTRIPLETS_OPTS")
              fi
              if [[ -n "$TMC_OPTS" ]]; then
                METHOD_ARGS+=(--tmc-opts "$TMC_OPTS")
              fi

              ./test-baseline-simulated.sh -r "R$i" $BASE_DIR_ARG -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${METHOD_ARGS[@]}" $FRESH_ARG
            fi
          done

        done
      done
    done
  done
done

echo "All runs finished."
