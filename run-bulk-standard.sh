#!/usr/bin/env bash
# Multi-Algorithm Dataset Runner Script (updated)
# Supports STELAR, ASTRAL, TreeQMC, and wQFMtree algorithms
# Usage: ./run-bulk-standard.sh [--base-dir /path/to/base] [--dataset-dir /path/to/datasets] [--fresh]
#   --base-dir, -b    Optional base directory (defaults to value below)
#   --dataset-dir, -d Optional dataset directory (defaults to BASE_DIR/datasets)
#   --fresh           Force re-run even if stat-<alg>.csv exists

set -euo pipefail

# =============================================================================
# DEFAULTS (edit these if you want different defaults)
# =============================================================================
BASE_DIR="/home/aaniksahaa/research"  # default; can be overridden with --base-dir or -b
DATASET_DIR=""                        # dataset directory; will be set to BASE_DIR/datasets if not specified
STELAR_ROOT=""                        # derived from BASE_DIR if not set explicitly
ASTRAL_ROOT=""                        # derived from BASE_DIR if not set explicitly
TREEQMC_ROOT=""                       # derived from BASE_DIR if not set explicitly
WQFMTREE_ROOT=""                      # derived from BASE_DIR if not set explicitly
FRESH=false

# Algorithm configuration
ALGORITHMS=("stelar" "astral" "treeqmc" "wqfmtree")

# Algorithm-specific options
STELAR_OPTS="GPU_PARALLEL NONE"
ASTRAL_OPTS=""  # ASTRAL doesn't need special options for basic runs
TREEQMC_OPTS=""  # TreeQMC doesn't need special options for basic runs
WQFMTREE_OPTS=""  # wQFMtree doesn't need special options for basic runs

# Algorithm command mappings - these will be used to construct the actual commands
declare -A ALG_COMMANDS
ALG_COMMANDS["stelar"]="stelar_command"
ALG_COMMANDS["astral"]="astral_command"
ALG_COMMANDS["treeqmc"]="treeqmc_command"
ALG_COMMANDS["wqfmtree"]="wqfmtree_command"

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# =============================================================================
# Dataset configuration (kept exactly as you provided)
# =============================================================================
# folders=("200-taxon")
folders=("200-taxon")

declare -A innerFolderNames
innerFolderNames["11-taxon"]="estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS"
innerFolderNames["15-taxon"]="100gene-100bp/estimated-genetrees 100gene-1000bp/estimated-genetrees 100gene-true 1000gene-100bp/estimated-genetrees 1000gene-1000bp/estimated-genetrees 1000gene-true"
innerFolderNames["37-taxon"]="estimated-genetrees/1X-200-500 estimated-genetrees/1X-200-1000 estimated-genetrees/1X-400-500 estimated-genetrees/1X-400-1000 estimated-genetrees/1X-800-500 estimated-genetrees/1X-800-1000 estimated-genetrees/2X-200-500"
innerFolderNames["48-taxon"]="estimated-genetrees/1X-25-500 estimated-genetrees/1X-50-500 estimated-genetrees/1X-100-500 estimated-genetrees/1X-200-500 estimated-genetrees/1X-500-500 estimated-genetrees/1X-1000-500 estimated-genetrees/2X-1000-500"
innerFolderNames["48-taxon-latest"]="estimated_genetrees/1X-1000-500"
innerFolderNames["100-taxon"]="inner100"
innerFolderNames["200-taxon"]="inner200"
innerFolderNames["500-taxon"]="model.500.2000000.0.000001"
innerFolderNames["biological"]="nuclear"
innerFolderNames["mammalian"]="424genes"
innerFolderNames["amniota"]="aa nt"

declare -A replicates
replicates["11-taxon"]=20
replicates["15-taxon"]=10
replicates["37-taxon"]=20
replicates["48-taxon"]=20
replicates["48-taxon-latest"]=20
replicates["100-taxon"]=10
replicates["200-taxon"]=10
replicates["500-taxon"]=1
replicates["biological"]=1
replicates["mammalian"]=1
replicates["amniota"]=1

# =============================================================================
# Helper functions
# =============================================================================

print_header() {
    echo "==============================================="
    echo "=== Multi-Algorithm Dataset Runner (updated) ==="
    echo "==============================================="
    echo "BASE_DIR: $BASE_DIR"
    echo "DATASET_DIR: $DATASET_DIR"
    echo "STELAR_ROOT: ${STELAR_ROOT:-(derived from BASE_DIR)}"
    echo "ASTRAL_ROOT: ${ASTRAL_ROOT:-(derived from BASE_DIR)}"
    echo "TREEQMC_ROOT: ${TREEQMC_ROOT:-(derived from BASE_DIR)}"
    echo "WQFMTREE_ROOT: ${WQFMTREE_ROOT:-(derived from BASE_DIR)}"
    echo "Algorithms: ${ALGORITHMS[*]}"
    echo "STELAR options: $STELAR_OPTS"
    echo "Fresh run: $FRESH"
    echo "Folders to process: ${folders[*]}"
    echo "==============================================="
    echo
}

validate_base_dir() {
    if [ ! -d "$BASE_DIR" ]; then
        echo -e "${RED}Error: BASE_DIR '$BASE_DIR' does not exist.${NC}"
        exit 1
    fi

    if [ ! -d "$DATASET_DIR" ]; then
        echo -e "${RED}Error: DATASET_DIR '$DATASET_DIR' does not exist.${NC}"
        exit 1
    fi

    if [ ! -d "$BASE_DIR/RF" ]; then
        # RF scripts may exist elsewhere; just warn (we have fallbacks)
        echo -e "${YELLOW}Warning: RF folder not found in '$BASE_DIR'. RF calculations may fail unless rf.py in STELAR root or other script present.${NC}"
    fi
}

validate_algorithm_binaries() {
    local errors_found=false
    
    # Check if stelar is in algorithms array
    if [[ " ${ALGORITHMS[*]} " =~ " stelar " ]]; then
        # If STELAR_ROOT/run.sh exists we'll use it; otherwise we expect target/stelar-mp-1.0-SNAPSHOT.jar in cwd of STELAR_ROOT
        if [ ! -d "$STELAR_ROOT" ]; then
            echo -e "${RED}Error: STELAR_ROOT '$STELAR_ROOT' does not exist. Please ensure STELAR-MP is at ${STELAR_ROOT} or set STELAR_ROOT variable inside the script.${NC}"
            errors_found=true
        elif [ ! -x "${STELAR_ROOT}/run.sh" ] && [ ! -f "${STELAR_ROOT}/target/stelar-mp-1.0-SNAPSHOT.jar" ]; then
            echo -e "${RED}Error: STELAR binaries not found in STELAR_ROOT ('run.sh' or 'target/stelar-mp-1.0-SNAPSHOT.jar' expected). Run build.sh or set STELAR_ROOT properly.${NC}"
            errors_found=true
        fi
    fi
    
    # Check if astral is in algorithms array
    if [[ " ${ALGORITHMS[*]} " =~ " astral " ]]; then
        if [ ! -d "$ASTRAL_ROOT" ]; then
            echo -e "${RED}Error: ASTRAL_ROOT '$ASTRAL_ROOT' does not exist. Please ensure ASTRAL is at ${ASTRAL_ROOT} or set ASTRAL_ROOT variable inside the script.${NC}"
            errors_found=true
        elif [ ! -x "${ASTRAL_ROOT}/run_astral.sh" ]; then
            echo -e "${RED}Error: ASTRAL binary not found in ASTRAL_ROOT ('run_astral.sh' expected). Please ensure ASTRAL is properly installed.${NC}"
            errors_found=true
        fi
    fi
    
    # Check if treeqmc is in algorithms array
    if [[ " ${ALGORITHMS[*]} " =~ " treeqmc " ]]; then
        if [ ! -d "$TREEQMC_ROOT" ]; then
            echo -e "${RED}Error: TREEQMC_ROOT '$TREEQMC_ROOT' does not exist. Please ensure TREE-QMC is at ${TREEQMC_ROOT} or set TREEQMC_ROOT variable inside the script.${NC}"
            errors_found=true
        elif [ ! -x "${TREEQMC_ROOT}/tree-qmc" ]; then
            echo -e "${RED}Error: TreeQMC binary not found in TREEQMC_ROOT ('tree-qmc' expected). Please ensure TreeQMC is properly installed.${NC}"
            errors_found=true
        fi
    fi
    
    # Check if wqfmtree is in algorithms array
    if [[ " ${ALGORITHMS[*]} " =~ " wqfmtree " ]]; then
        if [ ! -d "$WQFMTREE_ROOT" ]; then
            echo -e "${RED}Error: WQFMTREE_ROOT '$WQFMTREE_ROOT' does not exist. Please ensure wQFM-TREE is at ${WQFMTREE_ROOT} or set WQFMTREE_ROOT variable inside the script.${NC}"
            errors_found=true
        elif [ ! -x "${WQFMTREE_ROOT}/run.sh" ]; then
            echo -e "${RED}Error: wQFMtree binary not found in WQFMTREE_ROOT ('run.sh' expected). Please ensure wQFMtree is properly installed.${NC}"
            errors_found=true
        fi
    fi
    
    if [ "$errors_found" = true ]; then
        exit 1
    fi
}

# Run algorithm for one input -> output, capture resource usage and produce stat file
run_algorithm_and_write_stats() {
    local ALL_GT_FILE="$1"
    local TRUE_SPECIES_TREE="$2"
    local OUT_DIR="$3"           # directory where output-alg.tre and stat-alg.csv will be written
    local REPLICATE_NAME="$4"    # e.g., R1
    local FOLDER_NAME="$5"
    local INNER_FOLDER_NAME="$6"
    local ALGORITHM="$7"         # algorithm name (stelar, astral, etc.)

    mkdir -p "$OUT_DIR"

    local OUT_FILE="${OUT_DIR%/}/output-${ALGORITHM}.tre"
    local STAT_FILE="${OUT_DIR%/}/stat-${ALGORITHM}.csv"

    # If stat file exists and not FRESH -> skip
    if [[ "$FRESH" = false && -f "$STAT_FILE" ]]; then
        echo "      SKIPPING: ${STAT_FILE} already exists. Use --fresh to force rerun."
        return 0
    fi

    echo "      Running ${ALGORITHM^^} (output -> $OUT_FILE)"
    # prepare tempfiles
    local TIME_TMP
    local MON_TMP
    TIME_TMP=$(mktemp)
    MON_TMP=$(mktemp)

    local START_NS
    START_NS=$(date +%s%N)

    # run algorithm under /usr/bin/time -v and capture stderr to TIME_TMP
    (
      case "$ALGORITHM" in
        "stelar")
          cd "$STELAR_ROOT"
          if [ -x "./run.sh" ]; then
            /usr/bin/time -v ./run.sh "$ALL_GT_FILE" "$OUT_FILE" $STELAR_OPTS
          else
            # fallback to java invocation (Main class) if run.sh missing
            /usr/bin/time -v java -Djava.library.path="$(pwd)/cuda" -Djna.debug_load=false -Djna.platform.library.path="$(pwd)/cuda" -cp target/stelar-mp-1.0-SNAPSHOT.jar Main -i "$ALL_GT_FILE" -o "$OUT_FILE" -m "$STELAR_OPTS"
          fi
          ;;
        "astral")
          cd "$ASTRAL_ROOT"
          /usr/bin/time -v ./run_astral.sh -i "$ALL_GT_FILE" -o "$OUT_FILE" $ASTRAL_OPTS
          ;;
        "treeqmc")
          cd "$TREEQMC_ROOT"
          # Ensure output directory exists
          mkdir -p "$(dirname "$OUT_FILE")"
          /usr/bin/time -v ./tree-qmc -i "$ALL_GT_FILE" -o "$OUT_FILE" $TREEQMC_OPTS
          ;;
        "wqfmtree")
          cd "$WQFMTREE_ROOT"
          # Ensure output directory exists
          mkdir -p "$(dirname "$OUT_FILE")"
          /usr/bin/time -v ./run.sh "$ALL_GT_FILE" "$OUT_FILE" $WQFMTREE_OPTS
          ;;
        *)
          echo "Unknown algorithm: $ALGORITHM"
          exit 1
          ;;
      esac
    ) 2> "$TIME_TMP" &
    local ALGORITHM_WRAPPER_PID=$!

    # GPU monitor (if available)
    local MON_PID=""
    if command -v nvidia-smi >/dev/null 2>&1; then
      (
        curmax=0
        # sample at 0.1s to catch short spikes
        while kill -0 "$ALGORITHM_WRAPPER_PID" 2>/dev/null; do
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
    fi

    # wait for algorithm and capture exit code
    wait "$ALGORITHM_WRAPPER_PID"
    local ALGORITHM_EXIT_CODE=$?

    local END_NS
    END_NS=$(date +%s%N)
    local ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
    local RUNNING_TIME
    RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

    # ensure monitor stopped and read its result
    if [[ -n "${MON_PID}" ]]; then
      wait "$MON_PID" 2>/dev/null || true
    fi
    local MAX_GPU_VAL
    MAX_GPU_VAL=$(cat "$MON_TMP" 2>/dev/null || echo "NA")

    # convert GPU MiB -> MB decimal
    local MAX_GPU_MB
    if [[ "$MAX_GPU_VAL" =~ ^[0-9]+$ ]]; then
      MAX_GPU_MB=$(awk "BEGIN{printf \"%.3f\", ${MAX_GPU_VAL} * 1.024}")
    else
      MAX_GPU_MB="NA"
    fi

    # parse /usr/bin/time -v output to get Maximum resident set size (kbytes)
    local MAX_CPU_MB="NA"
    if grep -qi "Maximum resident set size" "$TIME_TMP" 2>/dev/null; then
      MAX_RSS_KB=$(grep -i "Maximum resident set size" "$TIME_TMP" | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | awk '{print int($1)}')
    elif grep -qi "Maximum resident set size (kbytes)" "$TIME_TMP" 2>/dev/null; then
      MAX_RSS_KB=$(grep -i "Maximum resident set size (kbytes)" "$TIME_TMP" | awk -F: '{gsub(/^[ \t]+/,"",$2); print $2}' | awk '{print int($1)}')
    else
      MAX_RSS_KB=""
    fi

    if [[ -n "${MAX_RSS_KB}" && "${MAX_RSS_KB}" =~ ^[0-9]+$ ]]; then
      MAX_CPU_MB=$(awk "BEGIN{printf \"%.3f\", ${MAX_RSS_KB}/1024}")
    fi

    # clean up tempfiles
    rm -f "$TIME_TMP" "$MON_TMP" 2>/dev/null || true

    # RF calculation
    local RF_RATE="NA"
    if [[ -f "$OUT_FILE" ]]; then
      # Prefer rf.py inside STELAR_ROOT if present
      if [[ -f "${STELAR_ROOT%/}/rf.py" && -x "$(command -v python3)" ]]; then
        echo "      Calculating RF using ${STELAR_ROOT%/}/rf.py"
        rf_output=$(python3 "${STELAR_ROOT%/}/rf.py" "$OUT_FILE" "$TRUE_SPECIES_TREE" 2>&1) || rf_output="$rf_output"
        # try to find a number in output
        rf_candidate=$(echo "$rf_output" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
        if [[ -n "$rf_candidate" ]]; then
          RF_RATE="$rf_candidate"
        else
          RF_RATE="NA"
        fi
      elif [[ -f "${BASE_DIR%/}/RF/getFpFn.py" && -x "$(command -v python3)" ]]; then
        echo "      Calculating RF using ${BASE_DIR%/}/RF/getFpFn.py"
        # getFpFn.py returns a tuple; attempt to extract the same metric as earlier script
        tuple=$(python3 "${BASE_DIR%/}/RF/getFpFn.py" -e "$OUT_FILE" -t "$TRUE_SPECIES_TREE" 2>/dev/null) || tuple=""
        if [[ -n "$tuple" ]]; then
          # try to extract first numeric token
          rf_candidate=$(echo "$tuple" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
          if [[ -n "$rf_candidate" ]]; then
            RF_RATE="$rf_candidate"
          else
            RF_RATE="NA"
          fi
        else
          RF_RATE="NA"
        fi
      else
        echo "      RF calculation skipped: no rf.py or getFpFn.py available."
        RF_RATE="NA"
      fi
    else
      echo "      ${ALGORITHM^^} output missing; skipping RF calculation."
      RF_RATE="NA"
    fi

    # Write stat CSV (overwrite per run)
    # Header: alg,folder,inner_folder,replicate,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb,exit-code
    {
      echo "alg,folder,inner_folder,replicate,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb,exit-code"
      echo "${ALGORITHM},${FOLDER_NAME},${INNER_FOLDER_NAME},${REPLICATE_NAME},${RF_RATE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${ALGORITHM_EXIT_CODE}"
    } > "$STAT_FILE"

    echo -e "      ${GREEN}Wrote stats to $STAT_FILE${NC}"
    return 0
}

# =============================================================================
# MAIN processing loop
# =============================================================================

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-dir|-b)
      BASE_DIR="$2"
      shift 2
      ;;
    --dataset-dir|-d)
      DATASET_DIR="$2"
      shift 2
      ;;
    --stelar-root)
      STELAR_ROOT="$2"
      shift 2
      ;;
    --fresh)
      FRESH=true
      shift
      ;;
    --help|-h)
      cat <<EOF
Usage: $0 [--base-dir /path/to/base] [--dataset-dir /path/to/datasets] [--stelar-root /path/to/stelar] [--fresh]

Multi-algorithm dataset runner supporting STELAR, ASTRAL, TreeQMC, and wQFMtree.

--base-dir, -b     Base directory containing RF/ (overrides default)
--dataset-dir, -d  Dataset directory (overrides default BASE_DIR/datasets)
--stelar-root      STELAR root directory (overrides default derived from base dir)
--fresh            Force rerun even if stat-<alg>.csv exists
--help, -h         Show this help

Algorithms run: ${ALGORITHMS[*]}
Algorithm root directories (auto-derived from base dir):
  STELAR: \${BASE_DIR}/STELAR-MP
  ASTRAL: \${BASE_DIR}/ASTRAL
  TreeQMC: \${BASE_DIR}/TREE-QMC
  wQFMtree: \${BASE_DIR}/wQFM-TREE/wQFM-TREE
EOF
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# derive DATASET_DIR from BASE_DIR if not set
if [[ -z "${DATASET_DIR}" ]]; then
  DATASET_DIR="${BASE_DIR%/}/phylo-datasets"
fi

# derive STELAR_ROOT from BASE_DIR if not set
if [[ -z "${STELAR_ROOT}" ]]; then
  STELAR_ROOT="${BASE_DIR%/}/STELAR-MP"
fi

# derive ASTRAL_ROOT from BASE_DIR if not set
if [[ -z "${ASTRAL_ROOT}" ]]; then
  ASTRAL_ROOT="${BASE_DIR%/}/ASTRAL"
fi

# derive TREEQMC_ROOT from BASE_DIR if not set
if [[ -z "${TREEQMC_ROOT}" ]]; then
  TREEQMC_ROOT="${BASE_DIR%/}/TREE-QMC"
fi

# derive WQFMTREE_ROOT from BASE_DIR if not set
if [[ -z "${WQFMTREE_ROOT}" ]]; then
  WQFMTREE_ROOT="${BASE_DIR%/}/wQFM-TREE/wQFM-TREE"
fi

print_header

echo -e "${YELLOW}Performing validation checks...${NC}"
validate_base_dir
validate_algorithm_binaries
echo -e "${GREEN}✓ Validation checks passed${NC}"
echo

echo -e "${YELLOW}Starting dataset processing...${NC}"

# Loop through each taxa folder
for folder in "${folders[@]}"; do
    echo -e "${YELLOW}Processing folder: $folder${NC}"

    if [ ! -d "${DATASET_DIR%/}/$folder" ]; then
        echo -e "${RED}Warning: Folder '${DATASET_DIR%/}/$folder' does not exist, skipping...${NC}"
        continue
    fi

    IFS=' ' read -r -a inner_folders <<< "${innerFolderNames[$folder]}"
    local_R=${replicates[$folder]:-1}

    for inner_folder in "${inner_folders[@]}"; do
        echo "  Processing inner folder: $inner_folder"

        for (( j=1; j<=local_R; j++ )); do
            # For 100-taxon, skip R1 and R3
            if [[ "$folder" == "100-taxon" && ($j -eq 1 || $j -eq 3 || $j -eq 8 || $j -eq 9) ]]; then
                continue
            fi
            
            REPL="R${j}"
            echo "    Processing replicate ${REPL}"

            GT_FOLDER="${inner_folder}/${REPL}"
            ALL_GT_FILE="${DATASET_DIR%/}/$folder/$GT_FOLDER/all_gt.tre.rooted"
            TRUE_TREE="${DATASET_DIR%/}/$folder/true_tree_trimmed"

            # special handling for 100/200-taxon true tree path (kept from your original)
            if [[ "$folder" == "100-taxon" || "$folder" == "200-taxon" ]]; then
                TRUE_TREE="${DATASET_DIR%/}/$folder/true-species-trees/${REPL}/sp-cleaned"
            fi

            if [[ ! -f "$ALL_GT_FILE" ]]; then
                echo -e "    ${RED}Warning: Input file does not exist: $ALL_GT_FILE${NC}"
                continue
            fi

            # Loop through each algorithm
            for algorithm in "${ALGORITHMS[@]}"; do
                echo "      Processing algorithm: $algorithm"
                
                # output directory for this replicate and algorithm
                OUT_DIR="${DATASET_DIR%/}/$folder/$GT_FOLDER/${algorithm}_outputs"

                # run and write per-run stat file (skips if stat exists and not fresh)
                run_algorithm_and_write_stats "$ALL_GT_FILE" "$TRUE_TREE" "$OUT_DIR" "$REPL" "$folder" "$inner_folder" "$algorithm"
            done
        done
    done

    echo
done

echo -e "${GREEN}Dataset processing complete!${NC}"
echo "Check output directories for each algorithm ('stelar_outputs', 'astral_outputs', 'treeqmc_outputs', 'wqfmtree_outputs') for 'output-<alg>.tre' and 'stat-<alg>.csv' files."
