#!/usr/bin/env bash
# Multi-Algorithm Dataset Runner Script (updated)
# Supports STELAR-X and optional baseline algorithms
# Usage: ./scripts/run-bulk-standard.sh [--base-dir /path/to/base] [--dataset-dir /path/to/datasets] [--fresh]
#   --base-dir, -b    Optional base directory (defaults to value below)
#   --dataset-dir, -d Optional dataset directory (defaults to BASE_DIR/datasets)
#   --fresh           Force re-run even if stat-<alg>.csv exists

set -uo pipefail

# =============================================================================
# DEFAULTS (edit these if you want different defaults)
# =============================================================================
BASE_DIR="$HOME/phylogeny"  # default; can be overridden with --base-dir or -b
DATASET_DIR=""                        # dataset directory; defaults to BASE_DIR/datasets/standard
STELAR_X_ROOT=""                      # STELAR-X root (this project); derived from script location
ASTER_ROOT=""                         # ASTER root; derived from STELAR_X_ROOT if not set
ASTRAL_ROOT=""                        # derived from BASE_DIR if not set explicitly
TREEQMC_ROOT=""                       # derived from BASE_DIR if not set explicitly
WQFMTREE_ROOT=""                      # derived from BASE_DIR if not set explicitly
SUPERTRIPLETS_ROOT=""                 # SuperTriplets baseline; derived from STELAR_X_ROOT/baselines
TMC_ROOT=""                           # TMC baseline; derived from STELAR_X_ROOT/baselines
RUN_WITH_MONITOR_SCRIPT=""           # derived from STELAR_X_ROOT if not set
RUN_BASELINE_WITH_MONITOR_SCRIPT=""  # derived from STELAR_X_ROOT if not set
METHODS_ARG=""                       # optional semicolon-separated methods override
FOLDERS_ARG=""                       # optional semicolon-separated folders override
FRESH=false
NO_NOTIFY=false
GENERIC_OPTS=""
GENERIC_OPTS_SET=false
GENERIC_OPTS_LIST_RAW=""
GENERIC_OPTS_LIST_SET=false

# Algorithm configuration. Override with --method for a single method or a
# semicolon-separated sweep. Optional baselines are validated only if selected.
ALGORITHMS=("stelarx")



# Algorithm-specific options
# STELAR-X examples:
# GENERIC_OPTS="--search-space S2 --intersection-method I2 -vv"
# GENERIC_OPTS_LIST_RAW="--search-space S1 -vv;--search-space S2 -vv;--search-space S3 -vv"
# The setting-name encoder ignores verbosity. The single setting above becomes:
#   search-space_S2__intersection-method_I2
STELAR_OPTS="--search-space S2 -vv"
STELAR_OPTS_LIST_RAW=""
STELAR_OPTS_LIST=()
ASTER_OPTS="-t 16"  # ASTER thread count
ASTRAL_OPTS=""  # ASTRAL doesn't need special options for basic runs
TREEQMC_OPTS=""  # TreeQMC doesn't need special options for basic runs
WQFMTREE_OPTS=""  # wQFMtree doesn't need special options for basic runs
SUPERTRIPLETS_OPTS=""  # SuperTriplets options
TMC_OPTS=""  # TMC options

# Algorithm command mappings - these will be used to construct the actual commands
declare -A ALG_COMMANDS
ALG_COMMANDS["stelarx"]="stelarx_command"
ALG_COMMANDS["aster"]="aster_command"
ALG_COMMANDS["astral"]="astral_command"
ALG_COMMANDS["treeqmc"]="treeqmc_command"
ALG_COMMANDS["wqfmtree"]="wqfmtree_command"
ALG_COMMANDS["supertriplets"]="supertriplets_command"
ALG_COMMANDS["tmc"]="tmc_command"

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'
NTFY_CHANNEL_NAME="${NTFY_CHANNEL_NAME:-anik-phylo-stx}"

RUNNER_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${RUNNER_ROOT}/experiment-setting-name.sh"

# =============================================================================
# Dataset configuration
# =============================================================================
folders=("37-taxon")

declare -A innerFolderNames
innerFolderNames["11-taxon"]="estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS"
innerFolderNames["15-taxon"]="100gene-100bp/estimated-genetrees 100gene-1000bp/estimated-genetrees 100gene-true 1000gene-100bp/estimated-genetrees 1000gene-1000bp/estimated-genetrees 1000gene-true"


innerFolderNames["37-taxon"]="estimated-genetrees/0.5X-200-500 estimated-genetrees/1X-200-500 estimated-genetrees/1X-200-1000 estimated-genetrees/1X-400-500 estimated-genetrees/1X-400-1000 estimated-genetrees/1X-800-500 estimated-genetrees/1X-800-1000 estimated-genetrees/2X-200-500 true-genetrees/0.5X-200-true true-genetrees/1X-200-true true-genetrees/1X-400-true true-genetrees/1X-800-true true-genetrees/2X-200-true"


innerFolderNames["48-taxon"]="estimated-genetrees/1X-25-500 estimated-genetrees/1X-50-500 estimated-genetrees/1X-100-500 estimated-genetrees/1X-200-500 estimated-genetrees/1X-500-500 estimated-genetrees/1X-1000-500 estimated-genetrees/2X-1000-500"
innerFolderNames["48-taxon-latest"]="estimated_genetrees/1X-1000-500"
innerFolderNames["100-taxon"]="inner100"
innerFolderNames["200-taxon"]="inner200"

innerFolderNames["500-taxon"]="estimated-genetrees/model.500.2000000.0.000001/1000-gt true-genetrees/model.500.2000000.0.000001/1000-gt"
innerFolderNames["1000-taxon"]="estimated-genetrees/model.1000.2000000.0.000001/1000-gt true-genetrees/model.1000.2000000.0.000001/1000-gt"


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
replicates["500-taxon"]=20
replicates["1000-taxon"]=20
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
    echo "STELAR_X_ROOT: ${STELAR_X_ROOT:-(auto-detected)}"
    echo "ASTER_ROOT: ${ASTER_ROOT:-(derived from STELAR_X_ROOT)}"
    echo "ASTRAL_ROOT: ${ASTRAL_ROOT:-(derived from BASE_DIR)}"
    echo "TREEQMC_ROOT: ${TREEQMC_ROOT:-(derived from BASE_DIR)}"
    echo "WQFMTREE_ROOT: ${WQFMTREE_ROOT:-(derived from BASE_DIR)}"
    echo "SUPERTRIPLETS_ROOT: ${SUPERTRIPLETS_ROOT:-(derived from STELAR_X_ROOT)}"
    echo "TMC_ROOT: ${TMC_ROOT:-(derived from STELAR_X_ROOT)}"
    echo "Algorithms: ${ALGORITHMS[*]}"
    echo "STELAR_OPTS: ${STELAR_OPTS:-<empty>}"
    echo "ASTER_OPTS: ${ASTER_OPTS:-<empty>}"
    echo "ASTRAL_OPTS: ${ASTRAL_OPTS:-<empty>}"
    echo "TREEQMC_OPTS: ${TREEQMC_OPTS:-<empty>}"
    echo "WQFMTREE_OPTS: ${WQFMTREE_OPTS:-<empty>}"
    echo "SUPERTRIPLETS_OPTS: ${SUPERTRIPLETS_OPTS:-<empty>}"
    echo "TMC_OPTS: ${TMC_OPTS:-<empty>}"
    echo "Fresh run: $FRESH"
    echo "Notifications: $(if [[ "$NO_NOTIFY" = true ]]; then echo "disabled"; else echo "enabled"; fi)"
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

normalize_algorithm_name() {
    local name="${1,,}"
    case "$name" in
      stelarx|stelar|stelar-x) echo "stelarx" ;;
      aster) echo "aster" ;;
      astral) echo "astral" ;;
      treeqmc|tree-qmc) echo "treeqmc" ;;
      wqfmtree|wqfm-tree) echo "wqfmtree" ;;
      supertriplets) echo "supertriplets" ;;
      stp-nni) echo "stp-nni" ;;
      tmc) echo "tmc" ;;
      *) return 1 ;;
    esac
}

split_semicolon_list() {
    local raw="$1"
    local -n out_arr="$2"
    out_arr=()

    local -a raw_items=()
    local item=""
    IFS=';' read -r -a raw_items <<< "$raw"
    for item in "${raw_items[@]}"; do
      item="$(echo "$item" | sed 's/^ *//;s/ *$//')"
      if [[ -n "$item" ]]; then
        out_arr+=("$item")
      fi
    done
}

assign_generic_opts_to_algorithm() {
    local algorithm="$1"
    local value="$2"

    case "$algorithm" in
      stelarx) STELAR_OPTS="$value" ;;
      aster) ASTER_OPTS="$value" ;;
      astral) ASTRAL_OPTS="$value" ;;
      treeqmc) TREEQMC_OPTS="$value" ;;
      wqfmtree) WQFMTREE_OPTS="$value" ;;
      supertriplets|stp-nni) SUPERTRIPLETS_OPTS="$value" ;;
      tmc) TMC_OPTS="$value" ;;
      *)
        echo -e "${RED}Error: Unsupported algorithm '$algorithm' for generic --opts.${NC}"
        exit 1
        ;;
    esac
}

get_csv_column_value() {
    local csv_file="$1"
    local column_name="$2"
    awk -F',' -v col="$column_name" '
      NR==1 { for (i=1; i<=NF; i++) if ($i == col) c=i; next }
      NR==2 && c>0 { print $c; exit }
    ' "$csv_file" 2>/dev/null || true
}

csv_escape() {
    local raw="$1"
    raw="${raw//\"/\"\"}"
    printf '%s' "$raw"
}

strip_trivial_opts_for_csv() {
    local raw="$1"
    local -a tokens=()
    local -a kept=()
    local token

    if [[ -z "${raw// }" ]]; then
      printf ''
      return
    fi

    read -r -a tokens <<< "$raw"
    for token in "${tokens[@]}"; do
      case "$token" in
        -v|-vv|-vvv|-q|--quiet|--verbose)
          continue
          ;;
      esac
      kept+=("$token")
    done

    printf '%s' "${kept[*]}"
}

method_opts_for_algorithm() {
    local algorithm="$1"
    case "$algorithm" in
      stelarx) printf '%s' "$STELAR_OPTS" ;;
      aster) printf '%s' "$ASTER_OPTS" ;;
      astral) printf '%s' "$ASTRAL_OPTS" ;;
      treeqmc) printf '%s' "$TREEQMC_OPTS" ;;
      wqfmtree) printf '%s' "$WQFMTREE_OPTS" ;;
      supertriplets) printf '%s' "$SUPERTRIPLETS_OPTS" ;;
      stp-nni) printf '%s' "$SUPERTRIPLETS_OPTS" ;;
      tmc) printf '%s' "$TMC_OPTS" ;;
      *) printf '' ;;
    esac
}

send_run_notification() {
    local algorithm="$1"
    local folder="$2"
    local inner_folder="$3"
    local replicate="$4"
    local status="$5"
    local rf_rate="$6"
    local running_time="$7"
    local max_cpu_mb="$8"
    local max_gpu_mb="$9"
    local out_file="${10}"
    local method_opts="${11}"
    local stats_header="${12}"
    local stats_row="${13}"
    local triplet_score="${14:-NA}"

    if [[ "$NO_NOTIFY" = true ]]; then
      return 0
    fi
    if ! command -v curl >/dev/null 2>&1; then
      return 0
    fi

    local status_emoji="✅"
    if [[ "$status" != "SUCCESS" ]]; then
      status_emoji="❌"
    fi

    local score_line=""
    if [[ "$algorithm" == "stelarx" ]]; then
      score_line="Triplet score: ${triplet_score}"$'\n'
    fi
    local notify_msg
    notify_msg="${status_emoji} ${algorithm} ${status}

Dataset: ${folder}/${inner_folder}/${replicate}
Opts: ${method_opts}
RF rate: ${rf_rate}
${score_line}Time: ${running_time}s
Time: ${running_time}s
CPU RAM: ${max_cpu_mb} MB
GPU VRAM: ${max_gpu_mb} MB

Output: $(basename "$out_file")

Stats CSV row:
${stats_header}
${stats_row}"

    curl -s -d "$notify_msg" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
}

validate_algorithm_binaries() {
    local errors_found=false
    local needs_baseline_wrapper=false

    if [[ ! -d "$STELAR_X_ROOT" ]]; then
        echo -e "${RED}Error: STELAR_X_ROOT '$STELAR_X_ROOT' does not exist.${NC}"
        errors_found=true
    fi

    local alg
    for alg in "${ALGORITHMS[@]}"; do
        case "$alg" in
          stelarx)
            if [[ ! -x "$RUN_WITH_MONITOR_SCRIPT" ]]; then
              echo -e "${RED}Error: STELAR-X monitor wrapper not found: $RUN_WITH_MONITOR_SCRIPT${NC}"
              errors_found=true
            fi
            if [[ ! -x "${STELAR_X_ROOT}/stelarx" ]]; then
              echo -e "${RED}Error: STELAR-X launcher (stelarx) not found in project root.${NC}"
              errors_found=true
            fi
            ;;
          aster|astral|treeqmc|wqfmtree|supertriplets|tmc)
            needs_baseline_wrapper=true
            ;;
          stp-nni)
            if [[ ! -x "$RUN_BASELINE_WITH_MONITOR_SCRIPT" ]]; then
              echo -e "${RED}Error: Baseline monitor wrapper not found: $RUN_BASELINE_WITH_MONITOR_SCRIPT${NC}"
              errors_found=true
            fi
            if [[ ! -x "${STELAR_X_ROOT}/run_supertriplets.sh" ]]; then
              echo -e "${RED}Error: run_supertriplets.sh not found in STELAR_X_ROOT.${NC}"
              errors_found=true
            fi
            if [[ ! -f "${SUPERTRIPLETS_ROOT}/SuperTriplets_v1.1.jar" ]]; then
              echo -e "${RED}Error: SuperTriplets JAR not found in SUPERTRIPLETS_ROOT.${NC}"
              errors_found=true
            fi
            ;;
        esac
    done

    if [[ "$needs_baseline_wrapper" = true && ! -d "${STELAR_X_ROOT%/}/baselines" ]]; then
      echo -e "${RED}Error: baselines directory not found in ${STELAR_X_ROOT%/}.${NC}"
      errors_found=true
    fi
    if [[ "$needs_baseline_wrapper" = true && ! -x "$RUN_BASELINE_WITH_MONITOR_SCRIPT" ]]; then
      echo -e "${RED}Error: Baseline monitor wrapper not found: $RUN_BASELINE_WITH_MONITOR_SCRIPT${NC}"
      errors_found=true
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
    local ALGORITHM="$7"         # algorithm name (stelarx, astral, etc.)

    local START_NS
    local END_NS
    local ELAPSED_MS
    local RUNNING_TIME
    local MAX_CPU_MB="NA"
    local MAX_GPU_MB="NA"
    local TRIPLET_SCORE="NA"
    local ALGORITHM_EXIT_CODE=0

    START_NS=$(date +%s%N)

    local TEMP_STP_OUTPUT=""
    local TEMP_WRAPPER_STATS_TO_DELETE=""
    local METHOD_OPTS_RAW
    METHOD_OPTS_RAW="$(method_opts_for_algorithm "$ALGORITHM")"
    local METHOD_OPTS_CSV_RAW
    METHOD_OPTS_CSV_RAW="$(strip_trivial_opts_for_csv "$METHOD_OPTS_RAW")"
    local METHOD_OPTS_ESCAPED
    METHOD_OPTS_ESCAPED="$(csv_escape "$METHOD_OPTS_CSV_RAW")"
    local SETTING_NAME
    SETTING_NAME="$(build_setting_name_from_opts "$METHOD_OPTS_RAW")"
    if [[ "$ALGORITHM" == "stelarx" ]]; then
      OUT_DIR="${OUT_DIR%/}/${SETTING_NAME}"
    fi

    mkdir -p "$OUT_DIR"

    local OUT_FILE="${OUT_DIR%/}/output-${ALGORITHM}.tre"
    local STAT_FILE="${OUT_DIR%/}/stat-${ALGORITHM}.csv"
    local WRAPPER_STATS_FILE="${OUT_FILE%.tre}_stats.csv"

    if [[ "$FRESH" = false && -f "$STAT_FILE" ]]; then
        echo "      SKIPPING: ${STAT_FILE} already exists. Use --fresh to force rerun."
        return 0
    fi

    echo "      Running ${ALGORITHM^^} [${SETTING_NAME}] (output -> $OUT_FILE)"
    case "$ALGORITHM" in
      "stelarx")
        if [[ -n "$STELAR_OPTS" ]]; then
          "$RUN_WITH_MONITOR_SCRIPT" \
            --stelar-root "$STELAR_X_ROOT" \
            --input "$ALL_GT_FILE" \
            --output "$OUT_FILE" \
            --opts "$STELAR_OPTS" \
            --no-notify
        else
          "$RUN_WITH_MONITOR_SCRIPT" \
            --stelar-root "$STELAR_X_ROOT" \
            --input "$ALL_GT_FILE" \
            --output "$OUT_FILE" \
            --no-notify
        fi
        ALGORITHM_EXIT_CODE=$?
        ;;
      "aster"|"astral"|"treeqmc"|"wqfmtree"|"supertriplets"|"tmc")
        local baseline_opts=()
        case "$ALGORITHM" in
          aster) [[ -n "$ASTER_OPTS" ]] && baseline_opts+=(--aster-opts "$ASTER_OPTS") ;;
          astral) [[ -n "$ASTRAL_OPTS" ]] && baseline_opts+=(--astral-opts "$ASTRAL_OPTS") ;;
          treeqmc) [[ -n "$TREEQMC_OPTS" ]] && baseline_opts+=(--treeqmc-opts "$TREEQMC_OPTS") ;;
          wqfmtree) [[ -n "$WQFMTREE_OPTS" ]] && baseline_opts+=(--wqfm-opts "$WQFMTREE_OPTS") ;;
          supertriplets) [[ -n "$SUPERTRIPLETS_OPTS" ]] && baseline_opts+=(--supertriplets-opts "$SUPERTRIPLETS_OPTS") ;;
          tmc) [[ -n "$TMC_OPTS" ]] && baseline_opts+=(--tmc-opts "$TMC_OPTS") ;;
        esac

        "$RUN_BASELINE_WITH_MONITOR_SCRIPT" \
          --stelar-root "$STELAR_X_ROOT" \
          --method "$ALGORITHM" \
          --input "$ALL_GT_FILE" \
          --output "$OUT_FILE" \
          --no-notify \
          "${baseline_opts[@]}"
        ALGORITHM_EXIT_CODE=$?
        ;;
      "stp-nni")
        TEMP_STP_OUTPUT=$(mktemp --suffix=.tre)
        WRAPPER_STATS_FILE="${TEMP_STP_OUTPUT%.tre}_stats.csv"
        TEMP_WRAPPER_STATS_TO_DELETE="$WRAPPER_STATS_FILE"

        local stp_opts=()
        [[ -n "$SUPERTRIPLETS_OPTS" ]] && stp_opts+=(--supertriplets-opts "$SUPERTRIPLETS_OPTS")

        "$RUN_BASELINE_WITH_MONITOR_SCRIPT" \
          --stelar-root "$STELAR_X_ROOT" \
          --method "supertriplets" \
          --input "$ALL_GT_FILE" \
          --output "$TEMP_STP_OUTPUT" \
          --no-notify \
          "${stp_opts[@]}"
        ALGORITHM_EXIT_CODE=$?

        if [[ $ALGORITHM_EXIT_CODE -eq 0 ]]; then
          local NNI_FILE="${TEMP_STP_OUTPUT}.nni"
          if [[ -f "$NNI_FILE" ]]; then
            cp "$NNI_FILE" "$OUT_FILE"
            echo "      Copied NNI tree to $OUT_FILE"
          elif [[ -f "$TEMP_STP_OUTPUT" ]]; then
            cp "$TEMP_STP_OUTPUT" "$OUT_FILE"
            echo "      Warning: NNI output not found, using SuperTriplets tree instead."
          fi
        fi
        rm -f "$TEMP_STP_OUTPUT" "${TEMP_STP_OUTPUT}.nni" 2>/dev/null || true
        ;;
      *)
        echo "Unknown algorithm: $ALGORITHM"
        return 1
        ;;
    esac

    END_NS=$(date +%s%N)
    ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
    RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", ${ELAPSED_MS}/1000}")

    if [[ -f "$WRAPPER_STATS_FILE" ]]; then
      local WRAPPER_RUNNING_TIME
      local WRAPPER_MAX_CPU
      local WRAPPER_MAX_GPU
      local WRAPPER_EXIT_CODE
      local WRAPPER_TRIPLET_SCORE
      WRAPPER_RUNNING_TIME=$(get_csv_column_value "$WRAPPER_STATS_FILE" "running_time_s")
      WRAPPER_MAX_CPU=$(get_csv_column_value "$WRAPPER_STATS_FILE" "max_cpu_mb")
      WRAPPER_MAX_GPU=$(get_csv_column_value "$WRAPPER_STATS_FILE" "max_gpu_mb")
      WRAPPER_TRIPLET_SCORE=$(get_csv_column_value "$WRAPPER_STATS_FILE" "optimal_triplet_score")
      WRAPPER_EXIT_CODE=$(get_csv_column_value "$WRAPPER_STATS_FILE" "exit_code")

      [[ -n "$WRAPPER_RUNNING_TIME" ]] && RUNNING_TIME="$WRAPPER_RUNNING_TIME"
      [[ -n "$WRAPPER_MAX_CPU" ]] && MAX_CPU_MB="$WRAPPER_MAX_CPU"
      [[ -n "$WRAPPER_MAX_GPU" ]] && MAX_GPU_MB="$WRAPPER_MAX_GPU"
      [[ -n "$WRAPPER_TRIPLET_SCORE" ]] && TRIPLET_SCORE="$WRAPPER_TRIPLET_SCORE"
      if [[ -n "$WRAPPER_EXIT_CODE" && "$WRAPPER_EXIT_CODE" =~ ^-?[0-9]+$ ]]; then
        ALGORITHM_EXIT_CODE="$WRAPPER_EXIT_CODE"
      fi
    fi

    if [[ -n "$TEMP_WRAPPER_STATS_TO_DELETE" ]]; then
      rm -f "$TEMP_WRAPPER_STATS_TO_DELETE" 2>/dev/null || true
    fi

    # Check if algorithm failed
    if [[ $ALGORITHM_EXIT_CODE -ne 0 ]]; then
      echo -e "      ${RED}ERROR: ${ALGORITHM^^} failed with exit code $ALGORITHM_EXIT_CODE${NC}"
      echo -e "      ${RED}Skipping CSV writing and RF calculation for this run${NC}"
      echo -e "      ${YELLOW}Continuing with next dataset...${NC}"
      local FAIL_STATS_HEADER="alg,setting,opts,folder,inner_folder,replicate,rf-rate,optimal-triplet-score,running-time-s,max-cpu-mb,max-gpu-mb,exit-code"
      local FAIL_STATS_ROW="${ALGORITHM},${SETTING_NAME},\"${METHOD_OPTS_ESCAPED}\",${FOLDER_NAME},${INNER_FOLDER_NAME},${REPLICATE_NAME},NA,${TRIPLET_SCORE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${ALGORITHM_EXIT_CODE}"
      send_run_notification "$ALGORITHM" "$FOLDER_NAME" "$INNER_FOLDER_NAME" "$REPLICATE_NAME" "FAILED" "NA" "$RUNNING_TIME" "$MAX_CPU_MB" "$MAX_GPU_MB" "$OUT_FILE" "${METHOD_OPTS_RAW:-<default/none>}" "$FAIL_STATS_HEADER" "$FAIL_STATS_ROW" "$TRIPLET_SCORE"
      return 1
    fi

    # RF calculation
    local RF_RATE="NA"
    if [[ -f "$OUT_FILE" ]]; then
      # Prefer rf.py inside STELAR_X_ROOT if present
      if [[ -f "${STELAR_X_ROOT%/}/scripts/rf.py" && -x "$PYTHON_BIN" ]]; then
        echo "      Calculating RF using ${STELAR_X_ROOT%/}/scripts/rf.py"
        rf_output=$("$PYTHON_BIN" "${STELAR_X_ROOT%/}/scripts/rf.py" "$OUT_FILE" "$TRUE_SPECIES_TREE" 2>&1) || rf_output="$rf_output"
        
        echo "$rf_output"
        
        # try to find a number in output
        rf_candidate=$(echo "$rf_output" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
        if [[ -n "$rf_candidate" ]]; then
          RF_RATE="$rf_candidate"
        else
          RF_RATE="NA"
        fi

        echo ""
        echo "RF rate: ${RF_RATE}"
        echo ""

      elif [[ -f "${BASE_DIR%/}/RF/getFpFn.py" && -x "$PYTHON_BIN" ]]; then
        echo "      Calculating RF using ${BASE_DIR%/}/RF/getFpFn.py"
        # getFpFn.py returns a tuple; attempt to extract the same metric as earlier script
        tuple=$("$PYTHON_BIN" "${BASE_DIR%/}/RF/getFpFn.py" -e "$OUT_FILE" -t "$TRUE_SPECIES_TREE" 2>/dev/null) || tuple=""
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
    # Header: alg,opts,folder,inner_folder,replicate,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb,exit-code
    local STATS_HEADER="alg,setting,opts,folder,inner_folder,replicate,rf-rate,optimal-triplet-score,running-time-s,max-cpu-mb,max-gpu-mb,exit-code"
    local STATS_ROW="${ALGORITHM},${SETTING_NAME},\"${METHOD_OPTS_ESCAPED}\",${FOLDER_NAME},${INNER_FOLDER_NAME},${REPLICATE_NAME},${RF_RATE},${TRIPLET_SCORE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},${ALGORITHM_EXIT_CODE}"
    {
      echo "$STATS_HEADER"
      echo "$STATS_ROW"
    } > "$STAT_FILE"

    echo -e "      ${GREEN}Wrote stats to $STAT_FILE${NC}"
    send_run_notification "$ALGORITHM" "$FOLDER_NAME" "$INNER_FOLDER_NAME" "$REPLICATE_NAME" "SUCCESS" "$RF_RATE" "$RUNNING_TIME" "$MAX_CPU_MB" "$MAX_GPU_MB" "$OUT_FILE" "${METHOD_OPTS_RAW:-<default/none>}" "$STATS_HEADER" "$STATS_ROW" "$TRIPLET_SCORE"
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
    --stelar-x-root)
      STELAR_X_ROOT="$2"
      shift 2
      ;;
    --method|-m)
      METHODS_ARG="$2"
      shift 2
      ;;
    --folder|-f)
      FOLDERS_ARG="$2"
      shift 2
      ;;
    --opts|--alg-opts)
      GENERIC_OPTS="$2"
      GENERIC_OPTS_SET=true
      shift 2
      ;;
    --opts=*|--alg-opts=*)
      GENERIC_OPTS="${1#*=}"
      GENERIC_OPTS_SET=true
      shift
      ;;
    --opts-list|--alg-opts-list)
      GENERIC_OPTS_LIST_RAW="$2"
      GENERIC_OPTS_LIST_SET=true
      shift 2
      ;;
    --opts-list=*|--alg-opts-list=*)
      GENERIC_OPTS_LIST_RAW="${1#*=}"
      GENERIC_OPTS_LIST_SET=true
      shift
      ;;
    --stelar-opts|--stelarx-opts)
      STELAR_OPTS="$2"
      shift 2
      ;;
    --stelar-opts=*|--stelarx-opts=*)
      STELAR_OPTS="${1#*=}"
      shift
      ;;
    --stelar-opts-list|--stelarx-opts-list)
      STELAR_OPTS_LIST_RAW="$2"
      shift 2
      ;;
    --stelar-opts-list=*|--stelarx-opts-list=*)
      STELAR_OPTS_LIST_RAW="${1#*=}"
      shift
      ;;
    --aster-opts)
      ASTER_OPTS="$2"
      shift 2
      ;;
    --aster-opts=*)
      ASTER_OPTS="${1#*=}"
      shift
      ;;
    --astral-opts)
      ASTRAL_OPTS="$2"
      shift 2
      ;;
    --astral-opts=*)
      ASTRAL_OPTS="${1#*=}"
      shift
      ;;
    --treeqmc-opts)
      TREEQMC_OPTS="$2"
      shift 2
      ;;
    --treeqmc-opts=*)
      TREEQMC_OPTS="${1#*=}"
      shift
      ;;
    --wqfm-opts|--wqfmtree-opts)
      WQFMTREE_OPTS="$2"
      shift 2
      ;;
    --wqfm-opts=*|--wqfmtree-opts=*)
      WQFMTREE_OPTS="${1#*=}"
      shift
      ;;
    --supertriplets-opts)
      SUPERTRIPLETS_OPTS="$2"
      shift 2
      ;;
    --supertriplets-opts=*)
      SUPERTRIPLETS_OPTS="${1#*=}"
      shift
      ;;
    --tmc-opts)
      TMC_OPTS="$2"
      shift 2
      ;;
    --tmc-opts=*)
      TMC_OPTS="${1#*=}"
      shift
      ;;
    --fresh)
      FRESH=true
      shift
      ;;
    --no-notify|-nn)
      NO_NOTIFY=true
      shift
      ;;
    --help|-h)
      cat <<EOF
Usage: $0 [--base-dir /path/to/base] [--dataset-dir /path/to/datasets] [--method "m1;m2"] [--folder "f1;f2"] [--fresh] [--no-notify]

Multi-algorithm dataset runner supporting STELAR-X, ASTER, ASTRAL, TreeQMC, wQFMtree, SuperTriplets, and TMC.

--base-dir, -b      Base directory containing RF/ and external tools (overrides default)
--dataset-dir, -d   Dataset directory (overrides default BASE_DIR/datasets)
--stelar-x-root     STELAR-X root directory (overrides auto-detection from script location)
--method, -m        Optional semicolon-separated methods (e.g. "stelar;astral;aster")
--folder, -f        Optional semicolon-separated folders (e.g. "37-taxon;48-taxon")
--opts, --alg-opts  Override the default option string for the selected algorithm
--opts-list, --alg-opts-list
                    Semicolon-separated list of option strings to loop over for the selected algorithm
--stelar-opts       Compatibility alias for STELAR-X-specific --opts
--stelar-opts-list  Compatibility alias for STELAR-X-specific --opts-list
--aster-opts        Override default ASTER_OPTS
--astral-opts       Override default ASTRAL_OPTS
--treeqmc-opts      Override default TREEQMC_OPTS
--wqfm-opts         Override default WQFMTREE_OPTS
--supertriplets-opts Override default SUPERTRIPLETS_OPTS
--tmc-opts          Override default TMC_OPTS
--fresh             Force rerun even if stat-<alg>.csv exists
--no-notify, -nn    Disable bulk-level ntfy notifications
--help, -h          Show this help

Algorithms available: stelarx, aster, astral, treeqmc, wqfmtree, supertriplets, stp-nni, tmc
Example STELAR-X setting sweep:
  --method "stelarx" --opts-list "--search-space S1 -vv;--search-space S2 -vv;--search-space S3 -vv"
Algorithm root directories:
  STELAR-X:       Auto-detected from script location
  ASTER:          \${STELAR_X_ROOT}/baselines/ASTER
  ASTRAL:         \${STELAR_X_ROOT}/baselines/ASTRAL
  TreeQMC:        \${STELAR_X_ROOT}/baselines/TREE-QMC
  wQFMtree:       \${STELAR_X_ROOT}/baselines/wQFM-TREE
  SuperTriplets:  \${STELAR_X_ROOT}/baselines/SuperTriplets
  TMC:            \${STELAR_X_ROOT}/baselines/TMC
EOF
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

if [[ -n "$METHODS_ARG" ]]; then
  split_semicolon_list "$METHODS_ARG" user_methods_raw
  if [[ ${#user_methods_raw[@]} -eq 0 ]]; then
    echo -e "${RED}Error: --method was provided but no valid values were found.${NC}"
    exit 1
  fi

  ALGORITHMS=()
  local_method=""
  for local_method in "${user_methods_raw[@]}"; do
    normalized_method=$(normalize_algorithm_name "$local_method") || {
      echo -e "${RED}Error: Unsupported method '$local_method'.${NC}"
      echo "Supported methods: stelarx, aster, astral, treeqmc, wqfmtree, supertriplets, stp-nni, tmc"
      exit 1
    }
    ALGORITHMS+=("$normalized_method")
  done
fi

if [[ "$GENERIC_OPTS_SET" == true || "$GENERIC_OPTS_LIST_SET" == true ]]; then
  if [[ ${#ALGORITHMS[@]} -ne 1 ]]; then
    echo -e "${RED}Error: --opts and --opts-list require exactly one selected algorithm in --method.${NC}"
    exit 1
  fi

  assign_generic_opts_to_algorithm "${ALGORITHMS[0]}" "$GENERIC_OPTS"

  if [[ "$GENERIC_OPTS_LIST_SET" == true ]]; then
    case "${ALGORITHMS[0]}" in
      stelarx)
        STELAR_OPTS_LIST_RAW="$GENERIC_OPTS_LIST_RAW"
        ;;
      *)
        echo -e "${RED}Error: --opts-list is currently supported only for stelarx in this script.${NC}"
        exit 1
        ;;
    esac
  fi
fi

if [[ -n "$STELAR_OPTS_LIST_RAW" ]]; then
  split_semicolon_list "$STELAR_OPTS_LIST_RAW" STELAR_OPTS_LIST
fi
if [[ ${#STELAR_OPTS_LIST[@]} -eq 0 ]]; then
  STELAR_OPTS_LIST=("$STELAR_OPTS")
fi

if [[ -n "$FOLDERS_ARG" ]]; then
  split_semicolon_list "$FOLDERS_ARG" user_folders
  if [[ ${#user_folders[@]} -eq 0 ]]; then
    echo -e "${RED}Error: --folder was provided but no valid values were found.${NC}"
    exit 1
  fi
  folders=("${user_folders[@]}")
fi

# derive STELAR_X_ROOT from script location if not set
if [[ -z "${STELAR_X_ROOT}" ]]; then
  STELAR_X_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

PYTHON_BIN="${STELARX_PYTHON:-${STELAR_X_ROOT%/}/.venv/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || true)"
fi
if [[ ! -x "$PYTHON_BIN" ]]; then
  echo -e "${RED}Error: Python 3 is required for RF calculations.${NC}"
  exit 1
fi

# derive DATASET_DIR from BASE_DIR if not set
if [[ -z "${DATASET_DIR}" ]]; then
  DATASET_DIR="${BASE_DIR%/}/datasets/standard"
fi

# derive ASTER_ROOT from STELAR_X_ROOT if not set
if [[ -z "${ASTER_ROOT}" ]]; then
  ASTER_ROOT="${STELAR_X_ROOT%/}/baselines/ASTER"
fi

# derive ASTRAL_ROOT from BASE_DIR if not set
if [[ -z "${ASTRAL_ROOT}" ]]; then
  ASTRAL_ROOT="${STELAR_X_ROOT%/}/baselines/ASTRAL"
fi

# derive TREEQMC_ROOT from BASE_DIR if not set
if [[ -z "${TREEQMC_ROOT}" ]]; then
  TREEQMC_ROOT="${STELAR_X_ROOT%/}/baselines/TREE-QMC"
fi

# derive WQFMTREE_ROOT from BASE_DIR if not set
if [[ -z "${WQFMTREE_ROOT}" ]]; then
  WQFMTREE_ROOT="${STELAR_X_ROOT%/}/baselines/wQFM-TREE"
fi

# derive SUPERTRIPLETS_ROOT from STELAR_X_ROOT if not set
if [[ -z "${SUPERTRIPLETS_ROOT}" ]]; then
  SUPERTRIPLETS_ROOT="${STELAR_X_ROOT%/}/baselines/SuperTriplets"
fi

# derive TMC_ROOT from STELAR_X_ROOT if not set
if [[ -z "${TMC_ROOT}" ]]; then
  TMC_ROOT="${STELAR_X_ROOT%/}/baselines/TMC"
fi

# derive monitor wrappers from STELAR_X_ROOT
if [[ -z "${RUN_WITH_MONITOR_SCRIPT}" ]]; then
  RUN_WITH_MONITOR_SCRIPT="${STELAR_X_ROOT%/}/scripts/run-stelarx-with-monitor.sh"
fi
if [[ -z "${RUN_BASELINE_WITH_MONITOR_SCRIPT}" ]]; then
  RUN_BASELINE_WITH_MONITOR_SCRIPT="${STELAR_X_ROOT%/}/baselines/run-baseline-with-monitor.sh"
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

    if [[ -z "${innerFolderNames[$folder]+x}" ]]; then
        echo -e "${RED}Warning: Folder key '$folder' is not configured in innerFolderNames, skipping...${NC}"
        continue
    fi

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
            
            if [[ "$folder" == "37-taxon" && "$inner_folder" == true-genetrees* ]]; then
                ALL_GT_FILE="${DATASET_DIR%/}/$folder/$GT_FOLDER/all_gt.tre"
            fi

            if [[ "$folder" == "48-taxon" ]]; then
                ALL_GT_FILE="${DATASET_DIR%/}/$folder/$GT_FOLDER/all_gt.tre"
            fi

            # special handling for 100/200-taxon true tree path (kept from your original)
            if [[ "$folder" == "100-taxon" || "$folder" == "200-taxon" ]]; then
                TRUE_TREE="${DATASET_DIR%/}/$folder/true-species-trees/${REPL}/sp-cleaned"
            fi

            # # special handling for 500/1000-taxon true tree path (kept from your original)
            # if [[ "$folder" == "500-taxon" || "$folder" == "1000-taxon" ]]; then
            #     # here 500 or 1000 should be conditional

            #     TRUE_TREE="${DATASET_DIR%/}/$folder/true-species-trees/model.1000.2000000.0.000001/${REPL}/sp-cleaned"

            #     # ALL_GT_FILE="${DATASET_DIR%/}/$folder/$GT_FOLDER/gt-cleaned"

            #     # ./scripts/process_unrooted.sh -i ${ALL_GT_FILE} -o ${ALL_GT_FILE}-rooted-og-0.tre -ogs "0"

            #     ALL_GT_FILE="${DATASET_DIR%/}/$folder/$GT_FOLDER/gt-cleaned-rooted-og-0.tre"

            #     # CURRENT_GT_FOLDER="${DATASET_DIR%/}/$folder/${inner_folder}/${REPL}"

            #     # # NEW_50_GT_FOLDER="${DATASET_DIR%/}/$folder/${inner_folder}/50-gt/${REPL}"
            #     # # NEW_200_GT_FOLDER="${DATASET_DIR%/}/$folder/${inner_folder}/200-gt/${REPL}"
            #     # NEW_1000_GT_FOLDER="${DATASET_DIR%/}/$folder/${inner_folder}/1000-gt/${REPL}"

            #     # # mkdir -p ${NEW_50_GT_FOLDER} ${NEW_200_GT_FOLDER} ${NEW_1000_GT_FOLDER}

            #     # mkdir -p ${NEW_1000_GT_FOLDER}

            #     # cp ${CURRENT_GT_FOLDER}/gt-cleaned ${NEW_1000_GT_FOLDER}/gt-cleaned
            #     # # cp ${CURRENT_GT_FOLDER}/gt-cleaned-50 ${NEW_50_GT_FOLDER}/gt-cleaned
            #     # # cp ${CURRENT_GT_FOLDER}/gt-cleaned-200 ${NEW_200_GT_FOLDER}/gt-cleaned

            #     # continue
            # fi

            # special handling for 500/1000-taxon true tree path (kept from your original)
            if [[ "$folder" == "500-taxon" ]]; then
                TRUE_TREE="${DATASET_DIR%/}/$folder/true-species-trees/model.500.2000000.0.000001/${REPL}/sp-cleaned"
                ALL_GT_FILE="${DATASET_DIR%/}/$folder/$GT_FOLDER/gt-cleaned-rooted-og-0.tre"
            fi

            if [[ "$folder" == "1000-taxon" ]]; then
                TRUE_TREE="${DATASET_DIR%/}/$folder/true-species-trees/model.1000.2000000.0.000001/${REPL}/sp-cleaned"
                ALL_GT_FILE="${DATASET_DIR%/}/$folder/$GT_FOLDER/gt-cleaned-rooted-og-0.tre"
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

                if [[ "$algorithm" == "stelarx" ]]; then
                    for STELARX_OPTS_ITEM in "${STELAR_OPTS_LIST[@]}"; do
                        STELAR_OPTS="$STELARX_OPTS_ITEM"
                        if ! run_algorithm_and_write_stats "$ALL_GT_FILE" "$TRUE_TREE" "$OUT_DIR" "$REPL" "$folder" "$inner_folder" "$algorithm"; then
                            echo -e "      ${RED}Failed to process $algorithm for ${folder}/${inner_folder}/${REPL}${NC}"
                            echo -e "      ${YELLOW}Continuing with next setting or algorithm...${NC}"
                            continue
                        fi
                    done
                else
                    if ! run_algorithm_and_write_stats "$ALL_GT_FILE" "$TRUE_TREE" "$OUT_DIR" "$REPL" "$folder" "$inner_folder" "$algorithm"; then
                        echo -e "      ${RED}Failed to process $algorithm for ${folder}/${inner_folder}/${REPL}${NC}"
                        echo -e "      ${YELLOW}Continuing with next algorithm...${NC}"
                        continue
                    fi
                fi
            done
        done
    done

    echo
done

echo -e "${GREEN}Dataset processing complete!${NC}"
echo "Check output directories for 'output-<alg>.tre' and 'stat-<alg>.csv' files."
echo "Output directories: stelarx_outputs, aster_outputs, astral_outputs, treeqmc_outputs, wqfmtree_outputs, supertriplets_outputs, stp-nni_outputs, tmc_outputs"
