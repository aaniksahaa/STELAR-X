#!/usr/bin/env bash
# test-baseline-simulated.sh
# Unified simulated-data runner for baseline methods (not STELAR-X).
# Uses run-baseline-with-monitor.sh and writes stat-<method>.csv per replicate.
#
# Usage:
#   ./test-baseline-simulated.sh -m aster -t 1000 -g 500
#   ./test-baseline-simulated.sh -m astral -t 1000 -g 500 --astral-opts "-t 8"
#
set -euo pipefail

# Defaults
METHOD=""
TAXA_NUM=""
GENE_TREES=""
REPLICATE="R1"
BASE_DIR=".."
SIMPHY_DIR=""
SIMPHY_DIR_SET=false
SIMPHY_DATA_DIR=""
SIMPHY_DATA_DIR_SET=false
STELAR_ROOT=""
STELAR_ROOT_SET=false

# Defaults that match run_simulator.sh
SB="0.000001"
SPMIN="500000"
SPMAX="1500000"

USE_LEGACY_LAYOUT=false
FRESH=false

# Method options (passed to run-baseline-with-monitor.sh)
ASTER_OPTS=""
ASTER_BIN=""
ASTRAL_OPTS=""
ASTRAL_XMS=""
ASTRAL_XMX=""
TREEQMC_OPTS=""
WQFM_OPTS=""
SUPERTRIPLETS_OPTS=""
TMC_OPTS=""

# Monitoring options (DEFAULT: ON)
TIME_MONITOR=true
GPU_MONITOR=true
NO_NOTIFY=false
DEBUG=0

print_help() {
  cat <<EOF
test-baseline-simulated.sh

Required:
  --method, -m      Method: aster, astral, treeqmc, wqfmtree, supertriplets, tmc
  --taxa_num, -t    Number of taxa (e.g. 1000)
  --gene_trees, -g  Number of gene trees (e.g. 500)

Optional:
  --replicate, -r       Replicate name (default: R1)
  --base-dir, -b        Base directory (default: ${BASE_DIR})
  --simphy-dir          Path to simphy dir (overrides --base-dir)
  --simphy-data-dir     Custom directory for simphy data storage
  --stelar-root         Path to STELAR-X root (overrides --base-dir)
  --sb                  Substitution/birthrate parameter (default: ${SB})
  --spmin               Population size minimum (default: ${SPMIN})
  --spmax               Population size maximum (default: ${SPMAX})
  --use-legacy-layout   Use legacy simphy layout
  --fresh               Force rerun even if stat-<method>.csv exists
  --no-time-monitor     Disable time-monitoring
  --no-gpu-monitor      Disable GPU-monitoring
  --no-notify, -nn      Disable ntfy.sh notifications
  --debug               Enable shell tracing

Method Options:
  --aster-opts "..."         Extra options for ASTER
  --aster-bin PATH          Path to ASTER binary (relative to ASTER root or absolute)
  --astral-opts "..."        Extra options for ASTRAL
  --astral-xms SIZE          ASTRAL Java Xms (e.g., 4g)
  --astral-xmx SIZE          ASTRAL Java Xmx (e.g., 128g)
  --astral-Xms SIZE          Same as --astral-xms
  --astral-Xmx SIZE          Same as --astral-xmx
  --treeqmc-opts "..."       Extra options for TreeQMC
  --wqfm-opts "..."          Extra options for wQFM-TREE
  --supertriplets-opts "..." Extra options for SuperTriplets wrapper
  --tmc-opts "..."           Extra options for TMC wrapper

EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --method|-m) METHOD="$2"; shift 2 ;;
    --taxa_num|-t) TAXA_NUM="$2"; shift 2 ;;
    --gene_trees|-g) GENE_TREES="$2"; shift 2 ;;
    --replicate|-r) REPLICATE="$2"; shift 2 ;;
    --simphy-dir) SIMPHY_DIR="$2"; SIMPHY_DIR_SET=true; shift 2 ;;
    --simphy-data-dir) SIMPHY_DATA_DIR="$2"; SIMPHY_DATA_DIR_SET=true; shift 2 ;;
    --stelar-root) STELAR_ROOT="$2"; STELAR_ROOT_SET=true; shift 2 ;;
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --sb) SB="$2"; shift 2 ;;
    --spmin) SPMIN="$2"; shift 2 ;;
    --spmax) SPMAX="$2"; shift 2 ;;
    --use-legacy-layout) USE_LEGACY_LAYOUT=true; shift ;;
    --fresh) FRESH=true; shift ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
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
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$METHOD" || -z "$TAXA_NUM" || -z "$GENE_TREES" ]]; then
  echo "Error: --method, --taxa_num and --gene_trees are required."
  print_help
  exit 2
fi

# Normalize method name
case "$METHOD" in
  aster|astral|treeqmc|tree-qmc|wqfmtree|wqfm-tree|supertriplets|tmc) ;;
  *)
    echo "Error: unknown method '$METHOD'"
    exit 3
    ;;
esac
if [[ "$METHOD" == "tree-qmc" ]]; then
  METHOD="treeqmc"
fi
if [[ "$METHOD" == "wqfm-tree" ]]; then
  METHOD="wqfmtree"
fi

# Derive SIMPHY_DIR/STELAR_ROOT from BASE_DIR if not explicitly set
if [[ "$SIMPHY_DIR_SET" = false ]]; then
  SIMPHY_DIR="${BASE_DIR%/}/STELAR-X/simphy"
fi
if [[ "$STELAR_ROOT_SET" = false ]]; then
  STELAR_ROOT="${BASE_DIR%/}/STELAR-X"
fi

PAIR="${TAXA_NUM}_${GENE_TREES}"

# Choose simphy run dir using SIMPHY_DATA_DIR if provided, otherwise use SIMPHY_DIR/data
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

STAT_FILE="${SIMPHY_RUN_DIR%/}/stat-${METHOD}.csv"
ALL_GT_FILE="${SIMPHY_RUN_DIR%/}/all_gt.tre"
TRUE_SPECIES_TREE="${SIMPHY_RUN_DIR%/}/s_tree.trees"
OUT_BASELINE="${SIMPHY_RUN_DIR%/}/out-${METHOD}.tre"
LOCK_FILE="${SIMPHY_RUN_DIR%/}/.${METHOD}.lock"

# Debug/tracing
if [[ "${DEBUG:-0}" = "1" ]]; then
  set -x
fi

# checkpoint: if lock file exists and --fresh not provided, skip everything
if [[ "$FRESH" = false && -f "${LOCK_FILE}" ]]; then
  echo "SKIPPING: ${LOCK_FILE} already exists. Use --fresh to force rerun."
  exit 0
fi

echo "Parameters:"
echo "  method:         $METHOD"
echo "  taxa_num:       $TAXA_NUM"
echo "  gene_trees:     $GENE_TREES"
echo "  replicate:      $REPLICATE"
if [[ "$SIMPHY_DATA_DIR_SET" = true ]]; then
  echo "  simphy data dir: $SIMPHY_DATA_DIR (custom)"
else
  echo "  simphy data dir: ${SIMPHY_DIR%/}/data (default)"
fi
echo "  simphy run dir: $SIMPHY_RUN_DIR"
echo "  out baseline:   $OUT_BASELINE"
echo "  stat file:      $STAT_FILE"
echo

if [[ ! -f "$ALL_GT_FILE" ]]; then
  echo "Error: gene-tree file not found at $ALL_GT_FILE"
  exit 6
fi

# Pick ASTER bin for large taxa if not specified
if [[ "$METHOD" == "aster" && -z "$ASTER_BIN" ]]; then
  if [[ "$TAXA_NUM" -ge 5000 ]]; then
    ASTER_BIN="bin/astral4_int128"
  fi
fi

mkdir -p "${SIMPHY_RUN_DIR%/}"

BASELINE_CMD=("${STELAR_ROOT%/}/run-baseline-with-monitor.sh" -m "$METHOD" -i "$ALL_GT_FILE" -o "$OUT_BASELINE" --stelar-root "$STELAR_ROOT")

if [[ "$TIME_MONITOR" = false ]]; then
  BASELINE_CMD+=(--no-time-monitor)
fi
if [[ "$GPU_MONITOR" = false ]]; then
  BASELINE_CMD+=(--no-gpu-monitor)
fi
if [[ "$NO_NOTIFY" = true ]]; then
  BASELINE_CMD+=(--no-notify)
fi
if [[ "$DEBUG" = 1 ]]; then
  BASELINE_CMD+=(--debug)
fi

if [[ -n "$ASTER_OPTS" ]]; then
  BASELINE_CMD+=(--aster-opts "$ASTER_OPTS")
fi
if [[ -n "$ASTER_BIN" ]]; then
  BASELINE_CMD+=(--aster-bin "$ASTER_BIN")
fi
if [[ -n "$ASTRAL_OPTS" ]]; then
  BASELINE_CMD+=(--astral-opts "$ASTRAL_OPTS")
fi
if [[ -n "$ASTRAL_XMS" ]]; then
  BASELINE_CMD+=(--astral-xms "$ASTRAL_XMS")
fi
if [[ -n "$ASTRAL_XMX" ]]; then
  BASELINE_CMD+=(--astral-xmx "$ASTRAL_XMX")
fi
if [[ -n "$TREEQMC_OPTS" ]]; then
  BASELINE_CMD+=(--treeqmc-opts "$TREEQMC_OPTS")
fi
if [[ -n "$WQFM_OPTS" ]]; then
  BASELINE_CMD+=(--wqfm-opts "$WQFM_OPTS")
fi
if [[ -n "$SUPERTRIPLETS_OPTS" ]]; then
  BASELINE_CMD+=(--supertriplets-opts "$SUPERTRIPLETS_OPTS")
fi
if [[ -n "$TMC_OPTS" ]]; then
  BASELINE_CMD+=(--tmc-opts "$TMC_OPTS")
fi

echo "==> Running baseline via: ${BASELINE_CMD[*]}"
set +e
"${BASELINE_CMD[@]}"
BASELINE_EXIT_CODE=$?
set -e

STATS_FROM_BASELINE="${OUT_BASELINE%.tre}_stats.csv"

if [[ ! -f "$STATS_FROM_BASELINE" ]]; then
  echo "Error: baseline stats file not found at $STATS_FROM_BASELINE"
  exit 7
fi

BASELINE_EXIT_CODE_FROM_STATS=$(awk -F, 'NR==2{print $7}' "$STATS_FROM_BASELINE")
RUNNING_TIME=$(awk -F, 'NR==2{print $4}' "$STATS_FROM_BASELINE")
MAX_CPU_MB=$(awk -F, 'NR==2{print $5}' "$STATS_FROM_BASELINE")
MAX_GPU_MB=$(awk -F, 'NR==2{print $6}' "$STATS_FROM_BASELINE")

if [[ "$BASELINE_EXIT_CODE" -ne 0 || "$BASELINE_EXIT_CODE_FROM_STATS" != "0" ]]; then
  echo "Baseline failed (exit code: ${BASELINE_EXIT_CODE}, stats exit: ${BASELINE_EXIT_CODE_FROM_STATS})"
  echo "Skipping stat/lock file creation."
  exit "${BASELINE_EXIT_CODE:-1}"
fi

# RF calculation
RF_RATE="NA"
if [[ -f "$OUT_BASELINE" && -f "$TRUE_SPECIES_TREE" ]]; then
  if [[ -f "${STELAR_ROOT%/}/rf.py" && -x "$(command -v python3)" ]]; then
    rf_output=$(python3 "${STELAR_ROOT%/}/rf.py" "$OUT_BASELINE" "$TRUE_SPECIES_TREE" 2>&1) || true
    rf_candidate=$(echo "$rf_output" | grep -Eo '[0-9]+(\.[0-9]+)?' | head -n1 || true)
    if [[ -n "$rf_candidate" ]]; then
      RF_RATE="$rf_candidate"
    fi
  fi
fi

# Write CSV (overwrite every run)
mkdir -p "$(dirname "$STAT_FILE")"
echo "alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb" > "$STAT_FILE"
CSV_ROW="${METHOD},${TAXA_NUM},${GENE_TREES},${REPLICATE},${SB},${SPMIN},${SPMAX},${RF_RATE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB}"
echo "$CSV_ROW" >> "$STAT_FILE"

echo "Wrote stats to $STAT_FILE"

# Create lock file to indicate successful completion
touch "${LOCK_FILE}"
echo -e "\\033[1;32m✅ Run completed successfully. Lock file created.\\033[0m"

# Send notification (ntfy) - only on success
if [[ "$NO_NOTIFY" = false ]] && command -v curl >/dev/null 2>&1; then
  curl -s -d "🎉 ${METHOD^^} completed for ${TAXA_NUM} taxa and ${GENE_TREES} gene trees!

⚙️ Config: sb=${SB} spmin=${SPMIN} spmax=${SPMAX} rep=${REPLICATE}

📊 Results:
RF: ${RF_RATE} | Time: ${RUNNING_TIME}s
CPU: ${MAX_CPU_MB} MB | GPU: ${MAX_GPU_MB} MB

📁 ${STAT_FILE}

📋 CSV Row:
${CSV_ROW}" https://ntfy.sh/anik-test || true
fi

echo "Done."
exit 0
