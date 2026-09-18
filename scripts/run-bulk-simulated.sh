#!/usr/bin/env bash
# run-bulk-simulated.sh
#
# Runs sim.sh and test-stelarx-simulated.sh or test-baseline-simulated.sh
# over all combinations of parameter lists.
#
# Usage:
#   ./scripts/run-bulk-simulated.sh -m stelar
#   ./scripts/run-bulk-simulated.sh --project-root /path/to/checkout

set -euo pipefail

STELARX_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${STELARX_ROOT}/scripts/phylogeny-data-dir.sh"
source "${STELARX_ROOT}/scripts/simphy-outputs-dir.sh"
source "${STELARX_ROOT}/scripts/experiment-setting-name.sh"

BASE_DIR=""
BASE_DIR_PROVIDED=false
METHOD="stelarx"  # default method
FRESH=false
NO_NOTIFY=false
GPU_MONITOR=true
SIMPHY_DATA_DIR=""
SIMPHY_OUTPUTS_DIR=""
OUTPUTS_MIRROR=true
NUM_REPLICATES=1
ASSUME_YES=false
DRY_RUN=false

T_LIST=(10)

# T_LIST=(1000 2500 5000 7500 10000)

G_LIST=(10)
SB_LIST=(0.000001)
SPMIN_LIST=(100000)
SPMAX_LIST=(200000)

# Optional list of per-replicate runs that have ALREADY been completed in an
# earlier invocation (for example on another machine or in a previous sweep),
# so that this sweep does not spend compute re-running them. This is purely a
# scheduling convenience: results for these replicates are produced exactly
# like every other run, just not again here. Each entry is an exact
# TAXA,GENE_TREES,SB,SPMIN,SPMAX,REPLICATE tuple with replicate names in R<n>
# form. sim.sh still prepares the surrounding dataset batch; only the listed
# per-replicate inference run is not repeated. Leave empty to run everything.
ALREADY_COMPLETED_SIMULATED_CONFIGS=(
)

IS_SIMULATED_CONFIG_ALREADY_COMPLETED() {
  local CANDIDATE_CONFIG="$1,$2,$3,$4,$5,$6"
  local COMPLETED_CONFIG
  for COMPLETED_CONFIG in "${ALREADY_COMPLETED_SIMULATED_CONFIGS[@]}"; do
    [[ "$CANDIDATE_CONFIG" == "$COMPLETED_CONFIG" ]] && return 0
  done
  return 1
}

# Method-specific options (passed through)
ASTER_OPTS=""
ASTER_BIN=""
ASTRAL_OPTS=""
STELARX_OPTS_LIST_RAW=""
ASTRAL_XMS=""
ASTRAL_XMX=""
TREEQMC_OPTS=""
WQFM_OPTS=""
SUPERTRIPLETS_OPTS=""
TMC_OPTS=""

# Permit the already-completed predicate and uppercase configuration array to be
# loaded by the isolated regression test without executing a simulated-data sweep.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  return 0
fi

print_help() {
  cat <<EOF
run-bulk-simulated.sh

Runs sim.sh and test-stelarx-simulated.sh or test-baseline-simulated.sh for all combinations of parameter lists.

Options:
  --method, -m      Method to use: stelarx (default: stelarx)
  --project-root    STELAR-X checkout root (default: this script's directory)
  --base-dir, -b    Compatibility alias for --project-root
  --num-replicates, -n  Number of replicates to run (default: 1)
  --taxa-list LIST       Comma/space-separated taxon counts (default: 10)
  --genes-list LIST      Comma/space-separated gene-tree counts (default: 10)
  --sb-list LIST         Comma/space-separated speciation rates
  --spmin-list LIST      Comma/space-separated minimum population sizes
  --spmax-list LIST      Comma/space-separated maximum population sizes
  --simphy-data-dir DIR  Store/read generated datasets under DIR
                         (default: \$PHYLOGENY_DATA_DIR/simphy/data)
  --simphy-outputs-dir DIR
                         Reproducibility mirror for run outputs and SimPhy commands
                         (default: \$PHYLOGENY_DATA_DIR/outputs/simphy)
  --no-outputs-mirror    Do not mirror run outputs
  --fresh           Pass --fresh to sim.sh and test scripts (recreate outputs)
  --no-gpu-monitor  Disable GPU-memory sampling
  --no-notify, -nn  Disable completion notifications
  --yes, -y         Start without the interactive confirmation
  --dry-run         Print the run plan (dataset / replicate / setting) and exit
  --opts, --alg-opts       Extra options for one STELAR-X simulated setting
  --opts-list, --alg-opts-list
                         Semicolon-separated list of STELAR-X option strings to loop over
  --help, -h        Show this message

Examples:
  ./scripts/run-bulk-simulated.sh --opts "--search-space S2 -vv"
  ./scripts/run-bulk-simulated.sh --taxa-list "10,20" --genes-list "10,50" --num-replicates 3
  ./scripts/run-bulk-simulated.sh --opts-list "--search-space S1 -vv;--search-space S2 -vv;--search-space S3 -vv"
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --method|-m) METHOD="$2"; shift 2 ;;
    --project-root|--base-dir|-b) BASE_DIR="$2"; BASE_DIR_PROVIDED=true; shift 2 ;;
    --num-replicates|-n) NUM_REPLICATES="$2"; shift 2 ;;
    --taxa-list) read -r -a T_LIST <<< "${2//,/ }"; shift 2 ;;
    --genes-list) read -r -a G_LIST <<< "${2//,/ }"; shift 2 ;;
    --sb-list) read -r -a SB_LIST <<< "${2//,/ }"; shift 2 ;;
    --spmin-list) read -r -a SPMIN_LIST <<< "${2//,/ }"; shift 2 ;;
    --spmax-list) read -r -a SPMAX_LIST <<< "${2//,/ }"; shift 2 ;;
    --simphy-data-dir) SIMPHY_DATA_DIR="$2"; shift 2 ;;
    --simphy-outputs-dir) SIMPHY_OUTPUTS_DIR="$2"; shift 2 ;;
    --no-outputs-mirror) OUTPUTS_MIRROR=false; shift ;;
    --opts|--alg-opts|--stelarx-opts) ASTRAL_OPTS="$2"; shift 2 ;;
    --opts=*|--alg-opts=*|--stelarx-opts=*) ASTRAL_OPTS="${1#*=}"; shift ;;
    --opts-list|--alg-opts-list|--stelarx-opts-list) STELARX_OPTS_LIST_RAW="$2"; shift 2 ;;
    --opts-list=*|--alg-opts-list=*|--stelarx-opts-list=*) STELARX_OPTS_LIST_RAW="${1#*=}"; shift ;;
    --fresh) FRESH=true; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --yes|-y) ASSUME_YES=true; shift ;;
    --dry-run|--plan-only) DRY_RUN=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Validate method
case "$METHOD" in
  stelarx|stelar) METHOD="stelarx" ;;
  *)
    echo "Error: --method must be stelarx."
    exit 1
    ;;
esac

# Keep the defaults deliberately small. Larger experiment matrices must be
# requested explicitly through the list options above.
for parameter_list in T_LIST G_LIST SB_LIST SPMIN_LIST SPMAX_LIST; do
  declare -n values="$parameter_list"
  if [[ ${#values[@]} -eq 0 ]]; then
    echo "Error: $parameter_list cannot be empty."
    exit 1
  fi
done
unset -n values

if [[ ! "$NUM_REPLICATES" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: --num-replicates must be a positive integer."
  exit 1
fi

SIMPHY_DATA_DIR="$(stelarx_prepare_simphy_data_dir "$SIMPHY_DATA_DIR")"

# -------------------------------
# execution
# -------------------------------

# Build base-dir argument if provided
BASE_DIR_ARGS=()
if $BASE_DIR_PROVIDED; then
  BASE_DIR_ARGS=(--project-root "$BASE_DIR")
  echo "Project root: $BASE_DIR"
else
  echo "Project root: $STELARX_ROOT"
fi

# Build fresh argument if provided
FRESH_ARGS=()
if $FRESH; then
  FRESH_ARGS=(--fresh)
  echo "Fresh:    yes"
else
  echo "Fresh:    no"
fi
SIM_DATA_ARGS=(--simphy-data-dir "$SIMPHY_DATA_DIR")
SHARED_TEST_ARGS=("${SIM_DATA_ARGS[@]}")
if [[ "$OUTPUTS_MIRROR" == false ]]; then
  SHARED_TEST_ARGS+=(--no-outputs-mirror)
  echo "Outputs mirror: disabled"
else
  if [[ -n "$SIMPHY_OUTPUTS_DIR" ]]; then
    SHARED_TEST_ARGS+=(--simphy-outputs-dir "$SIMPHY_OUTPUTS_DIR")
    echo "Outputs mirror: $SIMPHY_OUTPUTS_DIR"
  else
    echo "Outputs mirror: $(stelarx_default_simphy_outputs_dir "$SIMPHY_DATA_DIR")"
  fi
fi
if [[ "$GPU_MONITOR" == false ]]; then
  SHARED_TEST_ARGS+=(--no-gpu-monitor)
fi
if [[ "$NO_NOTIFY" == true ]]; then
  SHARED_TEST_ARGS+=(--no-notify)
fi
echo "Method:   $METHOD"
echo "Replicates: $NUM_REPLICATES"
echo "SimPhy data: $SIMPHY_DATA_DIR"

STELARX_OPTS_LIST=()
if [[ -n "$STELARX_OPTS_LIST_RAW" ]]; then
  IFS=';' read -r -a raw_opts_list <<< "$STELARX_OPTS_LIST_RAW"
  for opts in "${raw_opts_list[@]}"; do
    opts="$(echo "$opts" | sed 's/^ *//;s/ *$//')"
    [[ -n "$opts" ]] && STELARX_OPTS_LIST+=("$opts")
  done
fi
if [[ ${#STELARX_OPTS_LIST[@]} -eq 0 ]]; then
  STELARX_OPTS_LIST+=("${ASTRAL_OPTS}")
fi

# -------------------------------
# plan: one line per dataset / replicate / setting, then one confirmation
# -------------------------------
declare -a PLAN_LINES=()
declare -a PLAN_DATASETS=()   # "t g sb spmin spmax" per dataset, in run order
planned_runs=0
already_completed_runs=0

# Collapse sorted replicate numbers into "R1-R4, R6" style ranges.
format_replicate_ranges() {
  local -a nums=("$@")
  local out="" start="" prev=""
  local n
  for n in "${nums[@]}"; do
    if [[ -z "$start" ]]; then
      start=$n; prev=$n; continue
    fi
    if (( n == prev + 1 )); then
      prev=$n; continue
    fi
    out+="${out:+, }R${start}"; (( start != prev )) && out+="-R${prev}"
    start=$n; prev=$n
  done
  if [[ -n "$start" ]]; then
    out+="${out:+, }R${start}"; (( start != prev )) && out+="-R${prev}"
  fi
  printf '%s' "$out"
}

for t in "${T_LIST[@]}"; do
  for g in "${G_LIST[@]}"; do
    for sb in "${SB_LIST[@]}"; do
      for spmin in "${SPMIN_LIST[@]}"; do
        for spmax in "${SPMAX_LIST[@]}"; do
          DATASET_NAME="t_${t}_g_${g}_sb_${sb}_spmin_${spmin}_spmax_${spmax}"
          PLAN_DATASETS+=("$t $g $sb $spmin $spmax")
          included=()
          already_completed=()
          for ((i=1; i<=NUM_REPLICATES; i++)); do
            if IS_SIMULATED_CONFIG_ALREADY_COMPLETED "$t" "$g" "$sb" "$spmin" "$spmax" "R$i"; then
              already_completed+=("$i")
            else
              included+=("$i")
            fi
          done
          for STELARX_OPTS_ITEM in "${STELARX_OPTS_LIST[@]}"; do
            SETTING_NAME="$(build_setting_name_from_opts "$STELARX_OPTS_ITEM")"
            line="  ${DATASET_NAME} / "
            if [[ ${#included[@]} -gt 0 ]]; then
              line+="$(format_replicate_ranges "${included[@]}")"
            else
              line+="(none)"
            fi
            line+=" / ${SETTING_NAME}"
            if [[ ${#already_completed[@]} -gt 0 ]]; then
              line+="   (already completed earlier, not re-run: $(format_replicate_ranges "${already_completed[@]}"))"
            fi
            PLAN_LINES+=("$line")
            planned_runs=$((planned_runs + ${#included[@]}))
            already_completed_runs=$((already_completed_runs + ${#already_completed[@]}))
          done
        done
      done
    done
  done
done

echo
echo "Run plan (${#PLAN_DATASETS[@]} dataset(s), ${NUM_REPLICATES} replicate(s), ${#STELARX_OPTS_LIST[@]} setting(s)):"
echo "  <dataset> / <replicates> / <setting folder under ${METHOD}_outputs>"
printf '%s\n' "${PLAN_LINES[@]}"
echo
echo "Total: ${planned_runs} run(s) to execute, ${already_completed_runs} already completed earlier (not re-run)."
echo "Completed runs are skipped unless --fresh is given."

if [[ "$DRY_RUN" == true ]]; then
  echo "Dry run; nothing was simulated or executed."
  exit 0
fi

if [[ "$ASSUME_YES" == false ]]; then
  if [[ -t 0 ]]; then
    read -r -p "Proceed with ${planned_runs} run(s)? [y/N]: " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
      echo "Cancelled; nothing was executed."
      exit 0
    fi
  else
    echo "Non-interactive session: proceeding without confirmation (pass --yes to silence this note)."
  fi
fi

echo
echo "Starting bulk runs..."

for DATASET_SPEC in "${PLAN_DATASETS[@]}"; do
  read -r t g sb spmin spmax <<< "$DATASET_SPEC"

  echo ">>> Running: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax (method=$METHOD)"

  ./scripts/sim.sh -rs "$NUM_REPLICATES" "${BASE_DIR_ARGS[@]}" "${SIM_DATA_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}"

  # Run replicates
  for ((i=1; i<=NUM_REPLICATES; i++)); do
    REPLICATE_NAME="R$i"
    if IS_SIMULATED_CONFIG_ALREADY_COMPLETED \
        "$t" "$g" "$sb" "$spmin" "$spmax" "$REPLICATE_NAME"; then
      echo "  SKIPPING already-completed run (result exists from an earlier sweep): t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax replicate=$REPLICATE_NAME"
      continue
    fi
    echo "  Running replicate $REPLICATE_NAME with $METHOD"

    for STELARX_OPTS_ITEM in "${STELARX_OPTS_LIST[@]}"; do
      echo "  >>> t_${t}_g_${g}_sb_${sb}_spmin_${spmin}_spmax_${spmax} / ${REPLICATE_NAME} / $(build_setting_name_from_opts "$STELARX_OPTS_ITEM")"
      TEST_CMD=("${STELARX_ROOT}/scripts/test-stelarx-simulated.sh" -r "$REPLICATE_NAME" "${BASE_DIR_ARGS[@]}" "${SHARED_TEST_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}")
      if [[ -n "$STELARX_OPTS_ITEM" ]]; then
        TEST_CMD+=(--opts "$STELARX_OPTS_ITEM")
      fi
      "${TEST_CMD[@]}"
    done
  done
done

echo "All runs finished."














# T_LIST=(1000)
# G_LIST=(1000 2500 5000 7500 10000)
# SB_LIST=(0.000001)
# SPMIN_LIST=(100000)
# SPMAX_LIST=(200000)












# echo "Starting bulk runs... phase 2"

# for t in "${T_LIST[@]}"; do
#   for g in "${G_LIST[@]}"; do
#     for sb in "${SB_LIST[@]}"; do
#       for spmin in "${SPMIN_LIST[@]}"; do
#         for spmax in "${SPMAX_LIST[@]}"; do

#           echo ">>> Running: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax (method=$METHOD)"
          
#           ./scripts/sim.sh -rs "$NUM_REPLICATES" "${BASE_DIR_ARGS[@]}" "${SIM_DATA_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}"
          
#           # Run replicates
#           for ((i=1; i<=NUM_REPLICATES; i++)); do
#             REPLICATE_NAME="R$i"
#             if IS_SIMULATED_CONFIG_ALREADY_COMPLETED \
#                 "$t" "$g" "$sb" "$spmin" "$spmax" "$REPLICATE_NAME"; then
#               echo "  SKIPPING already-completed run (result exists from an earlier sweep): t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax replicate=$REPLICATE_NAME"
#               continue
#             fi
#             echo "  Running replicate $REPLICATE_NAME with $METHOD"
            
#             for STELARX_OPTS_ITEM in "${STELARX_OPTS_LIST[@]}"; do
#               TEST_CMD=("${STELARX_ROOT}/scripts/test-stelarx-simulated.sh" -r "$REPLICATE_NAME" "${BASE_DIR_ARGS[@]}" "${SHARED_TEST_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}")
#               if [[ -n "$STELARX_OPTS_ITEM" ]]; then
#                 TEST_CMD+=(--opts "$STELARX_OPTS_ITEM")
#               fi
#               "${TEST_CMD[@]}"
#             done
#           done

#         done
#       done
#     done
#   done
# done

# echo "All runs finished."
