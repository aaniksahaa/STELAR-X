#!/usr/bin/env bash
# collect-stats-simulated.sh (multi-algorithm version)
# Merges stat-*.csv files under simphy/data into one combined CSV file
# and appends gt-gt,gt-st from a stat-sim.csv in the same directory (if present).
# Handles header mismatches by taking the union of columns, filling missing ones with empty values.
# Outputs columns in the prescribed order: alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb,gt-gt,gt-st
#
# Usage:
#   ./collect-stats-simulated.sh
#   ./collect-stats-simulated.sh --base-dir /home/you/research --out /tmp/perf-combined.csv
#   ./collect-stats-simulated.sh --simphy-dir /home/you/research/STELAR-MP/simphy

set -euo pipefail

# Algorithm configuration - modify this to select which algorithms to collect
ALGORITHMS=("stelar" "astral" "tree-qmc" "wqfm-tree")

BASE_DIR="${HOME}/phylogeny"
SIMPHY_DIR=""
OUT_FILE="./perf-combined.csv"

print_help() {
  cat <<EOF
collect-stats-simulated.sh

Finds all stat-*.csv files for configured algorithms under <SIMPHY_DIR>/data (default derived from --base-dir),
takes the union of headers, merges them into --out with columns in prescribed order,
and appends gt-gt,gt-st from a stat-sim.csv located in the same directory as each stat file (if present).

Configured algorithms: ${ALGORITHMS[*]}

Options:
  --base-dir, -b    Base directory (default: ${BASE_DIR})
  --simphy-dir      Path to simphy dir (overrides --base-dir)
  --out, -o         Output CSV path (default: ${OUT_FILE})
  --help, -h        Show this message
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --simphy-dir) SIMPHY_DIR="$2"; shift 2 ;;
    --out|-o) OUT_FILE="$2"; shift 2 ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$SIMPHY_DIR" ]]; then
  SIMPHY_DIR="${BASE_DIR%/}/STELAR-MP/simphy"
fi

SIMPHY_DATA_DIR="${SIMPHY_DIR%/}/data"

if [[ ! -d "$SIMPHY_DATA_DIR" ]]; then
  echo "Error: simphy data directory not found at: $SIMPHY_DATA_DIR" >&2
  exit 2
fi

# find stat files for all configured algorithms
declare -a all_stat_files
total_files=0

for alg in "${ALGORITHMS[@]}"; do
  mapfile -t alg_files < <(find "$SIMPHY_DATA_DIR" -type f -name "stat-${alg}.csv" -print 2>/dev/null | sort)
  if [[ ${#alg_files[@]} -gt 0 ]]; then
    echo "Found ${#alg_files[@]} stat-${alg}.csv files"
    all_stat_files+=("${alg_files[@]}")
    total_files=$((total_files + ${#alg_files[@]}))
  else
    echo "No stat-${alg}.csv files found"
  fi
done

if [[ ${#all_stat_files[@]} -eq 0 ]]; then
  echo "No stat files found for any configured algorithms under $SIMPHY_DATA_DIR"
  echo "Configured algorithms: ${ALGORITHMS[*]}"
  exit 0
fi

echo "Found total ${total_files} stat files across all algorithms. Merging..."

# helper: normalize a single line (remove CR and trim)
norm_line() {
  local L
  IFS= read -r L || return 1
  L="${L//$'\r'/}"
  L="$(echo "$L" | awk '{$1=$1;print}')"
  printf "%s" "$L"
}

# Define the prescribed header order
PRESCRIBED_HEADER="alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb,gt-gt,gt-st"

# Write prescribed header to output
printf "%s\n" "$PRESCRIBED_HEADER" > "$OUT_FILE"

merged_rows=0
skipped_files=0

for stat_file in "${all_stat_files[@]}"; do
  if [[ ! -f "$stat_file" ]]; then
    echo "Warning: file disappeared: $stat_file" >&2
    ((skipped_files++))
    continue
  fi

  file_header=$(head -n 1 "$stat_file" | tr -d '\r' | awk '{$1=$1;print}')
  IFS=',' read -ra file_cols <<< "$file_header"

  # Create mapping of file columns to indices
  declare -A col_map
  for i in "${!file_cols[@]}"; do
    col_map["${file_cols[$i]}"]=$i
  done

  # Process stat-sim.csv for gt-gt, gt-st
  dir=$(dirname "$stat_file")
  stat_sim="${dir%/}/stat-sim.csv"
  gt_gt=""
  gt_st=""

  if [[ -f "$stat_sim" ]]; then
    data_line=$(awk 'BEGIN{FS=OFS=","} NR>0 {gsub(/\r/,""); if ($0 ~ /^[[:space:]]*$/) next; lines[++c]=$0} END{ if (c>=2) print lines[2]; else if (c==1) print lines[1] }' "$stat_sim" || true)
    if [[ -n "$data_line" ]]; then
      gt_gt=$(echo "$data_line" | awk -F, '{print $7}' | tr -d '\r' || true)
      gt_st=$(echo "$data_line" | awk -F, '{print $8}' | tr -d '\r' || true)
      gt_gt="$(echo "$gt_gt" | awk '{$1=$1;print}')"
      gt_st="$(echo "$gt_st" | awk '{$1=$1;print}')"
    fi
  fi

  # Process data rows, mapping to prescribed header
  awk -F, -v OFS=',' -v header="$file_header" -v prescribed="$PRESCRIBED_HEADER" -v gt_gt="$gt_gt" -v gt_st="$gt_st" '
  BEGIN {
    split(header, h, ",");
    split(prescribed, u, ",");
    for (i in h) col_map[h[i]]=i;
  }
  NR>1 {
    # Initialize output array with empty strings
    for (i=1; i<=length(u); i++) out[i]="";
    # Copy fields that exist in the file
    for (i in h) {
      for (j=1; j<=length(u); j++) {
        if (u[j] == h[i] && u[j] != "gt-gt" && u[j] != "gt-st") {
          out[j]=$i;
          break;
        }
      }
    }
    # Set gt-gt and gt-st
    out[length(u)-1]=gt_gt;
    out[length(u)]=gt_st;
    # Build output line
    line="";
    for (i=1; i<=length(u); i++) {
      line=(line ? line "," : "") out[i];
    }
    print line;
  }' "$stat_file" >> "$OUT_FILE"

  appended=$(awk 'END{print NR-1}' "$stat_file")
  appended=$(echo "$appended" | tr -d '[:space:]')
  merged_rows=$((merged_rows + appended))
done

echo "Merged $merged_rows data rows into $OUT_FILE"
if [[ $skipped_files -gt 0 ]]; then
  echo "Skipped $skipped_files files due to missing file."
fi

echo "Done."