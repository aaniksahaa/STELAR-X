#!/usr/bin/env bash
# collect-stelar-stats.sh (simplified)
# Merges stat-stelar.csv files under simphy/data into one perf-stelar.csv
# and appends gt-gt,gt-st from a stat-sim.csv in the same directory (if present).
# Handles header mismatches by taking the union of columns, filling missing ones with empty values.
# Outputs columns in the prescribed order: alg,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,running-time-s,max-cpu-mb,max-gpu-mb,gt-gt,gt-st
#
# Usage:
#   ./collect-stelar-stats.sh
#   ./collect-stelar-stats.sh --base-dir /home/you/research --out /tmp/perf-stelar.csv
#   ./collect-stelar-stats.sh --simphy-dir /home/you/research/STELAR-MP/simphy

set -euo pipefail

BASE_DIR="${HOME}/phylogeny"
SIMPHY_DIR=""
OUT_FILE="./perf-stelar.csv"

print_help() {
  cat <<EOF
collect-stelar-stats.sh

Finds all stat-stelar.csv under <SIMPHY_DIR>/data (default derived from --base-dir),
takes the union of headers, merges them into --out with columns in prescribed order,
and appends gt-gt,gt-st from a stat-sim.csv located in the same directory as each stat-stelar.csv (if present).

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

# find stat-stelar.csv files
mapfile -t stelar_files < <(find "$SIMPHY_DATA_DIR" -type f -name 'stat-stelar.csv' -print 2>/dev/null | sort)

if [[ ${#stelar_files[@]} -eq 0 ]]; then
  echo "No stat-stelar.csv files found under $SIMPHY_DATA_DIR"
  exit 0
fi

echo "Found ${#stelar_files[@]} stat-stelar.csv files. Merging..."

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

for stelar in "${stelar_files[@]}"; do
  if [[ ! -f "$stelar" ]]; then
    echo "Warning: file disappeared: $stelar" >&2
    ((skipped_files++))
    continue
  fi

  file_header=$(head -n 1 "$stelar" | tr -d '\r' | awk '{$1=$1;print}')
  IFS=',' read -ra file_cols <<< "$file_header"

  # Create mapping of file columns to indices
  declare -A col_map
  for i in "${!file_cols[@]}"; do
    col_map["${file_cols[$i]}"]=$i
  done

  # Process stat-sim.csv for gt-gt, gt-st
  dir=$(dirname "$stelar")
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
  }' "$stelar" >> "$OUT_FILE"

  appended=$(awk 'END{print NR-1}' "$stelar")
  appended=$(echo "$appended" | tr -d '[:space:]')
  merged_rows=$((merged_rows + appended))
done

echo "Merged $merged_rows data rows into $OUT_FILE"
if [[ $skipped_files -gt 0 ]]; then
  echo "Skipped $skipped_files files due to missing file."
fi

echo "Done."