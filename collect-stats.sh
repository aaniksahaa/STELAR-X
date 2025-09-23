#!/usr/bin/env bash
# collect-stelar-stats.sh (simplified)
# Merges stat-stelar.csv files under simphy/data into one perf-stelar.csv
# and appends gt-gt,gt-st from a stat-sim.csv in the same directory (if present).
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
verifies headers match, merges them into --out, and appends gt-gt,gt-st from a
stat-sim.csv located in the same directory as each stat-stelar.csv (if present).

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
  # read a line from stdin
  local L
  IFS= read -r L || return 1
  # remove CR, trim leading/trailing whitespace
  L="${L//$'\r'/}"
  # trim
  L="$(echo "$L" | awk '{$1=$1;print}')"
  printf "%s" "$L"
}

# canonical header from first file
FIRST="${stelar_files[0]}"
HEADER=$(head -n 1 "$FIRST" | tr -d '\r' | awk '{$1=$1;print}')
EXT_HEADER="${HEADER},gt-gt,gt-st"

# write header
printf "%s\n" "$EXT_HEADER" > "$OUT_FILE"

merged_rows=0
skipped_files=0

for stelar in "${stelar_files[@]}"; do
  if [[ ! -f "$stelar" ]]; then
    echo "Warning: file disappeared: $stelar" >&2
    ((skipped_files++))
    continue
  fi

  file_header=$(head -n 1 "$stelar" | tr -d '\r' | awk '{$1=$1;print}')
  if [[ "$file_header" != "$HEADER" ]]; then
    echo "Warning: header mismatch, skipping file: $stelar" >&2
    echo "  expected: $HEADER" >&2
    echo "  found:    $file_header" >&2
    ((skipped_files++))
    continue
  fi

  dir=$(dirname "$stelar")
  stat_sim="${dir%/}/stat-sim.csv"

  gt_gt=""
  gt_st=""

  if [[ -f "$stat_sim" ]]; then
    # Extract the second non-empty line (data line) and get 7th and 8th comma fields.
    # This assumes stat-sim has header then single data row.
    # We remove possible CRs before processing.
    data_line=$(awk 'BEGIN{FS=OFS=","} NR>0 {gsub(/\r/,""); if ($0 ~ /^[[:space:]]*$/) next; lines[++c]=$0} END{ if (c>=2) print lines[2]; else if (c==1) print lines[1] }' "$stat_sim" || true)
    if [[ -n "$data_line" ]]; then
      # get 7th and 8th fields (gt-gt, gt-st). Use awk field extraction to be robust for quotes etc.
      gt_gt=$(echo "$data_line" | awk -F, '{print $7}' | tr -d '\r' || true)
      gt_st=$(echo "$data_line" | awk -F, '{print $8}' | tr -d '\r' || true)
      # trim whitespace
      gt_gt="$(echo "$gt_gt" | awk '{$1=$1;print}')"
      gt_st="$(echo "$gt_st" | awk '{$1=$1;print}')"
    fi
  fi

  # Append all data rows from stat-stelar skipping header, and append gt fields
  # We use awk to print original line and then appended gt fields (preserve commas)
  # If gt_gt/gt_st are empty, we still append the commas (so columns align).
  if [[ -z "$gt_gt" && -z "$gt_st" ]]; then
    # append empty values
    awk 'NR>1{print $0 ",,"}' "$stelar" >> "$OUT_FILE"
    appended=$(awk 'END{print NR-1}' "$stelar")
  else
    # escape any embedded newlines are unlikely; just append values
    # ensure fields do not contain newlines
    gt_gt_safe=$(printf "%s" "$gt_gt" | tr -d '\n' | tr -d '\r')
    gt_st_safe=$(printf "%s" "$gt_st" | tr -d '\n' | tr -d '\r')
    awk -v a="$gt_gt_safe" -v b="$gt_st_safe" 'NR>1{print $0 "," a "," b}' "$stelar" >> "$OUT_FILE"
    appended=$(awk 'END{print NR-1}' "$stelar")
  fi

  appended=$(echo "$appended" | tr -d '[:space:]')
  merged_rows=$((merged_rows + appended))
done

echo "Merged $merged_rows data rows into $OUT_FILE"
if [[ $skipped_files -gt 0 ]]; then
  echo "Skipped $skipped_files files due to header mismatch or missing file."
fi

echo "Done."
