#!/usr/bin/env bash
# collect-stelar-stats.sh
# Collect all stat-stelar.csv files under simphy/data and merge into one perf-stelar.csv
#
# Usage:
#   ./collect-stelar-stats.sh
#   ./collect-stelar-stats.sh --base-dir /home/you/research --out /tmp/perf-stelar.csv
#   ./collect-stelar-stats.sh --simphy-dir /home/you/research/STELAR-MP/simphy
#
# Defaults:
#   base-dir = ${HOME}/phylogeny
#   simphy-dir = ${base-dir}/STELAR-MP/simphy
#   output = ./perf-stelar.csv

set -euo pipefail

BASE_DIR="${HOME}/phylogeny"
SIMPHY_DIR=""
OUT_FILE="./perf-stelar.csv"

print_help() {
  cat <<EOF
collect-stelar-stats.sh

Finds all stat-stelar.csv under <SIMPHY_DIR>/data (default derived from --base-dir),
verifies headers match, and writes a merged CSV to --out (overwrite).

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

# derive simphy dir if not provided
if [[ -z "$SIMPHY_DIR" ]]; then
  SIMPHY_DIR="${BASE_DIR%/}/STELAR-MP/simphy"
fi

SIMPHY_DATA_DIR="${SIMPHY_DIR%/}/data"

if [[ ! -d "$SIMPHY_DATA_DIR" ]]; then
  echo "Error: simphy data directory not found at: $SIMPHY_DATA_DIR" >&2
  exit 2
fi

# find stat-stelar.csv files
mapfile -t files < <(find "$SIMPHY_DATA_DIR" -type f -name 'stat-stelar.csv' -print 2>/dev/null | sort)

if [[ ${#files[@]} -eq 0 ]]; then
  echo "No stat-stelar.csv files found under $SIMPHY_DATA_DIR"
  exit 0
fi

echo "Found ${#files[@]} stat-stelar.csv files. Merging..."

# helper to normalize header: trim CR, leading/trailing whitespace
norm_header() {
  local f="$1"
  # read first non-empty line (protect in case of blank lines)
  head -n 1 "$f" | tr -d '\r' | awk '{$1=$1; print}'
}

FIRST="${files[0]}"
HEADER=$(norm_header "$FIRST")

# write header to output (overwrite)
printf "%s\n" "$HEADER" > "$OUT_FILE"

merged_rows=0
skipped_files=0

for f in "${files[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "Warning: file disappeared: $f" >&2
    ((skipped_files++))
    continue
  fi

  file_header=$(norm_header "$f")
  if [[ "$file_header" != "$HEADER" ]]; then
    echo "Warning: header mismatch, skipping file: $f" >&2
    echo "  expected: $HEADER" >&2
    echo "  found:    $file_header" >&2
    ((skipped_files++))
    continue
  fi

  # append data rows (skip header line)
  # count lines appended
  # use awk to skip exactly the first line and print rest (handles single-line files gracefully)
  appended=$(awk 'NR>1{print}' "$f" | tee -a "$OUT_FILE" | wc -l)
  appended=$(echo "$appended" | tr -d '[:space:]')
  merged_rows=$((merged_rows + appended))
done

echo "Merged $merged_rows data rows into $OUT_FILE"
if [[ $skipped_files -gt 0 ]]; then
  echo "Skipped $skipped_files files due to header mismatch or missing file."
fi

echo "Done."
