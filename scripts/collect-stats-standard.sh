#!/usr/bin/env bash
# Stat Collection Script for Standard Dataset Processing
# Collects all stat-*.csv files from dataset directories and merges them
# Usage: ./collect-stats-standard.sh [--base-dir /path/to/base] [--dataset-dir /path/to/datasets]
#   --base-dir, -b    Optional base directory (defaults to value below)
#   --dataset-dir, -d Optional dataset directory (defaults to BASE_DIR/datasets/standard)

set -euo pipefail

# =============================================================================
# DEFAULTS (edit these if you want different defaults)
# =============================================================================

BASE_DIR="${HOME}/phylogeny"

DATASET_DIR=""                        # dataset directory; defaults to BASE_DIR/datasets/standard
OUTPUT_FILE="stat-standard.csv"       # output merged CSV file

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# =============================================================================
# Helper functions
# =============================================================================

print_header() {
    echo "============================================="
    echo "=== Stat Collection Script (Standard) ==="
    echo "============================================="
    echo "BASE_DIR: $BASE_DIR"
    echo "DATASET_DIR: $DATASET_DIR"
    echo "OUTPUT_FILE: $OUTPUT_FILE"
    echo "============================================="
    echo
}

show_usage() {
    cat <<EOF
Usage: $0 [--base-dir /path/to/base] [--dataset-dir /path/to/datasets] [--output /path/to/output.csv]

Collects all stat-*.csv files from dataset directories and merges them.

--base-dir, -b     Base directory (overrides default)
--dataset-dir, -d  Dataset directory (overrides default BASE_DIR/datasets/standard)
--output, -o       Output CSV file (default: stat-standard.csv)
--help, -h         Show this help

The script will:
1. Find all stat-*.csv files recursively in the dataset directory
2. Check that all files have the same header format
3. Merge all data rows into a single CSV file
EOF
}

# =============================================================================
# Argument parsing
# =============================================================================

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
    --output|-o)
      OUTPUT_FILE="$2"
      shift 2
      ;;
    --help|-h)
      show_usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      show_usage
      exit 1
      ;;
  esac
done

# derive DATASET_DIR from BASE_DIR if not set
if [[ -z "${DATASET_DIR}" ]]; then
  DATASET_DIR="${BASE_DIR%/}/datasets/standard"
fi

print_header

# =============================================================================
# Validation
# =============================================================================

echo -e "${YELLOW}Performing validation checks...${NC}"

if [ ! -d "$BASE_DIR" ]; then
    echo -e "${RED}Error: BASE_DIR '$BASE_DIR' does not exist.${NC}"
    exit 1
fi

if [ ! -d "$DATASET_DIR" ]; then
    echo -e "${RED}Error: DATASET_DIR '$DATASET_DIR' does not exist.${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Validation checks passed${NC}"
echo

# =============================================================================
# Main collection logic
# =============================================================================

echo -e "${YELLOW}Searching for stat-*.csv files in $DATASET_DIR...${NC}"

# Find all stat-*.csv files
stat_files=($(find "$DATASET_DIR" -name "stat-*.csv" -type f 2>/dev/null || true))

if [ ${#stat_files[@]} -eq 0 ]; then
    echo -e "${RED}No stat-*.csv files found in $DATASET_DIR${NC}"
    exit 1
fi

echo -e "${GREEN}Found ${#stat_files[@]} stat files${NC}"

# Show first few files for confirmation
echo "Sample files found:"
for i in $(seq 0 $((${#stat_files[@]} < 5 ? ${#stat_files[@]} - 1 : 4))); do
    echo "  ${stat_files[i]}"
done
if [ ${#stat_files[@]} -gt 5 ]; then
    echo "  ... and $((${#stat_files[@]} - 5)) more files"
fi
echo

echo -e "${YELLOW}Checking header consistency...${NC}"

# Scan all files to find the canonical header (the one with the most columns).
# This handles schema evolution where newer runs added columns (e.g. "setting")
# that older output files don't have: we remap older files instead of aborting.
canonical_header=""
max_cols=0
skipped_empty=0
for stat_file in "${stat_files[@]}"; do
    if [ ! -f "$stat_file" ] || [ ! -s "$stat_file" ]; then
        skipped_empty=$((skipped_empty + 1))
        continue
    fi
    h=$(head -n1 "$stat_file")
    IFS=',' read -ra _cols <<< "$h"
    nc=${#_cols[@]}
    if [[ $nc -gt $max_cols ]]; then
        max_cols=$nc
        canonical_header="$h"
    fi
done

if [[ -z "$canonical_header" ]]; then
    echo -e "${RED}Error: No non-empty stat files found.${NC}"
    exit 1
fi

echo "Canonical header ($max_cols cols): $canonical_header"
[[ $skipped_empty -gt 0 ]] && echo -e "${YELLOW}Warning: $skipped_empty empty/missing files skipped${NC}"

# Report any files whose header differs from canonical (warn, don't abort)
remapped_count=0
for stat_file in "${stat_files[@]}"; do
    [ ! -f "$stat_file" ] || [ ! -s "$stat_file" ] && continue
    current_header=$(head -n1 "$stat_file")
    if [ "$current_header" != "$canonical_header" ]; then
        echo -e "${YELLOW}Warning: header mismatch (will remap): $stat_file${NC}"
        echo "  Canonical: $canonical_header"
        echo "  File:      $current_header"
        remapped_count=$((remapped_count + 1))
    fi
done

if [[ $remapped_count -gt 0 ]]; then
    echo -e "${YELLOW}$remapped_count file(s) have older headers and will be column-remapped (missing cols → empty).${NC}"
else
    echo -e "${GREEN}✓ All headers match canonical${NC}"
fi
echo

echo -e "${YELLOW}Merging stat files into $OUTPUT_FILE...${NC}"

# Helper: emit data rows from a file remapped to the canonical column order.
# Missing columns in the file get an empty string value.
remap_and_append() {
    local file="$1"
    local file_hdr="$2"
    awk -v fhdr="$file_hdr" -v chdr="$canonical_header" '
    BEGIN {
        nf = split(fhdr, fh, ",")
        nc = split(chdr, ch, ",")
        for (i = 1; i <= nc; i++) cidx[ch[i]] = i
        for (i = 1; i <= nf; i++) fmap[i] = cidx[fh[i]]
    }
    NR > 1 {
        nfields = split($0, fields, ",")
        for (i = 1; i <= nc; i++) out[i] = ""
        for (i = 1; i <= nfields; i++)
            if (fmap[i] > 0) out[fmap[i]] = fields[i]
        line = out[1]
        for (i = 2; i <= nc; i++) line = line "," out[i]
        print line
    }
    ' "$file"
}

# Create output file with canonical header
echo "$canonical_header" > "$OUTPUT_FILE"

# Counter for merged rows
total_rows=0

# Append data from all files (skip headers), remapping if needed
for stat_file in "${stat_files[@]}"; do
    if [ ! -f "$stat_file" ] || [ ! -s "$stat_file" ]; then
        continue
    fi

    file_header=$(head -n1 "$stat_file")
    rows_in_file=$(tail -n +2 "$stat_file" | wc -l)

    if [ "$rows_in_file" -gt 0 ]; then
        if [ "$file_header" = "$canonical_header" ]; then
            tail -n +2 "$stat_file" >> "$OUTPUT_FILE"
        else
            remap_and_append "$stat_file" "$file_header" >> "$OUTPUT_FILE"
        fi
        total_rows=$((total_rows + rows_in_file))
        echo "  Added $rows_in_file rows from $(basename "$stat_file")"
    fi
done

echo
echo -e "${GREEN}✓ Successfully merged ${#stat_files[@]} files${NC}"
echo -e "${GREEN}✓ Total data rows: $total_rows${NC}"
echo -e "${GREEN}✓ Output written to: $OUTPUT_FILE${NC}"

# Show summary statistics
echo
echo -e "${BLUE}=== Summary Statistics ===${NC}"
echo "Total stat files processed: ${#stat_files[@]}"
echo "Total data rows merged: $total_rows"
echo "Output file: $OUTPUT_FILE"
echo "Output file size: $(ls -lh "$OUTPUT_FILE" | awk '{print $5}')"

# Show first few rows of output for verification
echo
echo -e "${BLUE}=== First 5 rows of merged data ===${NC}"
head -n 6 "$OUTPUT_FILE" || true

echo
echo -e "${GREEN}Stat collection complete!${NC}"
