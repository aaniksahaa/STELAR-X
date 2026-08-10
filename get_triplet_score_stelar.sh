#!/bin/bash

# STELAR-X Triplet Score Calculator
# Usage: ./get_triplet_score_stelar.sh -i <gene_trees_file> -st <species_tree_file> [computation_mode_flags]

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default configuration
GENE_TREES_FILE=""
SPECIES_TREE_FILE=""
COMPUTATION_MODE="GPU_PARALLEL"
STELAR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION="$("$STELAR_ROOT/project-version.sh")"
JAR_PATH="$STELAR_ROOT/target/stelar-x-${VERSION}.jar"

# Function to display usage
usage() {
    echo "STELAR-X Triplet Score Calculator"
    echo ""
    echo "Calculates the triplet score between gene trees and a given species tree."
    echo ""
    echo "Usage: $0 -i <gene_trees_file> -st <species_tree_file> [--cpu | --cpu-single | --gpu]"
    echo ""
    echo "Arguments:"
    echo "  -i, --input        File containing gene trees in Newick format"
    echo "  -st, --species-tree File containing the species tree to score"
    echo "  --cpu              Run in CPU parallel mode"
    echo "  --cpu-single       Run in CPU single-threaded mode"
    echo "  --gpu              Run in GPU parallel mode (default)"
    echo "  -h, --help         Show this help message"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            GENE_TREES_FILE="$2"
            shift 2
            ;;
        -st|--species-tree)
            SPECIES_TREE_FILE="$2"
            shift 2
            ;;
        --cpu)
            COMPUTATION_MODE="CPU_PARALLEL"
            shift
            ;;
        --cpu-single)
            COMPUTATION_MODE="CPU_SINGLE"
            shift
            ;;
        --gpu|--gpu-parallel)
            COMPUTATION_MODE="GPU_PARALLEL"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo -e "${RED}Error: Unknown option: $1${NC}"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$GENE_TREES_FILE" || -z "$SPECIES_TREE_FILE" ]]; then
    echo -e "${RED}Error: Both gene trees file and species tree file are required.${NC}"
    usage
fi

# Check if files exist
if [[ ! -f "$GENE_TREES_FILE" ]]; then
    echo -e "${RED}Error: Gene trees file '$GENE_TREES_FILE' does not exist.${NC}"
    exit 1
fi

if [[ ! -f "$SPECIES_TREE_FILE" ]]; then
    echo -e "${RED}Error: Species tree file '$SPECIES_TREE_FILE' does not exist.${NC}"
    exit 1
fi

# Check if binaries exist
if [[ ! -f "$JAR_PATH" ]]; then
    echo -e "${RED}Error: JAR file not found. Please run build.sh first.${NC}"
    exit 1
fi

echo "=== STELAR-X Triplet Score Calculator ==="
echo "Gene trees file:   $GENE_TREES_FILE"
echo "Species tree file: $SPECIES_TREE_FILE"
echo "Computation mode:  $COMPUTATION_MODE"
echo

# Run the score calculation
echo -e "${YELLOW}Calculating triplet score...${NC}"
echo

java -Xms4g -Xmx128g \
    -Djava.library.path="$STELAR_ROOT/cuda" \
    -Djna.debug_load=false \
    -Djna.platform.library.path="$STELAR_ROOT/cuda" \
    -cp "$JAR_PATH" \
    Main -i "$GENE_TREES_FILE" -c "$SPECIES_TREE_FILE" -m "$COMPUTATION_MODE"

EXIT_CODE=$?

if [[ $EXIT_CODE -ne 0 ]]; then
    echo -e "${RED}Score calculation failed!${NC}"
    exit $EXIT_CODE
fi

echo -e "${GREEN}Score calculation completed successfully!${NC}"
