#!/bin/bash

# Default configuration
DEFAULT_INPUT_FILE="all_gt_bs_rooted_48.tre"
DEFAULT_OUTPUT_FILE="out.tre"
DEFAULT_COMPUTATION_MODE="CPU_PARALLEL"

# Get arguments or use defaults
INPUT_FILE=${1:-$DEFAULT_INPUT_FILE}
OUTPUT_FILE=${2:-$DEFAULT_OUTPUT_FILE}
COMPUTATION_MODE=${3:-$DEFAULT_COMPUTATION_MODE}

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Validate computation mode
VALID_MODES=("CPU_SINGLE" "CPU_PARALLEL")  # GPU_PARALLEL removed for now
if [[ ! " ${VALID_MODES[@]} " =~ " ${COMPUTATION_MODE} " ]]; then
    echo -e "${RED}Error: Invalid computation mode '$COMPUTATION_MODE'${NC}"
    echo "Valid modes: ${VALID_MODES[*]}"
    exit 1
fi

echo "=== STELAR-MP Run Script ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Computation mode: $COMPUTATION_MODE"
echo

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}"
    exit 1
fi

# Check if binaries exist
if [ ! -f "target/stelar-mp-1.0-SNAPSHOT.jar" ] || [ ! -d "bin" ]; then
    echo -e "${RED}Error: Binaries not found. Please run build.sh first.${NC}"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Run the program
echo -e "${YELLOW}Running STELAR-MP...${NC}"
java -cp target/stelar-mp-1.0-SNAPSHOT.jar Main -i "$INPUT_FILE" -o "$OUTPUT_FILE" -m "$COMPUTATION_MODE"

if [ $? -ne 0 ]; then
    echo -e "${RED}Program execution failed!${NC}"
    exit 1
fi

echo -e "${GREEN}Program completed successfully!${NC}"
echo "Output written to: $OUTPUT_FILE" 