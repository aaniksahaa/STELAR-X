#!/bin/bash

# Configuration
INPUT_FILE=$1
OUTPUT_FILE=$2
COMPUTATION_MODE=${3:-GPU_PARALLEL}  # Default to GPU mode if not specified

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Validate computation mode
VALID_MODES=("CPU_SINGLE" "CPU_PARALLEL" "GPU_PARALLEL")
if [[ ! " ${VALID_MODES[@]} " =~ " ${COMPUTATION_MODE} " ]]; then
    echo -e "${RED}Error: Invalid computation mode '$COMPUTATION_MODE'${NC}"
    echo "Valid modes: ${VALID_MODES[*]}"
    exit 1
fi

echo "=== STELAR-MP Build and Run Script ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Computation mode: $COMPUTATION_MODE"
echo

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Compile CUDA code only if in GPU mode
if [ "$COMPUTATION_MODE" = "GPU_PARALLEL" ]; then
    echo -e "${YELLOW}Compiling CUDA code...${NC}"
    cd cuda
    make clean
    make
    if [ $? -ne 0 ]; then
        echo -e "${RED}CUDA compilation failed!${NC}"
        exit 1
    fi
    cd ..
    echo -e "${GREEN}CUDA compilation successful${NC}"
    echo
fi

# Compile Java code
echo -e "${YELLOW}Compiling Java code...${NC}"
mvn clean package
if [ $? -ne 0 ]; then
    echo -e "${RED}Java compilation failed!${NC}"
    exit 1
fi
echo -e "${GREEN}Java compilation successful${NC}"
echo

# Run the program
echo -e "${YELLOW}Running STELAR-MP...${NC}"
if [ "$COMPUTATION_MODE" = "GPU_PARALLEL" ]; then
    java -Djava.library.path=cuda -cp target/stelar-mp-1.0-SNAPSHOT.jar Main -i "$INPUT_FILE" -o "$OUTPUT_FILE" -m "$COMPUTATION_MODE"
else
    java -cp target/stelar-mp-1.0-SNAPSHOT.jar Main -i "$INPUT_FILE" -o "$OUTPUT_FILE" -m "$COMPUTATION_MODE"
fi

if [ $? -ne 0 ]; then
    echo -e "${RED}Program execution failed!${NC}"
    exit 1
fi

echo -e "${GREEN}Program completed successfully!${NC}"
echo "Output written to: $OUTPUT_FILE" 