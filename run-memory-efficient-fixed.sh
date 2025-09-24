#!/bin/bash

# Memory-Efficient STELAR-MP Run Script
# Uses range-based processing during gene tree preprocessing to reduce memory usage
# while keeping the rest of the pipeline (weight calculation, DP) unchanged.

# Default configuration
DEFAULT_INPUT_FILE="all_gt_bs_rooted_37.tre"
DEFAULT_OUTPUT_FILE="memory_efficient_out.tre"
DEFAULT_COMPUTATION_MODE="CPU_PARALLEL"
DEFAULT_EXPANSION_METHOD="NONE"

# Get arguments or use defaults
INPUT_FILE=${1:-$DEFAULT_INPUT_FILE}
OUTPUT_FILE=${2:-$DEFAULT_OUTPUT_FILE}
COMPUTATION_MODE=${3:-$DEFAULT_COMPUTATION_MODE}
EXPANSION_METHOD=${4:-$DEFAULT_EXPANSION_METHOD}

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

echo "=== Memory-Efficient STELAR-MP Run Script ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Computation mode: $COMPUTATION_MODE"
echo "Expansion method: $EXPANSION_METHOD"
echo "Memory optimization: Gene tree processing only"
echo

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}"
    exit 1
fi

# Build the project first
echo -e "${YELLOW}Building project...${NC}"
if ! ./build.sh > build.log 2>&1; then
    echo -e "${RED}Build failed! Check build.log for details.${NC}"
    tail -20 build.log
    exit 1
fi
echo -e "${GREEN}Build completed successfully!${NC}"

# Check if JAR exists
if [ ! -f "target/stelar-mp-1.0-SNAPSHOT.jar" ]; then
    echo -e "${RED}Error: JAR file not found. Build may have failed.${NC}"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Debug information
echo -e "${YELLOW}Debug Information:${NC}"
echo "Current directory: $(pwd)"
echo "JAR file exists: $(if [ -f "target/stelar-mp-1.0-SNAPSHOT.jar" ]; then echo "Yes"; else echo "No"; fi)"
echo "Input file size: $(du -h "$INPUT_FILE" 2>/dev/null | cut -f1 || echo "Unknown")"
echo

# Run the memory-efficient program
echo -e "${YELLOW}Running Memory-Efficient STELAR-MP...${NC}"
echo -e "${YELLOW}Command: java ... MemoryEfficientMain ...${NC}"
echo

# Record start time
START_TIME=$(date +%s)

# Run with optimized JVM settings
java -Xms2g -Xmx8g \
     -XX:+UseG1GC \
     -XX:MaxGCPauseMillis=200 \
     -Djava.library.path="$(pwd)/cuda" \
     -cp target/stelar-mp-1.0-SNAPSHOT.jar MemoryEfficientMain \
     -i "$INPUT_FILE" \
     -o "$OUTPUT_FILE" \
     -m "$COMPUTATION_MODE" \
     -e "$EXPANSION_METHOD"

# Record end time and calculate duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

if [ $? -ne 0 ]; then
    echo -e "${RED}Program execution failed!${NC}"
    echo
    echo -e "${YELLOW}Checking for common issues:${NC}"
    
    # Check if output file was partially created
    if [ -f "$OUTPUT_FILE" ]; then
        echo "Partial output file created: $(du -h "$OUTPUT_FILE" | cut -f1)"
    fi
    
    exit 1
fi

echo
echo -e "${GREEN}Memory-Efficient Program completed successfully!${NC}"
echo "Execution time: ${DURATION} seconds"
echo "Output written to: $OUTPUT_FILE"

# Display output file information
if [ -f "$OUTPUT_FILE" ]; then
    OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
    TREE_COUNT=$(grep -c ';' "$OUTPUT_FILE" 2>/dev/null || echo "1")
    echo "Output file size: $OUTPUT_SIZE"
    echo "Number of trees: $TREE_COUNT"
    
    # Show first few characters of the output tree
    echo "Tree preview:"
    head -c 200 "$OUTPUT_FILE" && echo "..."
fi

echo
echo -e "${YELLOW}Memory Optimization Benefits:${NC}"
echo "• Range-based processing during gene tree preprocessing"
echo "• Identifies unique bipartitions before converting to BitSets"
echo "• Reduces memory usage during the most memory-intensive phase"
echo "• Keeps original inference pipeline unchanged for reliability"
echo "• Maintains same accuracy as original STELAR-MP"

# Cleanup
if [ -f "build.log" ]; then
    rm build.log
fi

echo
echo -e "${GREEN}Memory-efficient analysis complete!${NC}"
echo
echo "To compare with original STELAR-MP:"
echo "./run.sh \"$INPUT_FILE\" \"original_$OUTPUT_FILE\" \"$COMPUTATION_MODE\" \"$EXPANSION_METHOD\""
echo
echo "To run with different parameters:"
echo "$0 [input_file] [output_file] [CPU_SINGLE|CPU_PARALLEL|GPU_PARALLEL] [NONE|DISTANCE_ONLY|etc]"
