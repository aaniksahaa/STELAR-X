#!/bin/bash

# Memory-Optimized STELAR-MP Run Script
# Usage: ./run-memory-optimized.sh [input_file] [output_file] [computation_mode]
#
# Arguments (all optional):
#   input_file       - Gene trees file (default: all_gt_bs_rooted_37.tre)
#   output_file      - Output species tree file (default: memory_opt_out.tre)
#   computation_mode - CPU_SINGLE, CPU_PARALLEL (default: CPU_PARALLEL)

# Default configuration
DEFAULT_INPUT_FILE="all_gt_bs_rooted_37.tre"
DEFAULT_OUTPUT_FILE="memory_opt_out.tre"
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
VALID_MODES=("CPU_SINGLE" "CPU_PARALLEL")
if [[ ! " ${VALID_MODES[@]} " =~ " ${COMPUTATION_MODE} " ]]; then
    echo -e "${RED}Error: Invalid computation mode '$COMPUTATION_MODE'${NC}"
    echo "Valid modes for memory-optimized version: ${VALID_MODES[*]}"
    exit 1
fi

echo "=== Memory-Optimized STELAR-MP Run Script ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Computation mode: $COMPUTATION_MODE"
echo "Memory optimization: ENABLED"
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

# Run the memory-optimized program
echo -e "${YELLOW}Running Memory-Optimized STELAR-MP...${NC}"
echo -e "${YELLOW}Command: java ... MemoryOptimizedMain ...${NC}"
echo

# Record start time
START_TIME=$(date +%s)

# Run with memory profiling enabled
java -Xms2g -Xmx8g \
     -XX:+UseG1GC \
     -XX:MaxGCPauseMillis=200 \
     -cp target/stelar-mp-1.0-SNAPSHOT.jar MemoryOptimizedMain \
     -i "$INPUT_FILE" \
     -o "$OUTPUT_FILE" \
     -m "$COMPUTATION_MODE" \
     -profile

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
    
    # Check GC log for memory issues
    if [ -f "gc.log" ]; then
        echo "GC log created - checking for OutOfMemoryError..."
        if grep -q "OutOfMemoryError" gc.log; then
            echo -e "${RED}OutOfMemoryError detected! Try increasing heap size with -Xmx16g${NC}"
        fi
    fi
    
    exit 1
fi

echo
echo -e "${GREEN}Memory-Optimized Program completed successfully!${NC}"
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

# Memory usage summary
echo
echo -e "${YELLOW}Memory Usage Summary:${NC}"
if [ -f "gc.log" ]; then
    # Extract some basic GC statistics
    TOTAL_GC_TIME=$(grep "Total time for which application threads were stopped" gc.log | tail -1 | awk '{print $10}' || echo "N/A")
    MAX_HEAP=$(grep -o "PSYoungGen.*->" gc.log | tail -1 || echo "N/A")
    echo "GC log available - Total GC pause time: $TOTAL_GC_TIME seconds"
    echo "GC log location: gc.log"
else
    echo "No GC log generated"
fi

# Compare with traditional approach (if available)
echo
echo -e "${YELLOW}Memory Optimization Benefits:${NC}"
echo "• Reduced memory complexity from O(n²k) to O(nk)"
echo "• Range-based bipartition representation instead of BitSets"
echo "• Hash-based equality checking for efficient uniqueness detection"
echo "• Improved scalability for large datasets"

# Cleanup
if [ -f "build.log" ]; then
    rm build.log
fi

echo
echo -e "${GREEN}Memory-optimized analysis complete!${NC}"
echo
echo "To compare with original STELAR-MP:"
echo "./run.sh \"$INPUT_FILE\" \"original_$OUTPUT_FILE\" \"$COMPUTATION_MODE\" NONE"
echo
echo "To run with different parameters:"
echo "$0 [input_file] [output_file] [CPU_SINGLE|CPU_PARALLEL]"
