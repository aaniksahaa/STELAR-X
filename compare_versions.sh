#!/bin/bash

# Configuration
INPUT_FILE="all_gt_bs_rooted_37.tre"  # Change this to your input file
PARALLEL_OUTPUT="parallel_out.tre"
ORIGINAL_OUTPUT="original_out.tre"
ORIGINAL_JAR="STELAR.jar"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo "=== STELAR Version Comparison Script ==="
echo "Input file: $INPUT_FILE"
echo

# Function to format time
format_time() {
    local seconds=$1
    local minutes=$((seconds / 60))
    local remaining_seconds=$((seconds % 60))
    printf "%02d:%02d" $minutes $remaining_seconds
}

# Run parallel version
echo "Running parallel version..."
start_time=$(date +%s)
./run.sh "$INPUT_FILE" "$PARALLEL_OUTPUT"
parallel_time=$(( $(date +%s) - start_time ))
echo -e "${GREEN}Parallel version completed in $(format_time $parallel_time)${NC}"
echo

# Run original version
echo "Running original version..."
start_time=$(date +%s)
java -jar "$ORIGINAL_JAR" -i "$INPUT_FILE" -o "$ORIGINAL_OUTPUT"
original_time=$(( $(date +%s) - start_time ))
echo -e "${GREEN}Original version completed in $(format_time $original_time)${NC}"
echo

# Calculate speedup
speedup=$(echo "scale=2; $original_time / $parallel_time" | bc)
echo "=== Performance Comparison ==="
echo "Original version time: $(format_time $original_time)"
echo "Parallel version time: $(format_time $parallel_time)"
echo -e "Speedup: ${GREEN}${speedup}x${NC}"
echo

# Compare output trees
echo "=== Tree Comparison ==="
if cmp -s "$PARALLEL_OUTPUT" "$ORIGINAL_OUTPUT"; then
    echo -e "${GREEN}Output trees are identical!${NC}"
else
    echo -e "${RED}Output trees are different!${NC}"
    echo "Differences:"
    diff "$PARALLEL_OUTPUT" "$ORIGINAL_OUTPUT"
fi

# Cleanup
echo
echo "=== Output Files ==="
echo "Parallel output: $PARALLEL_OUTPUT"
echo "Original output: $ORIGINAL_OUTPUT" 