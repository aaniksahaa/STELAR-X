#!/bin/bash

# Configuration
INPUT_FILE="all_gt_bs_rooted_200.tre"  # Change this to your input file
OUTPUT_DIR="benchmark_results"
# ORIGINAL_JAR="/mnt/H/Research/STELAR-extension/STELAR/STELAR.jar"
ORIGINAL_JAR="/mnt/nvme/anik/phylogeny/STELAR-MP/original/STELAR/STELAR.jar"

# Computation modes to test
COMPUTATION_MODES=("CPU_SINGLE" "CPU_PARALLEL" "GPU_PARALLEL" "ORIGINAL")

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo
echo "=== STELAR-MP Benchmark Script ==="
echo "Input file: $INPUT_FILE"
echo

# Debug information for GPU
echo -e "${YELLOW}Debug Information:${NC}"
echo "Current directory: $(pwd)"
echo "Library path: $(pwd)/cuda"
echo "Library exists: $(if [ -f "$(pwd)/cuda/libweight_calc.so" ]; then echo "Yes"; else echo "No"; fi)"
echo "Library permissions: $(ls -l "$(pwd)/cuda/libweight_calc.so" 2>/dev/null || echo "Not found")"
echo

# Get current time in milliseconds
get_time_ms() {
    echo $(($(date +%s%3N)))
}

# Format milliseconds to MM:SS.mmm
format_time() {
    local total_ms=$1
    local minutes=$((total_ms / 60000))
    local seconds=$(((total_ms % 60000) / 1000))
    local milliseconds=$((total_ms % 1000))
    printf "%02d:%02d.%03d" "$minutes" "$seconds" "$milliseconds"
}

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file not found at $INPUT_FILE${NC}"
    exit 1
fi

# Check if original STELAR jar exists
if [ ! -f "$ORIGINAL_JAR" ]; then
    echo -e "${RED}Error: Original STELAR jar not found at $ORIGINAL_JAR${NC}"
    exit 1
fi

# Check if binaries exist
if [ ! -f "target/stelar-mp-1.0-SNAPSHOT.jar" ] || [ ! -d "bin" ]; then
    echo -e "${RED}Error: Binaries not found. Please run build.sh first.${NC}"
    exit 1
fi

# Run benchmark for each mode
declare -A results
for mode in "${COMPUTATION_MODES[@]}"; do
    echo -e "\n${YELLOW}Running $mode version...${NC}"
    output_file="$OUTPUT_DIR/${mode,,}_out.tre"
    
    start_time=$(get_time_ms)
    if [ "$mode" = "ORIGINAL" ]; then
        java -jar "$ORIGINAL_JAR" -i "$INPUT_FILE" -o "$output_file"
    else
        if [ "$mode" = "GPU_PARALLEL" ]; then
            java -Djava.library.path="$(pwd)/cuda" \
                 -Djna.debug_load=true \
                 -Djna.debug_load.jna=true \
                 -Djna.platform.library.path="$(pwd)/cuda" \
                 -Djna.memory.contiguous=true \
                 -Djna.memory.contiguous.alignment=8 \
                 -Djna.memory.contiguous.size=1024 \
                 -cp target/stelar-mp-1.0-SNAPSHOT.jar Main -i "$INPUT_FILE" -o "$output_file" -m "$mode"
        else
            java -cp target/stelar-mp-1.0-SNAPSHOT.jar Main -i "$INPUT_FILE" -o "$output_file" -m "$mode"
        fi
    fi
    end_time=$(get_time_ms)
    execution_time=$((end_time - start_time))
    
    if [ $? -eq 0 ]; then
        results[$mode]=$execution_time
        echo -e "${GREEN}$mode version completed in $(format_time $execution_time)${NC}"
    else
        echo -e "${RED}$mode version failed!${NC}"
    fi
done

# Print performance comparison
echo -e "\n=== Performance Comparison ==="
baseline=${results["ORIGINAL"]}
if [ ! -z "$baseline" ]; then
    echo "ORIGINAL (baseline): $(format_time $baseline)"
    for mode in "${COMPUTATION_MODES[@]}"; do
        if [ "$mode" != "ORIGINAL" ] && [ ! -z "${results[$mode]}" ]; then
            speedup=$(awk "BEGIN { printf \"%.2f\", $baseline / ${results[$mode]} }")
            echo "$mode: $(format_time ${results[$mode]}) (Speedup: ${speedup}x)"
        fi
    done
fi

# Compare output trees
echo -e "\n=== Tree Comparison ==="
for mode in "${COMPUTATION_MODES[@]}"; do
    if [ "$mode" != "ORIGINAL" ]; then
        output_file="$OUTPUT_DIR/${mode,,}_out.tre"
        original_file="$OUTPUT_DIR/original_out.tre"
        if [ -f "$output_file" ] && [ -f "$original_file" ]; then
            echo -e "\nComparing $mode with ORIGINAL..."
            java -cp bin CompareTrees "$output_file" "$original_file"
            if [ $? -eq 0 ]; then
                echo -e "${GREEN}Trees are identical!${NC}"
            else
                echo -e "${RED}Trees are different!${NC}"
            fi
        fi
    fi
done

# Output files
echo -e "\n=== Output Files ==="
for mode in "${COMPUTATION_MODES[@]}"; do
    output_file="$OUTPUT_DIR/${mode,,}_out.tre"
    if [ -f "$output_file" ]; then
        echo "$mode output: $output_file"
    fi
done
