#!/bin/bash

# Memory-Optimized STELAR Run Script
# Usage: ./run_memory_optimized.sh [input_file] [output_file] [computation_mode] [expansion_method]
#
# Arguments (compatible with original run.sh):
#   input_file       - Gene trees file (default: all_gt_bs_rooted_100.tre)
#   output_file      - Output species tree file (default: out.tre) 
#   computation_mode - CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL (default: CPU_PARALLEL)
#   expansion_method - Currently ignored (memory-optimized version uses built-in expansion)
#
# Examples:
#   ./run_memory_optimized.sh 1kp.tre out.tre GPU_PARALLEL NONE
#   ./run_memory_optimized.sh all_gt_bs_rooted_37.tre result.tre CPU_PARALLEL
#   ./run_memory_optimized.sh input.tre                    # Uses defaults for other params

# Check for help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Memory-Optimized STELAR Run Script"
    echo ""
    echo "Usage: $0 [input_file] [output_file] [computation_mode] [expansion_method]"
    echo ""
    echo "Arguments:"
    echo "  input_file       - Gene trees file (default: all_gt_bs_rooted_100.tre)"
    echo "  output_file      - Output species tree file (default: out.tre)"
    echo "  computation_mode - CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL (default: CPU_PARALLEL)"
    echo "  expansion_method - Currently ignored (built-in expansion used)"
    echo ""
    echo "Examples:"
    echo "  $0 1kp.tre out.tre GPU_PARALLEL NONE"
    echo "  $0 all_gt_bs_rooted_37.tre result.tre CPU_PARALLEL"
    echo "  $0 input.tre"
    echo ""
    exit 0
fi

# Set default values (compatible with original run.sh)
INPUT_FILE=${1:-"all_gt_bs_rooted_100.tre"}
OUTPUT_FILE=${2:-"out.tre"}
COMPUTATION_MODE=${3:-"CPU_PARALLEL"}
EXPANSION_METHOD=${4:-"NONE"}

# Validate computation mode
case "$COMPUTATION_MODE" in
    CPU_SINGLE|CPU_PARALLEL|GPU_PARALLEL)
        ;;
    *)
        echo "Error: Invalid computation mode '$COMPUTATION_MODE'"
        echo "Valid modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL"
        exit 1
        ;;
esac

echo "Running Memory-Optimized STELAR..."
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Computation mode: $COMPUTATION_MODE"
echo "Expansion method: $EXPANSION_METHOD (currently ignored)"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

# Set library path for CUDA libraries
export LD_LIBRARY_PATH=./cuda:$LD_LIBRARY_PATH

# Get the path to the Maven repository
M2_REPO="$HOME/.m2/repository"

# Run the memory-optimized implementation with proper classpath
java -cp ".:bin:target/stelar-mp-1.0-SNAPSHOT.jar:$M2_REPO/net/java/dev/jna/jna/5.13.0/jna-5.13.0.jar:$M2_REPO/net/java/dev/jna/jna-platform/5.13.0/jna-platform-5.13.0.jar" \
     -Xmx8g MemoryOptimizedMain "$INPUT_FILE" "$OUTPUT_FILE" "$COMPUTATION_MODE"

if [ $? -eq 0 ]; then
    echo "Success! Output written to: $OUTPUT_FILE"
else
    echo "Error: Memory-optimized STELAR execution failed!"
    exit 1
fi
