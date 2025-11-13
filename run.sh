#!/bin/bash


# STELAR-MP Run Script with Bipartition Expansion Support
# Usage: ./run.sh [input_file] [output_file] [computation_mode] [expansion_method] [distance_method] [verbose]
#
# Arguments (all optional, will use defaults if not provided):
#   input_file       - Gene trees file (default: all_gt_bs_rooted_100.tre)
#   output_file      - Output species tree file (default: out.tre)
#   computation_mode - CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL (default: CPU_PARALLEL)
#   expansion_method - NONE, DISTANCE_ONLY, CONSENSUS_ONLY, DISTANCE_CONSENSUS, FULL (default: DISTANCE_CONSENSUS)
#   distance_method  - UPGMA, NEIGHBOR_JOINING, BOTH (default: UPGMA)
#   verbose         - true/false for verbose expansion output (default: false)
#
# Examples:
#   ./run.sh                                                    # Use all defaults
#   ./run.sh in.tre out.tre                              # Custom input/output files
#   ./run.sh in.tre out.tre CPU_PARALLEL DISTANCE_ONLY   # Distance expansion only
#   ./run.sh in.tre out.tre GPU_PARALLEL NONE            # No expansion (original STELAR)
#   ./run.sh in.tre out.tre GPU_PARALLEL FULL UPGMA true # Full expansion with verbose output


#  ./run.sh all_gt_bs_rooted_11.tre out.tre GPU_PARALLEL NONE
#  ./run.sh all_gt_bs_rooted_37.tre out.tre GPU_PARALLEL NONE
#  ./run.sh all_gt_bs_rooted_200.tre out.tre GPU_PARALLEL NONE
#  ./run.sh all_gt_bs_rooted_200.tre out.tre GPU_PARALLEL FULL UPGMA true

# Check for help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "STELAR-MP Run Script with Bipartition Expansion Support"
    echo ""
    echo "Usage: $0 [input_file] [output_file] [computation_mode] [expansion_method] [distance_method] [verbose]"
    echo ""
    echo "Arguments (all optional):"
    echo "  input_file       Gene trees file (default: all_gt_bs_rooted_100.tre)"
    echo "  output_file      Output species tree file (default: out.tre)"
    echo "  computation_mode CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL (default: CPU_PARALLEL)"
    echo "  expansion_method NONE, DISTANCE_ONLY, CONSENSUS_ONLY, DISTANCE_CONSENSUS, FULL (default: DISTANCE_CONSENSUS)"
    echo "  distance_method  UPGMA, NEIGHBOR_JOINING, BOTH (default: UPGMA)"
    echo "  verbose         true/false for verbose expansion output (default: false)"
    echo ""
    echo "Examples:"
    echo "  $0                                                    # Use all defaults"
    echo "  $0 in.tre out.tre                              # Custom input/output files"
    echo "  $0 in.tre out.tre CPU_PARALLEL DISTANCE_ONLY   # Distance expansion only"
    echo "  $0 in.tre out.tre CPU_PARALLEL NONE            # No expansion (original STELAR)"
    echo "  $0 in.tre out.tre GPU_PARALLEL FULL UPGMA true # Full expansion with verbose output"
    echo ""
    echo "Bipartition Expansion Methods:"
    echo "  NONE              - No expansion, original STELAR behavior"
    echo "  DISTANCE_ONLY     - Only distance-based expansion (UPGMA trees)"
    echo "  CONSENSUS_ONLY    - Only consensus-based expansion (greedy consensus)"
    echo "  DISTANCE_CONSENSUS- Both distance and consensus expansion (recommended)"
    echo "  FULL              - All expansion methods including polytomy resolution"
    exit 0
fi

# Default configuration
DEFAULT_INPUT_FILE="all_gt_bs_rooted_100.tre"
DEFAULT_OUTPUT_FILE="out.tre"
DEFAULT_COMPUTATION_MODE="CPU_PARALLEL"
DEFAULT_EXPANSION_METHOD="DISTANCE_CONSENSUS"
DEFAULT_DISTANCE_METHOD="UPGMA"
DEFAULT_VERBOSE_EXPANSION="false"

# Get arguments or use defaults
INPUT_FILE=${1:-$DEFAULT_INPUT_FILE}
OUTPUT_FILE=${2:-$DEFAULT_OUTPUT_FILE}
COMPUTATION_MODE=${3:-$DEFAULT_COMPUTATION_MODE}
EXPANSION_METHOD=${4:-$DEFAULT_EXPANSION_METHOD}
DISTANCE_METHOD=${5:-$DEFAULT_DISTANCE_METHOD}
VERBOSE_EXPANSION=${6:-$DEFAULT_VERBOSE_EXPANSION}

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

# Validate expansion method
VALID_EXPANSION_METHODS=("NONE" "DISTANCE_ONLY" "CONSENSUS_ONLY" "DISTANCE_CONSENSUS" "FULL")
if [[ ! " ${VALID_EXPANSION_METHODS[@]} " =~ " ${EXPANSION_METHOD} " ]]; then
    echo -e "${RED}Error: Invalid expansion method '$EXPANSION_METHOD'${NC}"
    echo "Valid expansion methods: ${VALID_EXPANSION_METHODS[*]}"
    exit 1
fi

# Validate distance method
VALID_DISTANCE_METHODS=("UPGMA" "NEIGHBOR_JOINING" "BOTH")
if [[ ! " ${VALID_DISTANCE_METHODS[@]} " =~ " ${DISTANCE_METHOD} " ]]; then
    echo -e "${RED}Error: Invalid distance method '$DISTANCE_METHOD'${NC}"
    echo "Valid distance methods: ${VALID_DISTANCE_METHODS[*]}"
    exit 1
fi

echo "=== STELAR-MP Run Script ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Computation mode: $COMPUTATION_MODE"
echo "Expansion method: $EXPANSION_METHOD"
echo "Distance method: $DISTANCE_METHOD"
echo "Verbose expansion: $VERBOSE_EXPANSION"
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

# Debug information
echo -e "${YELLOW}Debug Information:${NC}"
echo "Current directory: $(pwd)"
echo "Library path: $(pwd)/cuda"
echo "Library exists: $(if [ -f "$(pwd)/cuda/libweight_calc.so" ]; then echo "Yes"; else echo "No"; fi)"
echo "Library permissions: $(ls -l "$(pwd)/cuda/libweight_calc.so" 2>/dev/null || echo "Not found")"
echo

# Build command line arguments
JAVA_ARGS="-i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\" -m \"$COMPUTATION_MODE\""

# Add expansion method
if [ "$EXPANSION_METHOD" != "DISTANCE_CONSENSUS" ]; then
    JAVA_ARGS="$JAVA_ARGS -e \"$EXPANSION_METHOD\""
fi

# Add distance method if not default
if [ "$DISTANCE_METHOD" != "UPGMA" ]; then
    JAVA_ARGS="$JAVA_ARGS -d \"$DISTANCE_METHOD\""
fi

# Add verbose flag if enabled
if [ "$VERBOSE_EXPANSION" = "true" ]; then
    JAVA_ARGS="$JAVA_ARGS -v"
fi

# Add no-expansion flag if expansion is disabled
if [ "$EXPANSION_METHOD" = "NONE" ]; then
    JAVA_ARGS="$JAVA_ARGS --no-expansion"
fi

# Run the program with the library path set for this run only
echo -e "${YELLOW}Running STELAR-MP...${NC}"
echo -e "${YELLOW}Command: java ... Main $JAVA_ARGS${NC}"
echo

eval "java -Xms4g -Xmx44g -Djava.library.path=\"$(pwd)/cuda\" \
     -Djna.debug_load=true \
     -Djna.debug_load.jna=true \
     -Djna.platform.library.path=\"$(pwd)/cuda\" \
     -Djna.memory.contiguous=true \
     -Djna.memory.contiguous.alignment=8 \
     -Djna.memory.contiguous.size=1024 \
     -cp target/stelar-mp-1.0-SNAPSHOT.jar Main $JAVA_ARGS"

if [ $? -ne 0 ]; then
    echo -e "${RED}Program execution failed!${NC}"
    exit 1
fi

echo -e "${GREEN}Program completed successfully!${NC}"
echo "Output written to: $OUTPUT_FILE" 