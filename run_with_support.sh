#!/bin/bash

# STELAR-MP Run Script with Branch Support Calculation
# Usage: ./run_with_support.sh [input_file] [output_file] [computation_mode] [expansion_method] [support_type] [lambda] [verbose]
#
# Arguments (all optional, will use defaults if not provided):
#   input_file       - Gene trees file (default: all_gt_bs_rooted_100.tre)
#   output_file      - Output species tree file (default: out_support.tre)
#   computation_mode - CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL (default: CPU_PARALLEL)
#   expansion_method - NONE, DISTANCE_ONLY, CONSENSUS_ONLY, DISTANCE_CONSENSUS, FULL (default: DISTANCE_CONSENSUS)
#   support_type     - NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL (default: POSTERIOR)
#   lambda           - Lambda parameter for branch support (default: 0.5)
#   verbose         - true/false for verbose output (default: false)
#
# Examples:
#   ./run_with_support.sh                                           # Use all defaults with posterior support
#   ./run_with_support.sh in.tre out.tre                          # Custom files with posterior support
#   ./run_with_support.sh in.tre out.tre CPU_PARALLEL NONE BOTH   # Both posterior and branch lengths
#   ./run_with_support.sh in.tre out.tre CPU_PARALLEL NONE ALL 1.0 # All annotations with lambda=1.0
#   ./run_with_support.sh in.tre out.tre CPU_PARALLEL NONE DETAILED 0.5 true # Detailed with verbose

# Check for help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "STELAR-MP Run Script with Branch Support Calculation"
    echo ""
    echo "Usage: $0 [input_file] [output_file] [computation_mode] [expansion_method] [support_type] [lambda] [verbose]"
    echo ""
    echo "Arguments (all optional):"
    echo "  input_file       Gene trees file (default: all_gt_bs_rooted_100.tre)"
    echo "  output_file      Output species tree file (default: out_support.tre)"
    echo "  computation_mode CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL (default: CPU_PARALLEL)"
    echo "  expansion_method NONE, DISTANCE_ONLY, CONSENSUS_ONLY, DISTANCE_CONSENSUS, FULL (default: DISTANCE_CONSENSUS)"
    echo "  support_type     NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL (default: POSTERIOR)"
    echo "  lambda           Lambda parameter for branch support (default: 0.5)"
    echo "  verbose         true/false for verbose output (default: false)"
    echo ""
    echo "Examples:"
    echo "  $0                                           # Use all defaults with posterior support"
    echo "  $0 in.tre out.tre                          # Custom files with posterior support"
    echo "  $0 in.tre out.tre CPU_PARALLEL NONE BOTH   # Both posterior and branch lengths"
    echo "  $0 in.tre out.tre CPU_PARALLEL NONE ALL 1.0 # All annotations with lambda=1.0"
    echo "  $0 in.tre out.tre CPU_PARALLEL NONE DETAILED 0.5 true # Detailed with verbose"
    echo ""
    echo "Branch Support Types:"
    echo "  NONE      - No branch support calculation"
    echo "  POSTERIOR - Posterior probabilities only (ASTRAL-style)"
    echo "  DETAILED  - Detailed quartet frequencies and effective N"
    echo "  LENGTH    - Branch lengths only (coalescent units)"
    echo "  BOTH      - Both posterior probabilities and branch lengths"
    echo "  PVALUE    - P-values for polytomy testing"
    echo "  ALL       - All available information"
    echo ""
    echo "Bipartition Expansion Methods:"
    echo "  NONE              - No expansion, original STELAR behavior"
    echo "  DISTANCE_ONLY     - Only distance-based expansion (UPGMA trees)"
    echo "  CONSENSUS_ONLY    - Only consensus-based expansion (greedy consensus)"
    echo "  DISTANCE_CONSENSUS- Both distance and consensus expansion (recommended)"
    echo "  FULL              - All expansion methods including polytomy resolution"
    echo ""
    echo "Lambda Parameter:"
    echo "  Controls the prior in the Beta-Binomial model for branch support"
    echo "  Typical values: 0.5 (default), 1.0, 2.0"
    echo "  Lower values give more conservative support estimates"
    exit 0
fi

# Default configuration
DEFAULT_INPUT_FILE="all_gt_bs_rooted_100.tre"
DEFAULT_OUTPUT_FILE="out_support.tre"
DEFAULT_COMPUTATION_MODE="CPU_PARALLEL"
DEFAULT_EXPANSION_METHOD="DISTANCE_CONSENSUS"
DEFAULT_SUPPORT_TYPE="POSTERIOR"
DEFAULT_LAMBDA="0.5"
DEFAULT_VERBOSE="false"

# Get arguments or use defaults
INPUT_FILE=${1:-$DEFAULT_INPUT_FILE}
OUTPUT_FILE=${2:-$DEFAULT_OUTPUT_FILE}
COMPUTATION_MODE=${3:-$DEFAULT_COMPUTATION_MODE}
EXPANSION_METHOD=${4:-$DEFAULT_EXPANSION_METHOD}
SUPPORT_TYPE=${5:-$DEFAULT_SUPPORT_TYPE}
LAMBDA=${6:-$DEFAULT_LAMBDA}
VERBOSE=${7:-$DEFAULT_VERBOSE}

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

# Validate support type
VALID_SUPPORT_TYPES=("NONE" "POSTERIOR" "DETAILED" "LENGTH" "BOTH" "PVALUE" "ALL")
if [[ ! " ${VALID_SUPPORT_TYPES[@]} " =~ " ${SUPPORT_TYPE} " ]]; then
    echo -e "${RED}Error: Invalid support type '$SUPPORT_TYPE'${NC}"
    echo "Valid support types: ${VALID_SUPPORT_TYPES[*]}"
    exit 1
fi

# Validate lambda parameter
if ! [[ "$LAMBDA" =~ ^[0-9]*\.?[0-9]+$ ]] || (( $(echo "$LAMBDA <= 0" | bc -l) )); then
    echo -e "${RED}Error: Lambda must be a positive number, got '$LAMBDA'${NC}"
    exit 1
fi

echo "=== STELAR-MP with Branch Support ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Computation mode: $COMPUTATION_MODE"
echo "Expansion method: $EXPANSION_METHOD"
echo -e "${BLUE}Support type: $SUPPORT_TYPE${NC}"
echo -e "${BLUE}Lambda parameter: $LAMBDA${NC}"
echo "Verbose output: $VERBOSE"
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
echo "Java version: $(java -version 2>&1 | head -n 1)"
echo

# Build command line arguments
JAVA_ARGS="-i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\" -m \"$COMPUTATION_MODE\""

# Add expansion method if not default
if [ "$EXPANSION_METHOD" != "DISTANCE_CONSENSUS" ]; then
    JAVA_ARGS="$JAVA_ARGS -e \"$EXPANSION_METHOD\""
fi

# Add no-expansion flag if expansion is disabled
if [ "$EXPANSION_METHOD" = "NONE" ]; then
    JAVA_ARGS="$JAVA_ARGS --no-expansion"
fi

# Add branch support parameters
if [ "$SUPPORT_TYPE" != "NONE" ]; then
    JAVA_ARGS="$JAVA_ARGS -s \"$SUPPORT_TYPE\""
    JAVA_ARGS="$JAVA_ARGS --lambda \"$LAMBDA\""
fi

# Add verbose flag if enabled
if [ "$VERBOSE" = "true" ]; then
    JAVA_ARGS="$JAVA_ARGS -v"
fi

# Run the program
echo -e "${YELLOW}Running STELAR-MP with Branch Support...${NC}"
echo -e "${YELLOW}Command: java ... Main $JAVA_ARGS${NC}"
echo

# Record start time
START_TIME=$(date +%s)

# Use different execution based on whether we need GPU support
if [ "$COMPUTATION_MODE" = "GPU_PARALLEL" ]; then
    echo -e "${BLUE}Using GPU mode with CUDA library...${NC}"
    eval "java -Djava.library.path=\"$(pwd)/cuda\" \
         -Djna.debug_load=false \
         -Djna.platform.library.path=\"$(pwd)/cuda\" \
         -cp target/stelar-mp-1.0-SNAPSHOT.jar Main $JAVA_ARGS"
else
    echo -e "${BLUE}Using CPU mode...${NC}"
    eval "java -cp target/stelar-mp-1.0-SNAPSHOT.jar Main $JAVA_ARGS"
fi

# Check exit status
EXIT_CODE=$?
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo
if [ $EXIT_CODE -ne 0 ]; then
    echo -e "${RED}Program execution failed with exit code $EXIT_CODE!${NC}"
    exit $EXIT_CODE
fi

echo -e "${GREEN}Program completed successfully in ${DURATION} seconds!${NC}"
echo "Output written to: $OUTPUT_FILE"

# Display support statistics if support was calculated
if [ "$SUPPORT_TYPE" != "NONE" ]; then
    echo
    echo -e "${BLUE}=== Branch Support Information ===${NC}"
    echo "Support type: $SUPPORT_TYPE"
    echo "Lambda parameter: $LAMBDA"
    
    # Show a preview of the output tree
    echo
    echo -e "${BLUE}Output tree preview (first 200 characters):${NC}"
    head -c 200 "$OUTPUT_FILE" && echo "..."
    echo
    
    # Provide interpretation guide
    echo -e "${BLUE}=== Support Value Interpretation ===${NC}"
    case $SUPPORT_TYPE in
        "POSTERIOR")
            echo "Values represent posterior probabilities (0.0 to 1.0)"
            echo "  >= 0.95: Strong support"
            echo "  >= 0.75: Moderate support"
            echo "  <  0.75: Weak support"
            ;;
        "LENGTH")
            echo "Values represent branch lengths in coalescent units"
            echo "Longer branches indicate more evolutionary time"
            ;;
        "BOTH")
            echo "Format: posterior:length"
            echo "Posterior probabilities followed by branch lengths"
            ;;
        "DETAILED"|"ALL")
            echo "Format includes quartet frequencies and effective sample sizes"
            echo "See ASTRAL documentation for detailed interpretation"
            ;;
        "PVALUE")
            echo "P-values for polytomy test (lower = more significant)"
            echo "  < 0.05: Reject polytomy (binary split supported)"
            echo "  >= 0.05: Cannot reject polytomy"
            ;;
    esac
fi

echo
echo -e "${GREEN}=== Run Complete ===${NC}"
