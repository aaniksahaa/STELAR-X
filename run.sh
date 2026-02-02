#!/bin/bash


# STELAR-X Run Script with Mixed Bipartition Support
# Usage: ./run.sh --input <file> --output <file> [options]
#
# Options:
#   --input, -i         Gene trees file (required)
#   --output, -o        Output species tree file (required)
#   --cpu               Computation mode: CPU_SINGLE
#   --cpu-parallel       Computation mode: CPU_PARALLEL
#   --gpu               Computation mode: GPU_PARALLEL
#   --mode, -m          Computation mode: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL (default: GPU_PARALLEL)
#   --expansion, -e     Enable mixed bipartitions (default: OFF)
#   --verbose, -v       Verbose expansion output (default: off)
#
# Examples:
#   ./run.sh --input in.tre --output out.tre
#   ./run.sh --input in.tre --output out.tre --cpu-parallel
#   ./run.sh --input in.tre --output out.tre --gpu --expansion



# Check for help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "STELAR-X Run Script"
    echo ""
    echo "Usage: $0 --input <file> --output <file> [options]"
    echo ""
    echo "Options:"
    echo "  -i, --input        Gene trees file (required)"
    echo "  -o, --output       Output species tree file (required)"
    echo "  --cpu              Computation mode: CPU_SINGLE"
    echo "  --cpu-parallel      Computation mode: CPU_PARALLEL"
    echo "  --gpu              Computation mode: GPU_PARALLEL"
    echo "  -m, --mode <mode>   Computation mode: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL"
    echo "  --expansion, -e     Enable mixed bipartitions (default: OFF)"
    echo "  -v, --verbose       Verbose expansion output"
    echo ""
    echo "Scoring mode (choose one):"
    echo "  --triplet           Use STELAR-X triplet matching score (default)"
    echo "  --quartet           Use ASTRAL-style quartet matching score"
    echo "  --scoring <mode>    Scoring mode: TRIPLET or QUARTET"
    echo ""
    echo "Examples:"
    echo "  $0 --input in.tre --output out.tre"
    echo "  $0 --input in.tre --output out.tre --cpu-parallel"
    echo "  $0 --input in.tre --output out.tre --gpu --expansion"
    echo "  $0 --input in.tre --output out.tre --quartet --gpu"
    exit 0
fi

# Default configuration
DEFAULT_INPUT_FILE=""
DEFAULT_OUTPUT_FILE=""
DEFAULT_COMPUTATION_MODE="GPU_PARALLEL"
DEFAULT_VERBOSE_EXPANSION="false"
DEFAULT_USE_MIXED="false"

# Configurable values (flags override defaults)
INPUT_FILE=""
OUTPUT_FILE=""
COMPUTATION_MODE="$DEFAULT_COMPUTATION_MODE"
VERBOSE_EXPANSION="$DEFAULT_VERBOSE_EXPANSION"
USE_MIXED="$DEFAULT_USE_MIXED"
EXPANSION_ENABLED="false"
SCORING_MODE=""  # Empty means default (TRIPLET)

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input) INPUT_FILE="$2"; shift 2 ;;
        -o|--output) OUTPUT_FILE="$2"; shift 2 ;;
        --cpu) COMPUTATION_MODE="CPU_SINGLE"; shift ;;
        --cpu-parallel) COMPUTATION_MODE="CPU_PARALLEL"; shift ;;
        --gpu|--gpu-parallel) COMPUTATION_MODE="GPU_PARALLEL"; shift ;;
        -m|--mode) COMPUTATION_MODE="$2"; shift 2 ;;
        --expansion|-e) EXPANSION_ENABLED="true"; shift ;;
        -v|--verbose) VERBOSE_EXPANSION="true"; shift ;;
        --triplet) SCORING_MODE="TRIPLET"; shift ;;
        --quartet) SCORING_MODE="QUARTET"; shift ;;
        --scoring) SCORING_MODE="$2"; shift 2 ;;
        --) shift; POSITIONAL+=("$@"); break ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done

# Positional args are no longer supported
if [[ ${#POSITIONAL[@]} -gt 0 ]]; then
    echo -e "${RED}Error: Positional arguments are not supported. Use -i/--input and -o/--output.${NC}"
    exit 1
fi

# Apply defaults for input/output if still unset
INPUT_FILE=${INPUT_FILE:-$DEFAULT_INPUT_FILE}
OUTPUT_FILE=${OUTPUT_FILE:-$DEFAULT_OUTPUT_FILE}

if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
    echo -e "${RED}Error: Both --input and --output are required.${NC}"
    exit 1
fi

# Mixed bipartitions are enabled by --expansion
if [[ "$EXPANSION_ENABLED" = "true" ]]; then
    USE_MIXED="true"
else
    USE_MIXED="false"
fi

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

# Expansion method and distance method are fixed in this branch.

echo "=== STELAR-X Run Script ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Computation mode: $COMPUTATION_MODE"
echo "Expansion method: NONE (fixed)"
echo "Verbose expansion: $VERBOSE_EXPANSION"
echo "Cross-tree recombination (via --expansion): $USE_MIXED"
echo

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}"
    exit 1
fi

# Check if binaries exist
if [ ! -f "target/stelar-x-1.0.0-SNAPSHOT.jar" ] || [ ! -d "bin" ]; then
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

# Expansion is fixed to NONE; use --expansion only to enable mixed bipartitions
if [ "$EXPANSION_ENABLED" = "true" ]; then
    JAVA_ARGS="$JAVA_ARGS --expansion"
fi

# Add verbose flag if enabled
if [ "$VERBOSE_EXPANSION" = "true" ]; then
    JAVA_ARGS="$JAVA_ARGS --verbose"
fi

# Add scoring mode if specified (--triplet, --quartet, or --scoring)
if [ -n "$SCORING_MODE" ]; then
    JAVA_ARGS="$JAVA_ARGS --scoring $SCORING_MODE"
fi

# Mixed bipartitions are implied by expansion; no explicit flag needed

# Run the program with the library path set for this run only
echo -e "${YELLOW}Running STELAR-X...${NC}"
echo -e "${YELLOW}Command: java ... Main $JAVA_ARGS${NC}"
echo

eval "java -Xms4g -Xmx128g -Djava.library.path=\"$(pwd)/cuda\" \
     -Djna.debug_load=true \
     -Djna.debug_load.jna=true \
     -Djna.platform.library.path=\"$(pwd)/cuda\" \
     -Djna.memory.contiguous=true \
     -Djna.memory.contiguous.alignment=8 \
     -Djna.memory.contiguous.size=1024 \
     -cp target/stelar-x-1.0.0-SNAPSHOT.jar Main $JAVA_ARGS"

if [ $? -ne 0 ]; then
    echo -e "${RED}Program execution failed!${NC}"
    exit 1
fi

echo -e "${GREEN}Program completed successfully!${NC}"
echo "Output written to: $OUTPUT_FILE" 
