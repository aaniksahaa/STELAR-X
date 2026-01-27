#!/bin/bash

# STELAR-X Run Script with Bipartition Expansion Support
# Usage: ./run.sh -i <input_file> [options]
#
# For detailed help: ./run.sh -h

set -euo pipefail

# Default configuration
DEFAULT_OUTPUT_FILE="out.tre"
DEFAULT_EXPANSION_METHOD="NONE"
DEFAULT_DISTANCE_METHOD="UPGMA"
DEFAULT_XMS="4g"
DEFAULT_XMX="128g"

# Initialize variables with defaults
INPUT_FILE=""
OUTPUT_FILE="$DEFAULT_OUTPUT_FILE"
NUM_THREADS=""  # Empty means use all available
CPU_ONLY=false  # GPU is default, --cpu disables it
EXPANSION_METHOD="$DEFAULT_EXPANSION_METHOD"
DISTANCE_METHOD="$DEFAULT_DISTANCE_METHOD"
VERBOSE_EXPANSION=false
XMS="${STELAR_XMS:-$DEFAULT_XMS}"
XMX="${STELAR_XMX:-$DEFAULT_XMX}"
JAVA_MEM_OPTS="-Xms${XMS} -Xmx${XMX}"
JAVA_MEM_OVERRIDE=false

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_help() {
    cat <<EOF
STELAR-X Run Script
A species tree inference tool with bipartition expansion support.

Usage:
    $0 -i <input_file> [options]
    $0 --help

Required:
    -i, --input FILE        Input gene trees file (Newick format)

Options:
    -o, --output FILE       Output species tree file (default: out.tre)
    --threads NUM           Number of threads to use (default: all available)
    --cpu                   Use CPU only (GPU acceleration is enabled by default)
    --xms SIZE              Java Xms (default: ${DEFAULT_XMS} or env STELAR_XMS)
    --xmx SIZE              Java Xmx (default: ${DEFAULT_XMX} or env STELAR_XMX)
    --Xms SIZE              Same as --xms
    --Xmx SIZE              Same as --xmx
    --java-mem OPTS         Java memory opts (override all; legacy)

Expansion Options:
    -e, --expansion METHOD  Bipartition expansion method (default: NONE)
                            Options: NONE, DISTANCE_ONLY, CONSENSUS_ONLY,
                                     DISTANCE_CONSENSUS, FULL
    -d, --distance METHOD   Distance calculation method (default: UPGMA)
                            Options: UPGMA, NEIGHBOR_JOINING, BOTH
    -v, --verbose           Enable verbose expansion output

General:
    -h, --help              Show this help message and exit
    --version               Show version information

Bipartition Expansion Methods:
    NONE              No expansion, original STELAR behavior
    DISTANCE_ONLY     Only distance-based expansion (UPGMA trees)
    CONSENSUS_ONLY    Only consensus-based expansion (greedy consensus)
    DISTANCE_CONSENSUS  Both distance and consensus expansion (recommended)
    FULL              All expansion methods including polytomy resolution

Examples:
    # Basic usage with GPU acceleration (default)
    $0 -i genes.tre -o species.tre

    # Use 4 threads with GPU acceleration
    $0 -i genes.tre -o species.tre --threads 4

    # Use CPU only (no GPU)
    $0 -i genes.tre -o species.tre --cpu

    # CPU only with 4 threads
    $0 -i genes.tre -o species.tre --cpu --threads 4

    # With distance expansion
    $0 -i genes.tre -o species.tre -e DISTANCE_ONLY

    # Full expansion with verbose output
    $0 -i genes.tre -o species.tre -e FULL -v

EOF
}

print_version() {
    echo "STELAR-X version 1.0.0"
    echo "Species tree inference with bipartition expansion support"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            if [[ -z "${2:-}" ]]; then
                echo -e "${RED}Error: --input requires a file argument${NC}" >&2
                exit 1
            fi
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            if [[ -z "${2:-}" ]]; then
                echo -e "${RED}Error: --output requires a file argument${NC}" >&2
                exit 1
            fi
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --threads)
            if [[ -z "${2:-}" ]]; then
                echo -e "${RED}Error: --threads requires a number argument${NC}" >&2
                exit 1
            fi
            if ! [[ "$2" =~ ^[0-9]+$ ]] || [[ "$2" -lt 1 ]]; then
                echo -e "${RED}Error: --threads must be a positive integer${NC}" >&2
                exit 1
            fi
            NUM_THREADS="$2"
            shift 2
            ;;
        --cpu)
            CPU_ONLY=true
            shift
            ;;
        --java-mem)
            if [[ -z "${2:-}" ]]; then
                echo -e "${RED}Error: --java-mem requires an argument${NC}" >&2
                exit 1
            fi
            JAVA_MEM_OPTS="$2"
            JAVA_MEM_OVERRIDE=true
            shift 2
            ;;
        --xms|--Xms)
            if [[ -z "${2:-}" ]]; then
                echo -e "${RED}Error: --xms requires an argument${NC}" >&2
                exit 1
            fi
            XMS="$2"
            shift 2
            ;;
        --xmx|--Xmx)
            if [[ -z "${2:-}" ]]; then
                echo -e "${RED}Error: --xmx requires an argument${NC}" >&2
                exit 1
            fi
            XMX="$2"
            shift 2
            ;;
        -e|--expansion)
            if [[ -z "${2:-}" ]]; then
                echo -e "${RED}Error: --expansion requires a method argument${NC}" >&2
                exit 1
            fi
            EXPANSION_METHOD="$2"
            shift 2
            ;;
        -d|--distance)
            if [[ -z "${2:-}" ]]; then
                echo -e "${RED}Error: --distance requires a method argument${NC}" >&2
                exit 1
            fi
            DISTANCE_METHOD="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE_EXPANSION=true
            shift
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        --version)
            print_version
            exit 0
            ;;
        -*)
            echo -e "${RED}Error: Unknown option: $1${NC}" >&2
            echo "Use '$0 --help' for usage information."
            exit 1
            ;;
        *)
            # Handle legacy positional arguments for backward compatibility
            if [[ -z "$INPUT_FILE" ]]; then
                INPUT_FILE="$1"
            elif [[ "$OUTPUT_FILE" == "$DEFAULT_OUTPUT_FILE" ]]; then
                OUTPUT_FILE="$1"
            else
                echo -e "${RED}Error: Unexpected argument: $1${NC}" >&2
                exit 1
            fi
            shift
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_FILE" ]]; then
    echo -e "${RED}Error: Input file is required.${NC}" >&2
    echo "Usage: $0 -i <input_file> [options]"
    echo "Use '$0 --help' for detailed usage information."
    exit 1
fi

# Validate expansion method
VALID_EXPANSION_METHODS=("NONE" "DISTANCE_ONLY" "CONSENSUS_ONLY" "DISTANCE_CONSENSUS" "FULL")
EXPANSION_METHOD=$(echo "$EXPANSION_METHOD" | tr '[:lower:]' '[:upper:]')
if [[ ! " ${VALID_EXPANSION_METHODS[*]} " =~ " ${EXPANSION_METHOD} " ]]; then
    echo -e "${RED}Error: Invalid expansion method '$EXPANSION_METHOD'${NC}" >&2
    echo "Valid expansion methods: ${VALID_EXPANSION_METHODS[*]}"
    exit 1
fi

# Validate distance method
VALID_DISTANCE_METHODS=("UPGMA" "NEIGHBOR_JOINING" "BOTH")
DISTANCE_METHOD=$(echo "$DISTANCE_METHOD" | tr '[:lower:]' '[:upper:]')
if [[ ! " ${VALID_DISTANCE_METHODS[*]} " =~ " ${DISTANCE_METHOD} " ]]; then
    echo -e "${RED}Error: Invalid distance method '$DISTANCE_METHOD'${NC}" >&2
    echo "Valid distance methods: ${VALID_DISTANCE_METHODS[*]}"
    exit 1
fi

echo "=== STELAR-X Run Script ==="
echo "Input file:       $INPUT_FILE"
echo "Output file:      $OUTPUT_FILE"
echo "Threads:          ${NUM_THREADS:-auto}"
echo "Mode:             $(if [[ "$CPU_ONLY" == true ]]; then echo "CPU only"; else echo "GPU (with CPU fallback)"; fi)"
echo "Expansion method: $EXPANSION_METHOD"
echo "Distance method:  $DISTANCE_METHOD"
echo "Verbose:          $VERBOSE_EXPANSION"
if [[ "$JAVA_MEM_OVERRIDE" = true ]]; then
    echo "Java mem opts:    $JAVA_MEM_OPTS (override)"
else
    JAVA_MEM_OPTS="-Xms${XMS} -Xmx${XMX}"
    echo "Java mem opts:    $JAVA_MEM_OPTS"
fi
echo

# Check if input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}" >&2
    exit 1
fi

# Check if binaries exist
if [[ ! -f "target/stelar-x-1.0.0-SNAPSHOT.jar" ]] || [[ ! -d "bin" ]]; then
    echo -e "${RED}Error: Binaries not found. Please run build.sh first.${NC}" >&2
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
if [[ "$OUTPUT_DIR" != "." ]] && [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Debug information
echo -e "${YELLOW}Debug Information:${NC}"
echo "Current directory: $(pwd)"
echo "Library path: $(pwd)/cuda"
echo "Library exists: $(if [[ -f "$(pwd)/cuda/libweight_calc.so" ]]; then echo "Yes"; else echo "No"; fi)"
echo "Library permissions: $(ls -l "$(pwd)/cuda/libweight_calc.so" 2>/dev/null || echo "Not found")"
echo

# Build command line arguments for Java application
JAVA_ARGS="-i \"$INPUT_FILE\" -o \"$OUTPUT_FILE\""

# Add thread count if specified
if [[ -n "$NUM_THREADS" ]]; then
    JAVA_ARGS="$JAVA_ARGS --threads $NUM_THREADS"
fi

# Add CPU-only flag if specified (GPU is default)
if [[ "$CPU_ONLY" == true ]]; then
    JAVA_ARGS="$JAVA_ARGS --cpu"
fi

# Add expansion method
if [[ "$EXPANSION_METHOD" != "NONE" ]]; then
    JAVA_ARGS="$JAVA_ARGS -e \"$EXPANSION_METHOD\""
else
    JAVA_ARGS="$JAVA_ARGS --no-expansion"
fi

# Add distance method if not default
if [[ "$DISTANCE_METHOD" != "UPGMA" ]]; then
    JAVA_ARGS="$JAVA_ARGS -d \"$DISTANCE_METHOD\""
fi

# Add verbose flag if enabled
if [[ "$VERBOSE_EXPANSION" == true ]]; then
    JAVA_ARGS="$JAVA_ARGS -v"
fi

# Run the program with the library path set for this run only
echo -e "${YELLOW}Running STELAR-X...${NC}"
echo -e "${YELLOW}Command: java ... Main $JAVA_ARGS${NC}"
echo

eval "java ${JAVA_MEM_OPTS} -Djava.library.path=\"$(pwd)/cuda\" \
     -Djna.debug_load=true \
     -Djna.debug_load.jna=true \
     -Djna.platform.library.path=\"$(pwd)/cuda\" \
     -Djna.memory.contiguous=true \
     -Djna.memory.contiguous.alignment=8 \
     -Djna.memory.contiguous.size=1024 \
     -cp target/stelar-x-1.0.0-SNAPSHOT.jar Main $JAVA_ARGS"

if [[ $? -ne 0 ]]; then
    echo -e "${RED}Program execution failed!${NC}" >&2
    exit 1
fi

echo -e "${GREEN}Program completed successfully!${NC}"
echo "Output written to: $OUTPUT_FILE"
