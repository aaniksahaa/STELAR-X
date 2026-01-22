#!/bin/bash

# Wrapper script to run TMC (treeFromTriplets)
#
# Usage:
#   ./run_tmc.sh -i input.tre -o output.tre
#   ./run_tmc.sh --script-dir /path/to/dir -i input.tre -o output.tre
#   ./run_tmc.sh -i input.tre -o output.tre --log output.log
#   ./run_tmc.sh -i input.tre -o output.tre --nd 1000

# Default script directory is where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT=""
OUTPUT=""
LOG_FILE=""
ND_VALUE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i) INPUT="$2"; shift 2 ;;
        -o) OUTPUT="$2"; shift 2 ;;
        --script-dir) SCRIPT_DIR="$2"; shift 2 ;;
        --log) LOG_FILE="$2"; shift 2 ;;
        --nd) ND_VALUE="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 -i <input_tree> -o <output_tree> [options]"
            echo ""
            echo "Options:"
            echo "  -i <file>           Input tree file (required)"
            echo "  -o <file>           Output tree file (required)"
            echo "  --script-dir <dir>  Directory containing treeFromTriplets (default: script location)"
            echo "  --log <file>        Log file for reconstruction details"
            echo "  --nd <number>       Limit number of triplets"
            echo "  -h, --help          Show this help message"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Validate executable exists in script directory
EXECUTABLE="$SCRIPT_DIR/treeFromTriplets"

if [[ ! -f "$EXECUTABLE" ]]; then
    echo "Error: treeFromTriplets not found in $SCRIPT_DIR"
    exit 1
fi

if [[ -z "$INPUT" ]] || [[ -z "$OUTPUT" ]]; then
    echo "Error: Both -i and -o are required"
    exit 1
fi

# Build command
CMD="$EXECUTABLE -fit $INPUT -frt $OUTPUT"

if [[ -n "$LOG_FILE" ]]; then
    CMD="$CMD -flg $LOG_FILE"
fi

if [[ -n "$ND_VALUE" ]]; then
    CMD="$CMD -nd $ND_VALUE"
fi

# Run TMC
echo "Running TMC (treeFromTriplets)..."
echo "Command: $CMD"
$CMD

if [[ $? -eq 0 ]]; then
    echo "Done: $OUTPUT"
else
    echo "Error: TMC failed"
    exit 1
fi
