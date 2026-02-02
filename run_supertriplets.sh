#!/bin/bash

# Wrapper script to run SuperTriplets and optionally clean branch support values
#
# Usage:
#   ./run_supertriplets.sh -i input.tre -o output.tre
#   ./run_supertriplets.sh -i input.tre -o output.tre --branch-support
#   ./run_supertriplets.sh --script-dir /path/to/dir -i input.tre -o output.tre

# Default script directory is where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT=""
OUTPUT=""
BRANCH_SUPPORT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i) INPUT="$2"; shift 2 ;;
        -o) OUTPUT="$2"; shift 2 ;;
        --script-dir) SCRIPT_DIR="$2"; shift 2 ;;
        --branch-support) BRANCH_SUPPORT="--branch-support"; shift ;;
        -h|--help)
            echo "Usage: $0 -i <input_tree> -o <output_tree> [--script-dir <dir>] [--branch-support]"
            echo "  --script-dir      Directory containing SuperTriplets_v1.1.jar and clean_tree.py (default: script location)"
            echo "  --branch-support  Keep branch support values (default: remove them)"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Validate required files exist in script directory
JAR_FILE="$SCRIPT_DIR/SuperTriplets_v1.1.jar"
PY_FILE="$SCRIPT_DIR/clean_tree.py"

if [[ ! -f "$JAR_FILE" ]]; then
    echo "Error: SuperTriplets_v1.1.jar not found in $SCRIPT_DIR"
    exit 1
fi

if [[ ! -f "$PY_FILE" ]]; then
    echo "Error: clean_tree.py not found in $SCRIPT_DIR"
    exit 1
fi

if [[ -z "$INPUT" ]] || [[ -z "$OUTPUT" ]]; then
    echo "Error: Both -i and -o are required"
    exit 1
fi

# Temp file for SuperTriplets output
TEMP_OUTPUT=$(mktemp)
TEMP_STDOUT=$(mktemp)

# Run SuperTriplets and capture stdout
echo "Running SuperTriplets..."
java -jar -Xmx16g "$JAR_FILE" "$INPUT" "$TEMP_OUTPUT" | tee "$TEMP_STDOUT"

# Extract the NNI tree (critNNI line) and save to .nni file
NNI_OUTPUT="${OUTPUT}.nni"
NNI_TREE=$(grep "^critNNI@" "$TEMP_STDOUT" | sed 's/^critNNI@[^@]*@//')
if [[ -n "$NNI_TREE" ]]; then
    echo "$NNI_TREE" > "$TEMP_OUTPUT.nni"
    echo "Cleaning NNI tree..."
    python "$PY_FILE" "$TEMP_OUTPUT.nni" "$NNI_OUTPUT" $BRANCH_SUPPORT
    rm -f "$TEMP_OUTPUT.nni"
    echo "NNI tree saved to: $NNI_OUTPUT"
fi

# Clean the tree (SNNI tree from output file)
echo "Cleaning SNNI tree..."
python "$PY_FILE" "$TEMP_OUTPUT" "$OUTPUT" $BRANCH_SUPPORT

# Remove temp files
rm -f "$TEMP_OUTPUT" "$TEMP_STDOUT"

echo "Done: $OUTPUT"
