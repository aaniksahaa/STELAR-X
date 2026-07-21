#!/bin/bash

# ASTRAL Runner Script
# This script runs ASTRAL using relative paths from the current directory
# NOTE: This script must be run from the ASTRAL root directory
# Usage: ./run_astral.sh -i input_file.tre [-o output_file.tre] [other_options]

# ./run_astral.sh -i ../inputs/in200.tr -o ./out.tr
# ./run_astral.sh -i ../inputs/1kp-2.tre -o ./out.tr

# Define ASTRAL root directory (current working directory)
ASTRAL_ROOT=$(pwd)

# Java memory options (override with --java-mem or env ASTRAL_JAVA_MEM)
XMS="${ASTRAL_XMS:-4g}"
XMX="${ASTRAL_XMX:-128g}"
JAVA_MEM_OPTS="-Xms${XMS} -Xmx${XMX}"
JAVA_MEM_OVERRIDE=false
JAVA_LIBRARY_PATH="${ASTRAL_JAVA_LIBRARY_PATH:-}"
ASTRAL_USE_AVX2="${ASTRAL_USE_AVX2:-false}"

# Define paths based on ASTRAL_ROOT
MAIN_DIR="${ASTRAL_ROOT}/main"
LIB_DIR="${ASTRAL_ROOT}/lib"

# Required JAR files (using absolute paths based on ASTRAL_ROOT)
MAIN_JAR="${LIB_DIR}/main.jar"
COLT_JAR="${LIB_DIR}/colt.jar"
JSAP_JAR="${LIB_DIR}/JSAP-2.1.jar"
JOCL_JAR="${LIB_DIR}/jocl-2.0.0.jar"

# Build classpath (with absolute paths, but current directory as "." for compiled classes)
CLASSPATH=".:${MAIN_JAR}:${COLT_JAR}:${JSAP_JAR}:${JOCL_JAR}"

# Check if compiled classes exist
if [ ! -f "${MAIN_DIR}/phylonet/coalescent/CommandLine.class" ]; then
    echo "Error: Compiled classes not found. Please compile first:"
    echo "cd ${MAIN_DIR}"
    echo "javac -g -classpath ${MAIN_JAR}:${COLT_JAR}:${JSAP_JAR}:${JOCL_JAR} phylonet/util/BitSet*.java phylonet/coalescent/*.java phylonet/tree/model/sti/*.java phylonet/tree/io/NewickWriter.java"
    exit 1
fi

# Check if all required JAR files exist
for jar in "${MAIN_JAR}" "${COLT_JAR}" "${JSAP_JAR}" "${JOCL_JAR}"; do
    if [ ! -f "$jar" ]; then
        echo "Error: Required JAR file not found: $jar"
        exit 1
    fi
done

# Print usage if no arguments
if [ $# -eq 0 ]; then
    echo "ASTRAL Runner Script"
    echo "==================="
    echo "Usage: $0 -i input_file.tre [-o output_file.tre] [other_options]"
    echo ""
    echo "Common options:"
    echo "  -i FILE    Input gene trees file (required)"
    echo "  -o FILE    Output species tree file (recommended)"
    echo "  -C         CPU-only mode (disable GPU)"
    echo "  -T N       Number of CPU threads"
    echo "  -G LIST    GPU indices to use (comma-separated)"
    echo "  -t N       Branch annotation level (0-4, 8, 10, 16, 32)"
    echo "  --xms SIZE       Java Xms (default: 4g or env ASTRAL_XMS)"
    echo "  --xmx SIZE       Java Xmx (default: 128g or env ASTRAL_XMX)"
    echo "  --Xms SIZE       Same as --xms"
    echo "  --Xmx SIZE       Same as --xmx"
    echo "  --java-library-path DIR  Native library path (default: lib/no_avx2)"
    echo "  --use-avx2              Use lib/libAstral.so instead of lib/no_avx2"
    echo "  --java-mem OPTS  Java memory opts (override all; legacy)"
    echo ""
    echo "Examples:"
    echo "  $0 -i inputs/in200.tr -o out.tr"
    echo "  $0 -i inputs/gene_trees.tre -o species_tree.tre -C"
    echo "  $0 -i inputs/gene_trees.tre -o species_tree.tre -T 8 -G 1"
    echo ""
    exit 0
fi

# Print configuration
echo "ASTRAL Configuration:"
echo "===================="
echo "ASTRAL_ROOT: ${ASTRAL_ROOT}"
echo "MAIN_DIR:    ${MAIN_DIR}"
echo "LIB_DIR:     ${LIB_DIR}"
echo "CLASSPATH:   ${CLASSPATH}"
echo ""

# Run ASTRAL with all provided arguments
ASTRAL_ARGS=()
OUTPUT_FILE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --java-mem)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --java-mem requires an argument"
                exit 1
            fi
            JAVA_MEM_OPTS="$2"
            JAVA_MEM_OVERRIDE=true
            shift 2
            ;;
        --xms|--Xms)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --xms requires an argument"
                exit 1
            fi
            XMS="$2"
            shift 2
            ;;
        --xmx|--Xmx)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --xmx requires an argument"
                exit 1
            fi
            XMX="$2"
            shift 2
            ;;
        --java-library-path)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --java-library-path requires an argument"
                exit 1
            fi
            JAVA_LIBRARY_PATH="$2"
            shift 2
            ;;
        --use-avx2)
            ASTRAL_USE_AVX2=true
            shift
            ;;
        -o)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -o requires an argument"
                exit 1
            fi
            OUTPUT_FILE="$2"
            ASTRAL_ARGS+=("$1" "$2")
            shift 2
            ;;
        *)
            ASTRAL_ARGS+=("$1")
            shift
            ;;
    esac
done

echo "Running ASTRAL with arguments: ${ASTRAL_ARGS[*]}"
echo "================================"

# Change to main directory to match your working command pattern
cd "${MAIN_DIR}"

# Execute ASTRAL (following your exact working pattern)
if [[ "$JAVA_MEM_OVERRIDE" = false ]]; then
    JAVA_MEM_OPTS="-Xms${XMS} -Xmx${XMX}"
fi
if [[ -z "$JAVA_LIBRARY_PATH" ]]; then
    if [[ "$ASTRAL_USE_AVX2" = true ]]; then
        JAVA_LIBRARY_PATH="$LIB_DIR"
    else
        JAVA_LIBRARY_PATH="${LIB_DIR}/no_avx2"
    fi
fi
java ${JAVA_MEM_OPTS} -Djava.library.path="${JAVA_LIBRARY_PATH}" -classpath "${CLASSPATH}" phylonet.coalescent.CommandLine "${ASTRAL_ARGS[@]}"

# Capture exit code
exit_code=$?

if [[ "$exit_code" -eq 139 || "$exit_code" -eq 134 ]]; then
    if [[ -n "$OUTPUT_FILE" && -s "$OUTPUT_FILE" ]]; then
        echo ""
        echo "Warning: ASTRAL exited with native crash code $exit_code after writing output."
        echo "Treating this as success because output exists: $OUTPUT_FILE"
        exit_code=0
    fi
fi

echo ""
echo "ASTRAL finished with exit code: $exit_code"
exit $exit_code
