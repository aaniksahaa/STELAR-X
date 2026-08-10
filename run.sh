#!/bin/bash
#
# STELAR-X Runner
# ================
# Usage: ./run.sh -i <gene_trees> -o <output> [options]
#
# Options:
#   -i, --input <file>     Gene trees file (required)
#   -o, --output <file>    Output species tree file (required)
#   -c, --score <tree>     Score-only mode: score gene trees against given species tree
#   --cpu                  CPU single-threaded mode
#   --cpu-parallel         CPU parallel mode
#   --gpu                  GPU parallel mode (default)
#   -m, --mode <mode>      CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL
#   -T, --threads <n>      Max CPU threads/cores to expose to Java
#   -e, --expansion        Enable mixed bipartitions
#   -s, --support <type>   Branch support: NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL
#   --lambda <val>         Lambda parameter for branch support (default: 0.5)
#   -v, --verbose          Verbose expansion output
#   --xms <size>           Java min heap (default: 4g)
#   --xmx <size>           Java max heap (default: 128g)
#   -h, --help             Show this message
#
# Examples:
#   ./run.sh -i genes.tre -o species.tre
#   ./run.sh -i genes.tre -o species.tre --cpu-parallel
#   ./run.sh -i genes.tre -o species.tre --gpu --expansion
#   ./run.sh -i genes.tre -c known_species.tre    # score-only mode

# ── Resolve STELAR_ROOT (works from any directory, even via symlink) ──
STELAR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Colors ──
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# ── Defaults ──
INPUT_FILE=""
OUTPUT_FILE=""
COMPUTATION_MODE=""
VERBOSE_EXPANSION="false"
EXPANSION_ENABLED="false"
XMS="${STELAR_XMS:-4g}"
XMX="${STELAR_XMX:-128g}"
EXTRA_JAVA_ARGS=()
THREADS_REQUESTED=""

# ── Help ──
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '3,/^$/{ s/^# \{0,1\}//; p }' "${BASH_SOURCE[0]}"
  exit 0
fi

# ── Parse args ──
JAVA_PROGRAM_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)       INPUT_FILE="$2"; JAVA_PROGRAM_ARGS+=("-i" "$2"); shift 2 ;;
    -o|--output)      OUTPUT_FILE="$2"; JAVA_PROGRAM_ARGS+=("-o" "$2"); shift 2 ;;
    -c|--score|--species-tree) JAVA_PROGRAM_ARGS+=("-c" "$2"); shift 2 ;;
    --cpu)            COMPUTATION_MODE="CPU_SINGLE"; JAVA_PROGRAM_ARGS+=("--cpu"); shift ;;
    --cpu-parallel)   COMPUTATION_MODE="CPU_PARALLEL"; JAVA_PROGRAM_ARGS+=("--cpu-parallel"); shift ;;
    --gpu|--gpu-parallel) COMPUTATION_MODE="GPU_PARALLEL"; JAVA_PROGRAM_ARGS+=("--gpu"); shift ;;
    -m|--mode)        COMPUTATION_MODE="$2"; JAVA_PROGRAM_ARGS+=("-m" "$2"); shift 2 ;;
    -T|--threads)
      THREADS_REQUESTED="$2"
      if [[ ! "$THREADS_REQUESTED" =~ ^[1-9][0-9]*$ ]]; then
        echo -e "${RED}Error: --threads requires a positive integer.${NC}"
        exit 1
      fi
      shift 2
      ;;
    -e|--expansion)   EXPANSION_ENABLED="true"; JAVA_PROGRAM_ARGS+=("--expansion"); shift ;;
    -v|--verbose)     VERBOSE_EXPANSION="true"; JAVA_PROGRAM_ARGS+=("--verbose"); shift ;;
    -s|--support|--branch-support) JAVA_PROGRAM_ARGS+=("-s" "$2"); shift 2 ;;
    --lambda)         JAVA_PROGRAM_ARGS+=("--lambda" "$2"); shift 2 ;;
    --xms|--Xms)      XMS="$2"; shift 2 ;;
    --xmx|--Xmx)      XMX="$2"; shift 2 ;;
    *)
      echo -e "${RED}Error: Unknown option '$1'. Run ./run.sh --help for usage.${NC}"
      exit 1
      ;;
  esac
done

# ── Validate ──
if [[ -z "$INPUT_FILE" ]]; then
  echo -e "${RED}Error: --input is required. Run ./run.sh --help for usage.${NC}"
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}"
  exit 1
fi

# ── Find JAR ──
VERSION="$("$STELAR_ROOT/project-version.sh")"
JAR="$STELAR_ROOT/target/stelar-x-${VERSION}.jar"
if [[ ! -f "$JAR" ]]; then
  echo -e "${RED}Error: JAR not found at $JAR${NC}"
  echo -e "Run ${YELLOW}./install.sh${NC} first to build STELAR-X."
  exit 1
fi

# ── Determine CUDA library path ──
CUDA_LIB_DIR="$STELAR_ROOT/cuda"
HAS_CUDA_LIB=false
if [[ -f "$CUDA_LIB_DIR/libweight_calc.so" ]]; then
  HAS_CUDA_LIB=true
fi

# ── Auto-detect mode if not specified ──
if [[ -z "$COMPUTATION_MODE" ]]; then
  if [[ "$HAS_CUDA_LIB" == true ]] && command -v nvidia-smi &>/dev/null; then
    COMPUTATION_MODE="GPU_PARALLEL"
  else
    COMPUTATION_MODE="CPU_PARALLEL"
  fi
  JAVA_PROGRAM_ARGS+=("-m" "$COMPUTATION_MODE")
fi

# Warn if GPU requested but no library
if [[ "$COMPUTATION_MODE" == "GPU_PARALLEL" && "$HAS_CUDA_LIB" == false ]]; then
  echo -e "${YELLOW}Warning: GPU mode requested but CUDA library not found at $CUDA_LIB_DIR/libweight_calc.so${NC}"
  echo -e "${YELLOW}The program will attempt GPU and fall back to CPU if it fails.${NC}"
fi

# ── Cap Java-visible CPU count if requested ──
if [[ -n "$THREADS_REQUESTED" ]]; then
  AVAILABLE_THREADS="$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo "$THREADS_REQUESTED")"
  if [[ ! "$AVAILABLE_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    AVAILABLE_THREADS="$THREADS_REQUESTED"
  fi
  ACTIVE_THREADS="$THREADS_REQUESTED"
  if (( ACTIVE_THREADS > AVAILABLE_THREADS )); then
    ACTIVE_THREADS="$AVAILABLE_THREADS"
  fi
  EXTRA_JAVA_ARGS+=("-XX:ActiveProcessorCount=${ACTIVE_THREADS}")
fi

# ── Create output directory if needed ──
if [[ -n "$OUTPUT_FILE" ]]; then
  OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
  if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
  fi
fi

# ── Print summary ──
echo "=== STELAR-X ==="
echo "Input:       $INPUT_FILE"
[[ -n "$OUTPUT_FILE" ]] && echo "Output:      $OUTPUT_FILE"
echo "Mode:        $COMPUTATION_MODE"
if [[ -n "$THREADS_REQUESTED" ]]; then
  echo "Threads:     ${ACTIVE_THREADS} (requested ${THREADS_REQUESTED})"
fi
echo "Expansion:   $EXPANSION_ENABLED"
echo "Java heap:   -Xms${XMS} -Xmx${XMX}"
echo

# ── Run ──
exec java \
  -Xms"${XMS}" -Xmx"${XMX}" \
  "${EXTRA_JAVA_ARGS[@]}" \
  -Djava.library.path="$CUDA_LIB_DIR" \
  -Djna.platform.library.path="$CUDA_LIB_DIR" \
  -cp "$JAR" \
  Main "${JAVA_PROGRAM_ARGS[@]}"
