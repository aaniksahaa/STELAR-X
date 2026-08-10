#!/bin/bash
#
# STELAR-X Distribution Builder
# ===============================
# Packages STELAR-X into a standalone, ready-to-run archive.
# The output archive requires ONLY Java 11+ to run — no Maven, no source code, no compilation.
#
# Usage:
#   ./dist.sh                  # Build distribution (runs install.sh first if needed)
#   ./dist.sh --skip-build     # Package existing build without rebuilding
#
# Output: dist/stelar-x-<version>.tar.gz
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Colors ──
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

VERSION="$(./project-version.sh)"
DIST_NAME="stelar-x-${VERSION}"
DIST_DIR="dist/${DIST_NAME}"
SKIP_BUILD=false

for arg in "$@"; do
  case "$arg" in
    --skip-build) SKIP_BUILD=true ;;
    -h|--help)
      echo "Usage: ./dist.sh [--skip-build]"
      echo "  --skip-build   Package existing build without rebuilding"
      exit 0
      ;;
  esac
done

echo -e "${BOLD}${CYAN}"
echo "╔══════════════════════════════════════╗"
echo "║    STELAR-X Distribution Builder     ║"
echo "╚══════════════════════════════════════╝"
echo -e "${NC}"

# ── Step 1: Ensure build exists ──
JAR_PATH="target/stelar-x-${VERSION}.jar"

if [[ "$SKIP_BUILD" == false ]]; then
  echo -e "${YELLOW}[1/4] Building STELAR-X...${NC}"
  ./install.sh --clean
  echo -e "${YELLOW}Running regression suite...${NC}"
  ./mvnw -q test
  echo ""
elif [[ ! -f "$JAR_PATH" ]]; then
  echo -e "${RED}Error: JAR not found at $JAR_PATH. Run ./install.sh or remove --skip-build.${NC}"
  exit 1
else
  echo -e "${YELLOW}[1/4] Skipping build (using existing artifacts)${NC}"
fi

# ── Step 2: Create distribution directory ──
echo -e "${YELLOW}[2/4] Assembling distribution...${NC}"

rm -rf "$DIST_DIR"
mkdir -p "$DIST_DIR/lib"

# Copy the fat JAR (rename to clean name)
cp "$JAR_PATH" "$DIST_DIR/lib/stelar-x.jar"
echo "  Copied stelar-x.jar ($(du -h "$DIST_DIR/lib/stelar-x.jar" | cut -f1))"

# Copy release metadata and user documentation
printf '%s\n' "$VERSION" > "$DIST_DIR/VERSION"
cp "README.md" "$DIST_DIR/README.md"
echo "  Added VERSION and README.md"

# Copy CUDA library if available
if [[ -f "cuda/libweight_calc.so" ]]; then
  cp "cuda/libweight_calc.so" "$DIST_DIR/lib/"
  echo "  Copied libweight_calc.so ($(du -h "$DIST_DIR/lib/libweight_calc.so" | cut -f1))"
else
  echo -e "  ${YELLOW}No CUDA library found — distribution will be CPU-only${NC}"
fi

# Copy example data
mkdir -p "$DIST_DIR/examples"
if [[ -f "all_gt_bs_rooted_37.tre" ]]; then
  cp "all_gt_bs_rooted_37.tre" "$DIST_DIR/examples/"
  echo "  Copied examples/all_gt_bs_rooted_37.tre ($(du -h "$DIST_DIR/examples/all_gt_bs_rooted_37.tre" | cut -f1))"
else
  echo -e "  ${YELLOW}Warning: all_gt_bs_rooted_37.tre not found — examples/ will be empty${NC}"
fi

# ── Step 3: Create the standalone launcher ──
echo -e "${YELLOW}[3/4] Creating launcher scripts...${NC}"

cat > "$DIST_DIR/stelar-x" << 'LAUNCHER_EOF'
#!/bin/bash
#
# STELAR-X — Species Tree Estimation using Triplet Frequencies
# ==============================================================
#
# Usage:
#   stelar-x -i <gene_trees> -o <output> [options]
#
# Options:
#   -i, --input <file>     Gene trees file (required)
#   -o, --output <file>    Output species tree file (required)
#   -c, --score <tree>     Score-only mode: score gene trees against given species tree
#   --cpu                  CPU single-threaded mode
#   --cpu-parallel         CPU parallel mode
#   --gpu                  GPU parallel mode (default if GPU available)
#   -m, --mode <mode>      CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL
#   -e, --expansion        Enable mixed bipartitions
#   -s, --support <type>   Branch support: NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL
#   --lambda <val>         Lambda parameter for branch support (default: 0.5)
#   -v, --verbose          Verbose expansion output
#   --xms <size>           Java min heap (default: 4g)
#   --xmx <size>           Java max heap (default: 128g)
#   --version              Print the installed STELAR-X version
#   -h, --help             Show this message
#
# Examples:
#   stelar-x -i genes.tre -o species.tre
#   stelar-x -i genes.tre -o species.tre --cpu-parallel
#   stelar-x -i genes.tre -o species.tre --gpu --expansion
#   stelar-x -i genes.tre -c known_species.tre

# ── Resolve install location (works via symlink too) ──
SOURCE="${BASH_SOURCE[0]}"
while [[ -L "$SOURCE" ]]; do
  DIR="$(cd "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ "$SOURCE" != /* ]] && SOURCE="$DIR/$SOURCE"
done
STELAR_HOME="$(cd "$(dirname "$SOURCE")" && pwd)"

# ── Paths ──
JAR="$STELAR_HOME/lib/stelar-x.jar"
LIB_DIR="$STELAR_HOME/lib"

if [[ ! -f "$JAR" ]]; then
  echo "Error: stelar-x.jar not found at $JAR"
  echo "The installation appears to be incomplete."
  exit 1
fi

if [[ "${1:-}" == "--version" ]]; then
  printf 'STELAR-X %s\n' "$(cat "$STELAR_HOME/VERSION")"
  exit 0
fi

# ── Colors ──
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# ── Help ──
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '3,/^$/{ s/^# \{0,1\}//; p }' "${BASH_SOURCE[0]}"
  exit 0
fi

# ── Defaults ──
XMS="${STELAR_XMS:-4g}"
XMX="${STELAR_XMX:-128g}"
HAS_CUDA_LIB=false
[[ -f "$LIB_DIR/libweight_calc.so" ]] && HAS_CUDA_LIB=true

# ── Parse args (separate Java heap opts from program args) ──
JAVA_PROGRAM_ARGS=()
MODE_SET=false
INPUT_FILE=""
OUTPUT_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)       INPUT_FILE="$2"; JAVA_PROGRAM_ARGS+=("-i" "$2"); shift 2 ;;
    -o|--output)      OUTPUT_FILE="$2"; JAVA_PROGRAM_ARGS+=("-o" "$2"); shift 2 ;;
    -c|--score|--species-tree) JAVA_PROGRAM_ARGS+=("-c" "$2"); shift 2 ;;
    --cpu)            MODE_SET=true; JAVA_PROGRAM_ARGS+=("--cpu"); shift ;;
    --cpu-parallel)   MODE_SET=true; JAVA_PROGRAM_ARGS+=("--cpu-parallel"); shift ;;
    --gpu|--gpu-parallel) MODE_SET=true; JAVA_PROGRAM_ARGS+=("--gpu"); shift ;;
    -m|--mode)        MODE_SET=true; JAVA_PROGRAM_ARGS+=("-m" "$2"); shift 2 ;;
    -e|--expansion)   JAVA_PROGRAM_ARGS+=("--expansion"); shift ;;
    -v|--verbose)     JAVA_PROGRAM_ARGS+=("--verbose"); shift ;;
    -s|--support|--branch-support) JAVA_PROGRAM_ARGS+=("-s" "$2"); shift 2 ;;
    --lambda)         JAVA_PROGRAM_ARGS+=("--lambda" "$2"); shift 2 ;;
    --xms|--Xms)      XMS="$2"; shift 2 ;;
    --xmx|--Xmx)      XMX="$2"; shift 2 ;;
    *)
      echo -e "${RED}Error: Unknown option '$1'. Run stelar-x --help for usage.${NC}"
      exit 1
      ;;
  esac
done

# ── Validate ──
if [[ -z "$INPUT_FILE" ]]; then
  echo -e "${RED}Error: --input is required. Run stelar-x --help for usage.${NC}"
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo -e "${RED}Error: Input file '$INPUT_FILE' does not exist.${NC}"
  exit 1
fi

# ── Auto-detect mode if not set ──
if [[ "$MODE_SET" == false ]]; then
  if [[ "$HAS_CUDA_LIB" == true ]] && command -v nvidia-smi &>/dev/null; then
    JAVA_PROGRAM_ARGS+=("-m" "GPU_PARALLEL")
  else
    JAVA_PROGRAM_ARGS+=("-m" "CPU_PARALLEL")
  fi
fi

# ── Create output directory if needed ──
if [[ -n "$OUTPUT_FILE" ]]; then
  mkdir -p "$(dirname "$OUTPUT_FILE")" 2>/dev/null || true
fi

# ── Run ──
exec java \
  -Xms"${XMS}" -Xmx"${XMX}" \
  -Djava.library.path="$LIB_DIR" \
  -Djna.platform.library.path="$LIB_DIR" \
  -cp "$JAR" \
  Main "${JAVA_PROGRAM_ARGS[@]}"
LAUNCHER_EOF

chmod +x "$DIST_DIR/stelar-x"
echo "  Created stelar-x launcher"

# ── Step 4: Package ──
echo -e "${YELLOW}[4/4] Creating archive...${NC}"

(cd dist && tar czf "${DIST_NAME}.tar.gz" "${DIST_NAME}")
(cd dist && sha256sum "${DIST_NAME}.tar.gz" > "${DIST_NAME}.tar.gz.sha256")

ARCHIVE="dist/${DIST_NAME}.tar.gz"
CHECKSUM="${ARCHIVE}.sha256"
ARCHIVE_SIZE=$(du -h "$ARCHIVE" | cut -f1)

echo ""
echo -e "${BOLD}${GREEN}"
echo "╔══════════════════════════════════════╗"
echo "║   Distribution Ready!                ║"
echo "╚══════════════════════════════════════╝"
echo -e "${NC}"

echo -e "  Archive:  ${GREEN}${ARCHIVE}${NC} (${ARCHIVE_SIZE})"
echo -e "  SHA-256: ${GREEN}${CHECKSUM}${NC}"
echo ""
echo -e "  Contents:"
echo "    stelar-x/stelar-x                          (launcher script)"
echo "    stelar-x/VERSION                           (release version)"
echo "    stelar-x/README.md                         (documentation)"
echo "    stelar-x/lib/stelar-x.jar                  (Java application)"
if [[ -f "$DIST_DIR/lib/libweight_calc.so" ]]; then
  echo "    stelar-x/lib/libweight_calc.so             (CUDA GPU library)"
fi
if [[ -f "$DIST_DIR/examples/all_gt_bs_rooted_37.tre" ]]; then
  echo "    stelar-x/examples/all_gt_bs_rooted_37.tre  (example gene trees)"
fi
echo ""
echo -e "${BOLD}For end users — install & run:${NC}"
echo ""
echo "  tar xzf ${DIST_NAME}.tar.gz"
echo "  cd ${DIST_NAME}"
echo "  ./stelar-x -i examples/all_gt_bs_rooted_37.tre -o examples/out_37.tre"
echo ""
echo -e "${BOLD}Optional — add to PATH:${NC}"
echo ""
echo "  sudo ln -sf \$(pwd)/stelar-x /usr/local/bin/stelar-x"
echo "  stelar-x -i gene_trees.tre -o output.tre"
