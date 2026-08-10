#!/bin/bash
#
# STELAR-X Installer
# ===================
# One-command build for STELAR-X.
#
# Prerequisites: Java 11+ (JDK)
# Optional:      CUDA toolkit (nvcc) for GPU acceleration
#
# Usage:
#   ./install.sh              # Build everything (auto-detects CUDA)
#   ./install.sh --no-cuda    # Skip CUDA compilation entirely
#   ./install.sh --force-cuda # Fail if CUDA compilation fails
#   ./install.sh --clean      # Clean build (remove previous artifacts first)
#
set -euo pipefail

# ── Resolve script directory (works even when called via symlink) ──
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Colors ──
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# ── Defaults ──
CUDA_MODE="auto"   # auto | skip | required
CLEAN=false
QUICK=false        # skip Java checks (for fast rebuilds via build.sh)

# ── Parse args ──
for arg in "$@"; do
  case "$arg" in
    --no-cuda)    CUDA_MODE="skip" ;;
    --force-cuda) CUDA_MODE="required" ;;
    --clean)      CLEAN=true ;;
    --quick)      QUICK=true ;;
    -h|--help)
      echo "STELAR-X Installer / Builder"
      echo ""
      echo "Usage: ./install.sh [options]"
      echo "       ./build.sh [options]     (same thing, skips Java checks by default)"
      echo ""
      echo "Options:"
      echo "  --no-cuda      Skip CUDA compilation (CPU-only mode)"
      echo "  --force-cuda   Fail if CUDA cannot be compiled"
      echo "  --clean        Clean previous build artifacts first"
      echo "  --quick        Skip Java/JDK checks (fast rebuild)"
      echo "  -h, --help     Show this message"
      echo ""
      echo "Prerequisites:"
      echo "  - Java 11+ (JDK)   Required"
      echo "  - CUDA toolkit     Optional (for GPU acceleration)"
      exit 0
      ;;
    *)
      echo -e "${RED}Unknown option: $arg${NC}"
      echo "Run ./install.sh --help for usage."
      exit 1
      ;;
  esac
done

if [[ "$QUICK" == true ]]; then
  echo -e "${BOLD}${CYAN}"
  echo "╔══════════════════════════════════════╗"
  echo "║        STELAR-X Build                ║"
  echo "╚══════════════════════════════════════╝"
  echo -e "${NC}"
else
  echo -e "${BOLD}${CYAN}"
  echo "╔══════════════════════════════════════╗"
  echo "║        STELAR-X Installer            ║"
  echo "╚══════════════════════════════════════╝"
  echo -e "${NC}"

  # ── Step 1: Check Java ──
  echo -e "${YELLOW}[1/3] Checking Java...${NC}"

  if ! command -v java &>/dev/null; then
    echo -e "${RED}Error: Java is not installed.${NC}"
    echo ""
    echo "Please install Java 11+ (JDK):"
    echo "  Ubuntu/Debian:  sudo apt install -y openjdk-17-jdk"
    echo "  Fedora/RHEL:    sudo dnf install -y java-17-openjdk-devel"
    echo "  macOS:          brew install openjdk@17"
    echo "  Arch:           sudo pacman -S jdk17-openjdk"
    echo ""
    echo "Then re-run: ./install.sh"
    exit 1
  fi

  # Check Java version (need 11+)
  JAVA_VER=$(java -version 2>&1 | head -n1 | sed -E 's/.*"([0-9]+)\..*/\1/' | head -c 10)
  if [[ "$JAVA_VER" =~ ^[0-9]+$ ]] && (( JAVA_VER < 11 )); then
    echo -e "${RED}Error: Java $JAVA_VER detected, but Java 11+ is required.${NC}"
    echo "Please upgrade Java and re-run: ./install.sh"
    exit 1
  fi

  # Check for javac (JDK, not just JRE)
  if ! command -v javac &>/dev/null; then
    echo -e "${RED}Error: javac not found. You have Java runtime but not the JDK.${NC}"
    echo "Please install the full JDK (not just JRE)."
    echo "  Ubuntu/Debian:  sudo apt install -y openjdk-17-jdk"
    exit 1
  fi

  echo -e "  ${GREEN}Java $JAVA_VER found${NC}"
fi

# ── Build Java fat JAR ──
STEP_PREFIX=""
if [[ "$QUICK" == false ]]; then STEP_PREFIX="[2/3] "; fi
echo -e "\n${YELLOW}${STEP_PREFIX}Building STELAR-X JAR...${NC}"

if [[ "$CLEAN" == true ]]; then
  echo "  Cleaning previous build..."
  ./mvnw -q clean 2>/dev/null || true
fi

# Use Maven Wrapper (no Maven installation needed!)
echo "  Compiling with Maven Wrapper (no Maven installation needed)..."
./mvnw -q package -DskipTests 2>&1 | grep -v "^\[INFO\]" | grep -v "^Progress" || true

# Verify JAR was built
VERSION="$(./project-version.sh)"
JAR_PATH="target/stelar-x-${VERSION}.jar"
if [[ ! -f "$JAR_PATH" ]]; then
  echo -e "${RED}Error: JAR build failed. $JAR_PATH not found.${NC}"
  echo "Try running: ./mvnw package  (for detailed output)"
  exit 1
fi

JAR_SIZE=$(du -h "$JAR_PATH" | cut -f1)
echo -e "  ${GREEN}JAR built successfully ($JAR_SIZE)${NC}"

# ── Build CUDA library ──
STEP_PREFIX=""
if [[ "$QUICK" == false ]]; then STEP_PREFIX="[3/3] "; fi
echo -e "\n${YELLOW}${STEP_PREFIX}Building CUDA library...${NC}"

CUDA_SO="cuda/libweight_calc.so"
CUDA_STATUS="skipped"

if [[ "$CUDA_MODE" == "skip" ]]; then
  echo "  Skipping CUDA compilation (--no-cuda)"
  if [[ -f "$CUDA_SO" ]]; then
    echo -e "  ${GREEN}Pre-built CUDA library found — GPU mode available${NC}"
    CUDA_STATUS="pre-built"
  else
    echo -e "  ${YELLOW}No CUDA library — GPU mode will NOT be available${NC}"
    CUDA_STATUS="none"
  fi
elif command -v nvcc &>/dev/null; then
  echo "  nvcc found: $(nvcc --version 2>&1 | grep 'release' | sed 's/.*release //' | sed 's/,.*//')"
  echo "  Compiling CUDA kernels..."
  if (cd cuda && make clean >/dev/null 2>&1; make 2>&1 | tail -3); then
    echo -e "  ${GREEN}CUDA library compiled successfully${NC}"
    CUDA_STATUS="compiled"
  else
    if [[ "$CUDA_MODE" == "required" ]]; then
      echo -e "${RED}Error: CUDA compilation failed (--force-cuda was set)${NC}"
      exit 1
    fi
    echo -e "  ${YELLOW}CUDA compilation failed${NC}"
    if [[ -f "$CUDA_SO" ]]; then
      echo -e "  ${GREEN}Using pre-built CUDA library${NC}"
      CUDA_STATUS="pre-built"
    else
      echo -e "  ${YELLOW}No CUDA library — will run in CPU mode only${NC}"
      CUDA_STATUS="none"
    fi
  fi
else
  if [[ "$CUDA_MODE" == "required" ]]; then
    echo -e "${RED}Error: nvcc not found but --force-cuda was set${NC}"
    echo "Install CUDA toolkit: https://developer.nvidia.com/cuda-downloads"
    exit 1
  fi
  echo "  nvcc not found — skipping CUDA compilation"
  if [[ -f "$CUDA_SO" ]]; then
    echo -e "  ${GREEN}Pre-built CUDA library found — GPU mode available${NC}"
    CUDA_STATUS="pre-built"
  else
    echo -e "  ${YELLOW}No CUDA library — will run in CPU mode only${NC}"
    CUDA_STATUS="none"
  fi
fi

# ── Done! ──
echo ""
echo -e "${BOLD}${GREEN}"
if [[ "$QUICK" == true ]]; then
  echo "╔══════════════════════════════════════╗"
  echo "║        Build Complete!               ║"
  echo "╚══════════════════════════════════════╝"
else
  echo "╔══════════════════════════════════════╗"
  echo "║     Installation Complete!           ║"
  echo "╚══════════════════════════════════════╝"
fi
echo -e "${NC}"

echo -e "  JAR:   ${GREEN}$JAR_PATH${NC}"
if [[ "$CUDA_STATUS" == "none" ]]; then
  echo -e "  CUDA:  ${YELLOW}not available (CPU mode only)${NC}"
else
  echo -e "  CUDA:  ${GREEN}$CUDA_SO${NC}"
fi

echo ""
echo -e "${BOLD}Quick Start:${NC}"
echo "  ./run.sh -i <gene_trees.tre> -o <output.tre>               # GPU mode (default)"
echo "  ./run.sh -i <gene_trees.tre> -o <output.tre> --cpu-parallel # CPU parallel mode"
echo ""
echo -e "${BOLD}With monitoring:${NC}"
echo "  ./run-with-monitor.sh -i <gene_trees.tre> -o <output.tre>"
