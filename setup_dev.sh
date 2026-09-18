#!/usr/bin/env bash
# Prepare and verify a Linux STELAR-X development checkout.
set -euo pipefail

STELARX_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${STELARX_ROOT}/.venv"
CPU_ONLY=false
CHECK_ONLY=false
BUILD_PROJECT=true
RUN_TESTS=true
CUDA_ARCH_VALUE="all-major"

usage() {
  cat <<'EOF'
Usage: ./setup_dev.sh [options]

Options:
  --cpu-only          Skip CUDA compilation
  --cuda-arch VALUE   CUDA target passed to nvcc (default: all-major)
  --no-build          Create/check the Python environment without compiling
  --no-tests          Skip the CPU regression tests
  --check             Only verify the current environment; change nothing
  -h, --help          Show this help

The script creates .venv, installs requirements-dev.txt, builds Java, builds
CUDA libraries when nvcc is available, and runs the CPU tests. System packages
(JDK, Python, CUDA toolkit) are checked but never installed with sudo.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cpu-only) CPU_ONLY=true; shift ;;
    --cuda-arch) [[ $# -ge 2 ]] || { echo "--cuda-arch requires a value" >&2; exit 2; }; CUDA_ARCH_VALUE="$2"; shift 2 ;;
    --no-build) BUILD_PROJECT=false; RUN_TESTS=false; shift ;;
    --no-tests) RUN_TESTS=false; shift ;;
    --check) CHECK_ONLY=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

failures=0
warnings=0

ok()   { printf '  [OK]   %s\n' "$1"; }
warn() { printf '  [WARN] %s\n' "$1"; warnings=$((warnings + 1)); }
fail() { printf '  [FAIL] %s\n' "$1"; failures=$((failures + 1)); }

require_command() {
  local command_name="$1"
  local description="$2"
  if command -v "$command_name" >/dev/null 2>&1; then
    ok "$description: $(command -v "$command_name")"
  else
    fail "$description ('$command_name' not found)"
  fi
}

echo "=== STELAR-X developer environment ==="
echo "Repository: $STELARX_ROOT"
echo
echo "Required tools"
require_command bash "Bash"
require_command java "Java runtime"
require_command javac "Java compiler"
require_command jar "JAR tool"
require_command python3 "Python 3"
require_command realpath "GNU realpath"
require_command awk "awk"
require_command sed "sed"
require_command grep "grep"
require_command find "find"

if command -v javac >/dev/null 2>&1; then
  JAVA_MAJOR="$(javac -version 2>&1 | awk '{print $2}' | cut -d. -f1)"
  if [[ "$JAVA_MAJOR" =~ ^[0-9]+$ ]] && (( JAVA_MAJOR >= 21 )); then
    ok "JDK version: $(javac -version 2>&1)"
  else
    fail "JDK 21 or newer is required (found: $(javac -version 2>&1))"
  fi
fi

if command -v python3 >/dev/null 2>&1; then
  if python3 -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 9) else 1)'; then
    ok "Python version: $(python3 --version 2>&1)"
  else
    fail "Python 3.9 or newer is required"
  fi
fi

echo
echo "Optional integrations"
if [[ -x /usr/bin/time ]]; then ok "GNU time monitor: /usr/bin/time"; else warn "GNU time not found; monitoring will omit peak CPU RAM"; fi
if command -v curl >/dev/null 2>&1; then ok "curl notifications"; else warn "curl not found; notifications will be disabled"; fi
if command -v nvidia-smi >/dev/null 2>&1; then ok "NVIDIA driver tool: $(command -v nvidia-smi)"; else warn "nvidia-smi not found; GPU monitoring will be disabled"; fi
if command -v nvcc >/dev/null 2>&1; then ok "CUDA compiler: $(nvcc --version | tail -1)"; else warn "nvcc not found; only CPU development builds are available"; fi
if [[ -x "${STELARX_ROOT}/simphy/simphy_lnx64" ]]; then ok "Bundled SimPhy executable"; else warn "Bundled SimPhy executable is missing or not executable"; fi

if (( failures > 0 )); then
  echo
  echo "Environment check failed with ${failures} required item(s) missing." >&2
  echo "On Ubuntu/Debian, install the base tools with:" >&2
  echo "  sudo apt install openjdk-21-jdk python3 python3-venv time curl" >&2
  exit 1
fi

PYTHON_BIN="python3"
if [[ -x "${VENV_DIR}/bin/python" ]]; then
  PYTHON_BIN="${VENV_DIR}/bin/python"
fi

if [[ "$CHECK_ONLY" == true ]]; then
  if "$PYTHON_BIN" -c 'import dendropy' >/dev/null 2>&1; then
    ok "Python dependency: DendroPy"
  else
    fail "DendroPy is missing; run ./setup_dev.sh"
  fi
  if [[ -f "${STELARX_ROOT}/build/stelarx/Main.class" ]]; then ok "Java build output"; else warn "Java build output is absent; run ./setup_dev.sh"; fi
  if [[ "$CPU_ONLY" != true && -x "$(command -v nvcc 2>/dev/null || true)" ]]; then
    if [[ -f "${STELARX_ROOT}/native/libstelarx_weight.so" ]]; then ok "CUDA native libraries"; else warn "CUDA libraries are absent; run ./setup_dev.sh"; fi
  fi
  echo
  (( failures == 0 )) || exit 1
  echo "Developer environment is ready (${warnings} optional warning(s))."
  exit 0
fi

echo
if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
  echo "Creating Python environment: ${VENV_DIR}"
  if ! python3 -m venv "$VENV_DIR"; then
    echo "Failed to create .venv. Install the python3-venv system package." >&2
    exit 1
  fi
fi
PYTHON_BIN="${VENV_DIR}/bin/python"
"$PYTHON_BIN" -m pip install --disable-pip-version-check -r "${STELARX_ROOT}/requirements-dev.txt"

if [[ "$BUILD_PROJECT" == true ]]; then
  "${STELARX_ROOT}/build.sh"
  if [[ "$CPU_ONLY" != true ]]; then
    if command -v nvcc >/dev/null 2>&1; then
      CUDA_ARCH="$CUDA_ARCH_VALUE" "${STELARX_ROOT}/build_native.sh"
    else
      warn "Skipping CUDA build because nvcc is unavailable"
    fi
  fi
fi

if [[ "$RUN_TESTS" == true ]]; then
  bash "${STELARX_ROOT}/test/run_stelarx_tests.sh"
fi

echo
echo "Developer setup complete."
echo "  Check later : ./setup_dev.sh --check"
echo "  Run         : ./stelarx -i rooted_gene_trees.tre -o species_tree.tre --search-space S1"
echo "  Monitor     : ./run-stelarx-with-monitor.sh -i gene_trees.tre -o species_tree.tre --search-space S1"
