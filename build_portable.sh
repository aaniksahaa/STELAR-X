#!/usr/bin/env bash
# Build a self-contained STELAR-X application image and transport archive.
# Target machines need neither Java nor the CUDA toolkit. NVIDIA acceleration
# is optional; an unusable/missing driver falls back to CPU at runtime.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
DIST_DIR="${ROOT}/dist"
CUDA_MODE="auto"          # auto | off | required
CUDA_ARCH_VALUE="all-major"
MAKE_ARCHIVE=true
VERSION=""
FORCE=false

usage() {
  cat <<'EOF'
Usage: ./build_portable.sh [options]

Options:
  --without-cuda          Build without CUDA libraries
  --cpu-only              Alias for --without-cuda
  --with-cuda             Require CUDA instead of using the automatic default
  --cuda-arch VALUE       nvcc architecture (default: all-major)
  --version VERSION       Release version (default: source version)
  --output-dir DIR        Artifact destination (default: ./dist)
  --no-archive            Keep only the unpacked application image
  --force                 Replace this platform's existing artifact/version
  -h, --help              Show this help

The produced application is native to the build host's OS and CPU architecture.
CUDA is bundled automatically on Linux when nvcc is available; every CUDA build
also contains the complete CPU fallback.
Run this script on each target platform (normally through the release CI matrix).
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --without-cuda|--cpu-only) CUDA_MODE="off"; shift ;;
    --with-cuda) CUDA_MODE="required"; shift ;;
    --cuda-arch) [[ $# -ge 2 ]] || { echo "--cuda-arch requires a value" >&2; exit 2; }; CUDA_ARCH_VALUE="$2"; shift 2 ;;
    --version) [[ $# -ge 2 ]] || { echo "--version requires a value" >&2; exit 2; }; VERSION="$2"; shift 2 ;;
    --output-dir) [[ $# -ge 2 ]] || { echo "--output-dir requires a value" >&2; exit 2; }; DIST_DIR="$2"; shift 2 ;;
    --no-archive) MAKE_ARCHIVE=false; shift ;;
    --force) FORCE=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for tool in javac jar jlink jpackage; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Error: '$tool' is required to build the artifact (JDK 21 expected)." >&2
    exit 1
  fi
done

JAVA_MAJOR="$(javac -version 2>&1 | awk '{print $2}' | cut -d. -f1)"
if [[ "$JAVA_MAJOR" -lt 21 ]]; then
  echo "Error: JDK 21 or newer is required; found javac $(javac -version 2>&1)." >&2
  exit 1
fi

OS_RAW="$(uname -s)"
ARCH_RAW="$(uname -m)"
case "$OS_RAW" in
  Linux)  PLATFORM_OS="linux" ;;
  Darwin) PLATFORM_OS="macos" ;;
  *) echo "Error: $OS_RAW is not supported by this script; use build_portable.ps1 on Windows." >&2; exit 1 ;;
esac
case "$ARCH_RAW" in
  x86_64|amd64) PLATFORM_ARCH="x86_64" ;;
  arm64|aarch64) PLATFORM_ARCH="arm64" ;;
  *) PLATFORM_ARCH="$ARCH_RAW" ;;
esac

SOURCE_VERSION="$(sed -n 's/.*DEFAULT = "\([^"]*\)".*/\1/p' "${ROOT}/src/stelarx/Version.java" | head -1)"
if [[ -z "$SOURCE_VERSION" ]]; then
  echo "Error: could not determine the STELAR-X source version." >&2
  exit 1
fi
[[ -n "$VERSION" ]] || VERSION="$SOURCE_VERSION"
if [[ ! "$VERSION" =~ ^[0-9]+([.][0-9]+){0,2}$ ]]; then
  echo "Error: --version must contain one to three numeric components (for example 1.2.0)." >&2
  exit 2
fi

INCLUDE_CUDA=false
if [[ "$CUDA_MODE" != "off" && "$PLATFORM_OS" == "linux" ]] && command -v nvcc >/dev/null 2>&1; then
  INCLUDE_CUDA=true
elif [[ "$CUDA_MODE" == "required" ]]; then
  echo "Error: --with-cuda requires Linux and an nvcc installation." >&2
  exit 1
fi

CAPABILITY="cpu"
[[ "$INCLUDE_CUDA" == true ]] && CAPABILITY="cuda-with-cpu-fallback"
# There is exactly one public artifact name per OS/CPU family.  On CUDA-capable
# platforms the richer build still contains the complete CPU implementation and
# automatically falls back, so separate -cpu and -cuda downloads are needless.
ARTIFACT="stelarx-${VERSION}-${PLATFORM_OS}-${PLATFORM_ARCH}"
VERSION_DIR="${DIST_DIR}/${VERSION}"
FINAL_IMAGE="${VERSION_DIR}/${ARTIFACT}"
ARCHIVE="${VERSION_DIR}/${ARTIFACT}.tar.gz"
MANIFEST="${VERSION_DIR}/${ARTIFACT}.manifest.json"
if [[ "$FORCE" != true ]] && { [[ -e "$FINAL_IMAGE" ]] || [[ -e "$ARCHIVE" ]] || [[ -e "$MANIFEST" ]]; }; then
  echo "Error: artifact already exists for STELAR-X ${VERSION} on ${PLATFORM_OS}-${PLATFORM_ARCH}." >&2
  echo "Use --force to replace only this platform artifact, or choose another --version." >&2
  exit 1
fi
mkdir -p "$VERSION_DIR"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-portable.XXXXXX")"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

echo "=== STELAR-X portable build ==="
echo "  Version      : $VERSION"
echo "  Platform     : ${PLATFORM_OS}-${PLATFORM_ARCH}"
echo "  CUDA bundle  : $INCLUDE_CUDA"
echo "  Artifact     : $ARTIFACT"

"${ROOT}/build.sh"

APP_INPUT="${WORK}/input"
RUNTIME="${WORK}/runtime"
JPACKAGE_OUT="${WORK}/jpackage"
mkdir -p "$APP_INPUT" "$JPACKAGE_OUT"

JAR_MANIFEST="${WORK}/MANIFEST.MF"
printf 'Manifest-Version: 1.0\nMain-Class: stelarx.Main\nImplementation-Title: STELAR-X\nImplementation-Version: %s\n\n' \
    "$VERSION" > "$JAR_MANIFEST"
jar --create --file "${APP_INPUT}/stelarx.jar" \
    --manifest "$JAR_MANIFEST" -C "${ROOT}/build" .

if [[ "$INCLUDE_CUDA" == true ]]; then
  NATIVE_OUT_DIR="$APP_INPUT" CUDA_ARCH="$CUDA_ARCH_VALUE" "${ROOT}/build_native.sh"
fi

# The code depends only on java.base. A trimmed runtime removes the target
# machine's Java-installation requirement while remaining a normal HotSpot JVM.
jlink --add-modules java.base \
      --strip-debug --no-header-files --no-man-pages --compress=zip-6 \
      --output "$RUNTIME"

jpackage --type app-image \
    --name stelarx \
    --app-version "$VERSION" \
    --input "$APP_INPUT" \
    --main-jar stelarx.jar \
    --main-class stelarx.Main \
    --runtime-image "$RUNTIME" \
    --dest "$JPACKAGE_OUT" \
    --java-options '-Djava.library.path=$APPDIR' \
    --java-options '-Dfile.encoding=UTF-8' \
    --java-options '-XX:InitialRAMPercentage=2.0' \
    --java-options '-XX:MaxRAMPercentage=85.0' \
    --java-options '-XX:ErrorFile=crash_logs/stelarx-hotspot-crash-%p.log'

if [[ "$PLATFORM_OS" == "macos" ]]; then
  # jpackage produces a standard .app bundle on macOS.  Put it in a small
  # command-line-friendly distribution directory and expose one obvious entry
  # point next to it.  The symlink is preserved by the .tar.gz archive.
  RAW_IMAGE="${JPACKAGE_OUT}/stelarx.app"
  IMAGE="${WORK}/${ARTIFACT}"
  mkdir -p "$IMAGE"
  mv "$RAW_IMAGE" "${IMAGE}/stelarx.app"
  install -m 0755 "${ROOT}/packaging/stelarx-launcher" "${IMAGE}/stelarx"
  PACKAGED_LAUNCHER="${IMAGE}/stelarx"
else
  IMAGE="${JPACKAGE_OUT}/stelarx"
  install -m 0755 "${ROOT}/packaging/stelarx-launcher" "${IMAGE}/stelarx"
  PACKAGED_LAUNCHER="${IMAGE}/stelarx"
fi

# Linux portability is bounded by the newest glibc symbol used by any bundled
# ELF object. Record that auditable floor instead of making an unqualified
# "runs on every Linux" claim. Release runners should build on the oldest
# supported distribution.
MINIMUM_GLIBC=""
if [[ "$PLATFORM_OS" == "linux" ]] && command -v objdump >/dev/null 2>&1; then
  MINIMUM_GLIBC="$(
    find "$IMAGE" -type f -exec objdump -T {} + 2>/dev/null \
      | sed -n 's/.*GLIBC_\([0-9][0-9.]*\).*/\1/p' \
      | sort -V | tail -1 || true
  )"
fi

EXAMPLE_DIR="${IMAGE}/example"
mkdir -p "$EXAMPLE_DIR"
cp "${ROOT}/all_gt_bs_rooted_37.tre" "${EXAMPLE_DIR}/all_gt_37.tre"
cp "${ROOT}/true_37.tre" "${EXAMPLE_DIR}/true_37.tre"

cat > "${IMAGE}/README.txt" <<EOF
STELAR-X ${VERSION} — self-contained ${PLATFORM_OS}-${PLATFORM_ARCH} build

Run:
  ./stelarx --help
  ./stelarx --diagnose
  ./stelarx -i /path/to/rooted_gene_trees.tre -o /path/to/output_species_tree.tre

Ready-made 37-taxon example (run from this directory):
  ./stelarx -i example/all_gt_37.tre -o example/predicted_st_37.tre --search-space S1 -vv

The example directory initially contains:
  all_gt_37.tre   input gene trees
  true_37.tre     reference/true species tree

After the command finishes, example/predicted_st_37.tre contains the inferred
species tree. The reference tree is provided for comparison and is not used by
STELAR-X during inference.

No Java installation is required. This artifact includes CUDA acceleration: ${INCLUDE_CUDA}.
CUDA still requires a compatible NVIDIA GPU and installed NVIDIA driver. If CUDA cannot
be used, STELAR-X automatically explains why and falls back to CPU. Use --gpu-strict only
when falling back would be undesirable.

Unexpected failure reports are stored in crash_logs/ under the directory from
which you launch STELAR-X. Set STELARX_CRASH_DIR to override Java report storage.
EOF

{
  echo "version=${VERSION}"
  echo "platform=${PLATFORM_OS}-${PLATFORM_ARCH}"
  echo "capability=${CAPABILITY}"
  [[ -z "$MINIMUM_GLIBC" ]] || echo "minimum_glibc=${MINIMUM_GLIBC}"
  echo "built_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "java=$(java -version 2>&1 | head -1)"
  if [[ "$INCLUDE_CUDA" == true ]]; then
    echo "nvcc=$(nvcc --version | tail -1)"
    echo "cuda_arch=${CUDA_ARCH_VALUE}"
  fi
  if command -v git >/dev/null 2>&1; then
    echo "git_commit=$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
  fi
} > "${IMAGE}/BUILD-INFO.txt"

# Smoke tests use only the packaged launcher/runtime.
VERSION_OUTPUT="$(NO_COLOR=1 "$PACKAGED_LAUNCHER" --version)"
[[ "$VERSION_OUTPUT" == *"STELAR-X  v${VERSION}"* \
   && "$VERSION_OUTPUT" == *"Welcome to STELAR-X version ${VERSION}!"* ]] || {
  echo "Error: packaged version mismatch: ${VERSION_OUTPUT}" >&2
  exit 1
}
"$PACKAGED_LAUNCHER" --cpu --diagnose >/dev/null
"$PACKAGED_LAUNCHER" --cpu --search-space S2 -q \
    -i "${EXAMPLE_DIR}/all_gt_37.tre" \
    -o "${WORK}/smoke-species-tree.tre" \
    --log-file "${WORK}/smoke-run.log"
if [[ ! -s "${WORK}/smoke-species-tree.tre" ]]; then
  echo "Error: packaged end-to-end inference smoke test produced no tree." >&2
  exit 1
fi
if [[ ! -s "${WORK}/smoke-run.log" ]] || ! grep -q "Run Summary" "${WORK}/smoke-run.log"; then
  echo "Error: packaged --log-file smoke test did not capture the complete run." >&2
  exit 1
fi
if LC_ALL=C grep -q $'\r' "${WORK}/smoke-run.log"; then
  echo "Error: packaged --log-file retained transient progress-bar repaints." >&2
  exit 1
fi

if [[ -e "$FINAL_IMAGE" ]]; then rm -rf "$FINAL_IMAGE"; fi
mv "$IMAGE" "$FINAL_IMAGE"

ARCHIVE_NAME=""
ARCHIVE_SHA256=""
if [[ "$MAKE_ARCHIVE" == true ]]; then
  rm -f "$ARCHIVE" "${ARCHIVE}.sha256"
  tar -C "$VERSION_DIR" -czf "$ARCHIVE" "$ARTIFACT"
  if command -v sha256sum >/dev/null 2>&1; then
    (cd "$VERSION_DIR" && sha256sum "$(basename "$ARCHIVE")" > "$(basename "$ARCHIVE").sha256")
  else
    (cd "$VERSION_DIR" && shasum -a 256 "$(basename "$ARCHIVE")" > "$(basename "$ARCHIVE").sha256")
  fi
  ARCHIVE_NAME="$(basename "$ARCHIVE")"
  ARCHIVE_SHA256="$(awk '{print $1}' "${ARCHIVE}.sha256")"
  echo "  Archive      : $ARCHIVE"
  echo "  SHA-256      : ${ARCHIVE}.sha256"
fi

cat > "$MANIFEST" <<EOF
{
  "version": "${VERSION}",
  "platform": "${PLATFORM_OS}-${PLATFORM_ARCH}",
  "capability": "${CAPABILITY}",
  "minimum_glibc": "${MINIMUM_GLIBC}",
  "archive": "${ARCHIVE_NAME}",
  "sha256": "${ARCHIVE_SHA256}"
}
EOF

echo "Portable application ready: ${FINAL_IMAGE}/stelarx"
echo "Release manifest: $MANIFEST"
