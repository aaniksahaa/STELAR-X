#!/usr/bin/env bash
# STELAR-X command-line and identity smoke tests.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-cli-tests.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

"${ROOT}/build.sh" >/dev/null

[[ -f "${ROOT}/build/stelarx/Main.class" ]]
[[ "$(find "${ROOT}/build" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')" == "stelarx" ]]

VERSION_TEXT="$(NO_COLOR=1 "${ROOT}/stelarx" --version --no-build)"
[[ "$VERSION_TEXT" == *"STELAR-X  v2.0.0"* ]]
[[ "$VERSION_TEXT" == *"Welcome to STELAR-X version 2.0.0!"* ]]

HELP_TEXT="$(NO_COLOR=1 "${ROOT}/stelarx" --help 2>&1)"
[[ "$HELP_TEXT" == *"STELAR-X wrapper"* ]]
[[ "$HELP_TEXT" == *"--search-space"* ]]
[[ "$HELP_TEXT" == *"--intersection-method"* ]]
[[ "$HELP_TEXT" == *"--gpu-strict"* ]]

DIAG_TEXT="$(NO_COLOR=1 java -Djava.library.path="${ROOT}/native" \
  -cp "${ROOT}/build" stelarx.Main --cpu --diagnose 2>&1)"
[[ "$DIAG_TEXT" == *"STELAR-X DIAGNOSTICS"* ]]
[[ "$DIAG_TEXT" == *"Selected compute:            CPU"* ]]

OUTPUT_TREE="${WORK}/species-tree.tre"
STELARX_CRASH_DIR="${WORK}/launcher-crash_logs" \
  NO_COLOR=1 "${ROOT}/stelarx" --no-build --cpu -q \
  -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
  -o "$OUTPUT_TREE" --search-space S1 --intersection-method I2 >/dev/null
[[ -s "$OUTPUT_TREE" ]]
[[ -d "${WORK}/launcher-crash_logs" ]]

STELARX_SKIP_BUILD=1 "${ROOT}/test/run_stelarx_tests.sh" >/dev/null

echo "STELAR-X CLI and identity tests: PASS"
