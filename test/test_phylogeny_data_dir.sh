#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/scripts/phylogeny-data-dir.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-data-dir-test.XXXXXX")"
trap 'rm -rf -- "$WORK"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

DEFAULT_BASE="${WORK}/phylogeny data"
RESOLVED="$(PHYLOGENY_DATA_DIR="${DEFAULT_BASE}/" stelarx_prepare_simphy_data_dir "")"
[[ "$RESOLVED" == "${DEFAULT_BASE}/simphy/data" ]] || fail "unexpected default path: $RESOLVED"
[[ -d "${DEFAULT_BASE}/simphy/data" ]] || fail "default directory was not created"

OVERRIDE="${WORK}/explicit data"
RESOLVED="$(env -u PHYLOGENY_DATA_DIR bash -c \
  'source "$1/scripts/phylogeny-data-dir.sh"; stelarx_prepare_simphy_data_dir "$2"' \
  _ "$ROOT" "$OVERRIDE")"
[[ "$RESOLVED" == "$OVERRIDE" ]] || fail "explicit override was not preserved"
[[ -d "$OVERRIDE" ]] || fail "override directory was not created"

CONFLICT="${WORK}/not-a-directory"
touch "$CONFLICT"
if PHYLOGENY_DATA_DIR="$CONFLICT" stelarx_prepare_simphy_data_dir "$CONFLICT" >"${WORK}/conflict.out" 2>&1; then
  fail "a file was accepted as the SimPhy data directory"
fi
grep -q "not a directory" "${WORK}/conflict.out" || fail "file-conflict error was unclear"

DEFAULT_DATASET="${DEFAULT_BASE}/simphy/data/t_1_g_1_sb_0.000001_spmin_500000_spmax_1500000/R1"
mkdir -p "$DEFAULT_DATASET"
touch "${DEFAULT_DATASET}/stat-sim.csv"
PHYLOGENY_DATA_DIR="$DEFAULT_BASE" "${ROOT}/scripts/sim.sh" -t 1 -g 1 >"${WORK}/default-sim.out"
grep -Fq "SKIPPING: ${DEFAULT_DATASET}/stat-sim.csv" "${WORK}/default-sim.out" || \
  fail "sim.sh did not use the environment-derived checkpoint path"

OVERRIDE_DATASET="${OVERRIDE}/t_2_g_3_sb_0.000001_spmin_500000_spmax_1500000/R1"
mkdir -p "$OVERRIDE_DATASET"
touch "${OVERRIDE_DATASET}/stat-sim.csv"
env -u PHYLOGENY_DATA_DIR "${ROOT}/scripts/sim.sh" -t 2 -g 3 --simphy-data-dir "$OVERRIDE" \
  >"${WORK}/override-sim.out"
grep -Fq "SKIPPING: ${OVERRIDE_DATASET}/stat-sim.csv" "${WORK}/override-sim.out" || \
  fail "sim.sh did not honor the explicit override"

if env -u PHYLOGENY_DATA_DIR "${ROOT}/scripts/sim.sh" -t 4 -g 5 >"${WORK}/missing.out" 2>&1; then
  fail "sim.sh succeeded without an environment default or override"
fi
grep -q "PHYLOGENY_DATA_DIR is not set" "${WORK}/missing.out" || \
  fail "missing-environment error was unclear"

STATS_BASE="${WORK}/stats root"
PHYLOGENY_DATA_DIR="$STATS_BASE" "${ROOT}/scripts/collect-stats-simulated.sh" \
  --out "${WORK}/unused.csv" >"${WORK}/stats.out"
[[ -d "${STATS_BASE}/simphy/data" ]] || fail "stats collector did not create the default directory"
grep -q "No stat files" "${WORK}/stats.out" || fail "empty stats directory was not handled"

DOWNLOAD_BASE="${WORK}/download root"
PHYLOGENY_DATA_DIR="$DOWNLOAD_BASE" "${ROOT}/scripts/download-bulk-simulated.sh" \
  --dry-run --download-script /bin/true --taxa-list 1 --gene-trees-list 1 \
  --sb-list 0.1 --spmin-list 1 --spmax-list 1 >"${WORK}/download.out"
[[ -d "${DOWNLOAD_BASE}/simphy/data" ]] || fail "download tool did not create the default directory"

UPLOAD_BASE="${WORK}/upload root"
PHYLOGENY_DATA_DIR="$UPLOAD_BASE" "${ROOT}/scripts/upload-bulk-simulated.sh" \
  --dry-run --uploader /bin/true --python /bin/true >"${WORK}/upload.out"
[[ -d "${UPLOAD_BASE}/simphy/data" ]] || fail "upload tool did not create the default directory"

OUTPUTS_BASE="${WORK}/outputs root"
PHYLOGENY_DATA_DIR="$OUTPUTS_BASE" "${ROOT}/scripts/upload-bulk-simulated-outputs.sh" \
  --dry-run --uploader /bin/true --python /bin/true >"${WORK}/upload-outputs.out"
[[ -d "${OUTPUTS_BASE}/outputs/simphy" ]] || fail "outputs upload tool did not create the default outputs directory"
grep -q "No <method>_outputs/<dataset> directories matched" "${WORK}/upload-outputs.out" || \
  fail "empty outputs mirror was not handled"
PHYLOGENY_DATA_DIR="$OUTPUTS_BASE" "${ROOT}/scripts/sync-simulated-outputs.sh" --dry-run >"${WORK}/sync-outputs.out"
[[ -d "${OUTPUTS_BASE}/simphy/data" && -d "${OUTPUTS_BASE}/outputs/simphy" ]] || \
  fail "outputs sync tool did not create the default directories"

echo "PASS: PHYLOGENY_DATA_DIR defaults, creation, and overrides"
