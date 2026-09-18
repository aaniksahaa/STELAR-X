#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-clear-simulated-test.XXXXXX")"
trap 'rm -rf -- "$WORK"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

CLEAR_SCRIPT="${ROOT}/clear-bulk-simulated.sh"

"$CLEAR_SCRIPT" --help >"${WORK}/help.out"
grep -Fq '$PHYLOGENY_DATA_DIR/simphy/data' "${WORK}/help.out" || \
  fail "help does not state the exact deletion target"

if env -u PHYLOGENY_DATA_DIR "$CLEAR_SCRIPT" --yes >"${WORK}/unset.out" 2>&1; then
  fail "cleanup succeeded without PHYLOGENY_DATA_DIR"
fi
grep -q "PHYLOGENY_DATA_DIR is not set" "${WORK}/unset.out" || \
  fail "missing-environment error was unclear"

if PHYLOGENY_DATA_DIR="relative/path" "$CLEAR_SCRIPT" --yes \
    >"${WORK}/relative.out" 2>&1; then
  fail "cleanup accepted a relative PHYLOGENY_DATA_DIR"
fi
grep -q "must be an absolute path" "${WORK}/relative.out" || \
  fail "relative-path error was unclear"

BASE="${WORK}/phylogeny data"
TARGET="${BASE}/simphy/data"
mkdir -p "${TARGET}/dataset/results" "${BASE}/simphy/keep"
touch "${TARGET}/dataset/all_gt.tre"
touch "${TARGET}/dataset/results/out-stelarx.tre"
touch "${BASE}/simphy/keep/sentinel"

PHYLOGENY_DATA_DIR="$BASE" "$CLEAR_SCRIPT" --dry-run >"${WORK}/dry-run.out"
grep -Fq "Delete completely:   $TARGET" "${WORK}/dry-run.out" || \
  fail "dry run did not print the exact target"
[[ -f "${TARGET}/dataset/results/out-stelarx.tre" ]] || \
  fail "dry run removed data"

if PHYLOGENY_DATA_DIR="$BASE" "$CLEAR_SCRIPT" >"${WORK}/noninteractive.out" 2>&1; then
  fail "non-interactive cleanup succeeded without --yes"
fi
grep -q "without --yes" "${WORK}/noninteractive.out" || \
  fail "non-interactive refusal was unclear"
[[ -d "$TARGET" ]] || fail "refused cleanup still removed the target"

PHYLOGENY_DATA_DIR="$BASE" "$CLEAR_SCRIPT" --yes >"${WORK}/delete.out"
[[ ! -e "$TARGET" ]] || fail "--yes did not remove the complete data directory"
[[ -f "${BASE}/simphy/keep/sentinel" ]] || \
  fail "cleanup removed a sibling outside the data directory"

PHYLOGENY_DATA_DIR="$BASE" "$CLEAR_SCRIPT" --yes >"${WORK}/absent.out"
grep -q "nothing to remove" "${WORK}/absent.out" || \
  fail "an absent data directory was not handled idempotently"

OUTSIDE="${WORK}/outside"
mkdir -p "$OUTSIDE" "${BASE}/simphy"
touch "${OUTSIDE}/sentinel"
ln -s "$OUTSIDE" "$TARGET"
if PHYLOGENY_DATA_DIR="$BASE" "$CLEAR_SCRIPT" --yes >"${WORK}/symlink.out" 2>&1; then
  fail "cleanup accepted a symbolic-link target"
fi
grep -q "symbolic link" "${WORK}/symlink.out" || \
  fail "symbolic-link refusal was unclear"
[[ -f "${OUTSIDE}/sentinel" ]] || fail "cleanup followed the symbolic link"

echo "PASS: clear-bulk-simulated path safety and complete deletion"

