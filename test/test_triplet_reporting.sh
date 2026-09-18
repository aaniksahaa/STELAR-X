#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-triplet-reporting.XXXXXX")"
trap 'status=$?; rm -rf -- "$WORK"; exit "$status"' EXIT

MOCK_BIN="${WORK}/bin"
mkdir -p "$MOCK_BIN"
printf '%s\n' \
  '#!/usr/bin/env bash' \
  'printf "%s\n" "$@" > "$NTFY_CAPTURE"' \
  > "${MOCK_BIN}/curl"
chmod +x "${MOCK_BIN}/curl"

GENES="${ROOT}/test/input/test_5taxa.tre"
CANDIDATE="${ROOT}/test/input/stelar_candidate_5taxa.tre"
OUTPUT="${WORK}/inferred.tre"
MONITOR_CAPTURE="${WORK}/monitor-ntfy.txt"

PATH="${MOCK_BIN}:${PATH}" NTFY_CAPTURE="$MONITOR_CAPTURE" \
  "${ROOT}/run-stelarx-with-monitor.sh" \
  --input "$GENES" --output "$OUTPUT" \
  --opts "--cpu -q --no-build" --no-time-monitor --no-gpu-monitor \
  > "${WORK}/monitor.log" 2>&1

STATS="${OUTPUT%.tre}_stats.csv"
grep -q 'optimal_triplet_score' "$STATS"
! grep -qi 'quartet' "$STATS"
TRIPLET_SCORE="$(awk -F, 'NR==2 {print $7}' "$STATS")"
[[ "$TRIPLET_SCORE" =~ ^[0-9]+([.][0-9]+)?$ ]]
grep -Eq "Triplet score:[[:space:]]*${TRIPLET_SCORE}" "${WORK}/monitor.log"
grep -Eq "Triplet score:[[:space:]]*${TRIPLET_SCORE}" "$MONITOR_CAPTURE"
! grep -qi 'quartet' "${WORK}/monitor.log" "$MONITOR_CAPTURE"

SCORE_CAPTURE="${WORK}/score-only-ntfy.txt"
PATH="${MOCK_BIN}:${PATH}" NTFY_CAPTURE="$SCORE_CAPTURE" NO_COLOR=1 \
  "${ROOT}/run.sh" --no-build --input "$GENES" \
  --score-species-tree "$CANDIDATE" --cpu -q --xms 64m --xmx 1g \
  > "${WORK}/score-only.log" 2>&1

grep -q '^TRIPLET_SCORE: 21$' "${WORK}/score-only.log"
grep -q 'STELAR-X score-only completed' "$SCORE_CAPTURE"
grep -q 'Triplet score: 21' "$SCORE_CAPTURE"
! grep -qi 'quartet' "${WORK}/score-only.log" "$SCORE_CAPTURE"

echo "Triplet score reporting and ntfy: PASS"
