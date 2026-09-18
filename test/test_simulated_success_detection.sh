#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-success-test.XXXXXX")"
trap 'status=$?; rm -rf -- "$TMP"; exit "$status"' EXIT

DATA="$TMP/data"
RUN_DIR="$DATA/t_4_g_1_sb_0.000001_spmin_100000_spmax_200000/R1"
mkdir -p "$RUN_DIR"
printf '((a,b),(c,d));\n' > "$RUN_DIR/all_gt.tre"
printf '((a,b),(c,d));\n' > "$RUN_DIR/s_tree.trees"

COMMON=(--simphy-data-dir "$DATA" -t 4 -g 1 -r R1
  --sb 0.000001 --spmin 100000 --spmax 200000
  --opts '--search-space S1 --cpu -q'
  --no-time-monitor --no-gpu-monitor --no-notify)

"$ROOT/scripts/test-stelarx-simulated.sh" "${COMMON[@]}" >/dev/null
RESULTS_DIR=$(find "$RUN_DIR/stelarx_outputs" -mindepth 1 -maxdepth 1 -type d)
OUTPUT="$RESULTS_DIR/out-stelarx.tre"
SIDE="$RESULTS_DIR/out-stelarx_stats.csv"
SUCCESS="$RESULTS_DIR/.stelarx.success"
[[ -s "$OUTPUT" && -s "$SIDE" && -s "$SUCCESS" ]]
STAT_FILE="$RESULTS_DIR/stat-stelarx.csv"
grep -q 'optimal-triplet-score' "$STAT_FILE"
! grep -qi 'quartet' "$STAT_FILE"
[[ "$(awk -F, 'NR==2 {print $10}' "$STAT_FILE")" =~ ^[0-9]+([.][0-9]+)?$ ]]

# The reproducibility mirror sits beside the data tree and equals the results leaf.
MIRROR_LEAF="$TMP/outputs/stelarx_outputs/t_4_g_1_sb_0.000001_spmin_100000_spmax_200000/R1/$(basename "$RESULTS_DIR")"
[[ -s "$MIRROR_LEAF/out-stelarx.tre" && -s "$MIRROR_LEAF/stat-stelarx.csv" ]]
diff -r "$RESULTS_DIR" "$MIRROR_LEAF" >/dev/null
[[ ! -e "$TMP/outputs/stelarx_outputs/t_4_g_1_sb_0.000001_spmin_100000_spmax_200000/R1/all_gt.tre" ]]

# A failed sidecar plus a stale tree must never be accepted as completed.
rm -f "$SUCCESS"
sed -i '2s/,0$/,1/' "$SIDE"
rerun_log=$("$ROOT/scripts/test-stelarx-simulated.sh" "${COMMON[@]}" 2>&1)
[[ "$rerun_log" == *"Previous statistics exist but no successful output was recorded; rerunning."* ]]
[[ -s "$OUTPUT" && -s "$SUCCESS" ]]
[[ "$(awk -F, 'NR==2 {print $9}' "$SIDE")" == "0" ]]

skip_log=$("$ROOT/scripts/test-stelarx-simulated.sh" "${COMMON[@]}" 2>&1)
[[ "$skip_log" == *"SKIPPING: successful output already exists"* ]]

COMBINED="${TMP}/combined.csv"
"${ROOT}/scripts/collect-stats-simulated.sh" --simphy-data-dir "$DATA" --out "$COMBINED" >/dev/null
grep -q 'optimal-triplet-score' "$COMBINED"
! grep -qi 'quartet' "$COMBINED"
[[ "$(awk -F, 'NR==2 {print $10}' "$COMBINED")" =~ ^[0-9]+([.][0-9]+)?$ ]]

echo "Simulated-run success detection: PASS"
