#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-tests.XXXXXX")"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

if [[ "${STELARX_SKIP_BUILD:-0}" != 1 ]]; then
  "${ROOT}/build.sh" >/dev/null
fi
python3 "${ROOT}/test/test_stelarx_triplets.py"
"${ROOT}/test/test_triplet_reporting.sh"

for verifier in --verify-clusters --verify-partitions --verify-dp; do
  java -cp "${ROOT}/build" stelarx.Main --cpu -q \
    -i "${ROOT}/test/input/test_5taxa.tre" "$verifier" \
    >"${WORK}/${verifier#--}.log" 2>&1
  grep -q "ALL ASSERTIONS PASSED" "${WORK}/${verifier#--}.log"
done

for preset in S1 S2 S3; do
  java -cp "${ROOT}/build" stelarx.Main --cpu -q \
    -i "${ROOT}/test/input/test_incomplete.tre" --search-space "$preset" \
    -o "${WORK}/${preset}.tre" >"${WORK}/${preset}.log" 2>&1
  test -s "${WORK}/${preset}.tre"
  grep -q "Triplet score" "${WORK}/${preset}.log"
done

NO_COLOR=1 "${ROOT}/stelarx" --no-build --version >"${WORK}/version.log" 2>&1
grep -q "STELAR-X  v" "${WORK}/version.log"

java -Djava.library.path="${ROOT}/native" -cp "${ROOT}/build" \
  stelarx.Main --cpu --diagnose >"${WORK}/diagnose.log" 2>&1
grep -q "weight / CUDA probe:.*loaded" "${WORK}/diagnose.log"

"${ROOT}/test/test_phylogeny_data_dir.sh"
"${ROOT}/test/test_clear_bulk_simulated.sh"
"${ROOT}/test/test_bulk_simulated_already_completed.sh"
"${ROOT}/test/test_simulated_outputs_mirror.sh"
"${ROOT}/test/test_a10k_outputs_mirror.sh"

echo "STELAR-X focused suite: PASS"
