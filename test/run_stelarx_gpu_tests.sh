#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-gpu-tests.XXXXXX")"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

if [[ "${STELARX_SKIP_BUILD:-0}" != 1 ]]; then
  "${ROOT}/build.sh" >/dev/null
fi
if [[ ! -f "${ROOT}/native/libstelarx_weight.so" ]]; then
  "${ROOT}/build_native.sh" >/dev/null
fi

JAVA=(java -Djava.library.path="${ROOT}/native" -cp "${ROOT}/build" stelarx.Main)

"${JAVA[@]}" --gpu-strict --diagnose >"${WORK}/diagnose.log" 2>&1
grep -q "CUDA.*" "${WORK}/diagnose.log"
grep -q "Usable:.*yes" "${WORK}/diagnose.log"

score_gpu() {
  local genes="$1" species="$2" method="$3" expected="$4" log="$5"
  shift 5
  "${JAVA[@]}" --gpu-strict -q -i "$genes" \
    --score-species-tree "$species" --intersection-method "$method" \
    "$@" \
    >"$log" 2>&1
  grep -q "\[STELAR-X GPU\] weight" "$log"
  [[ "$(sed -n 's/^TRIPLET_SCORE: //p' "$log" | tail -1)" == "$expected" ]]
}

for method in I1 I2 I3 I4; do
  score_gpu "${ROOT}/test/input/test_5taxa.tre" \
    "${ROOT}/test/input/stelar_candidate_5taxa.tre" "$method" 21 \
    "${WORK}/small-${method}.log"
  score_gpu "${ROOT}/test/input/stelar_polytomy_5taxa.tre" \
    "${ROOT}/test/input/stelar_candidate_5taxa.tre" "$method" 11 \
    "${WORK}/polytomy-${method}.log"
  score_gpu "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
    "${ROOT}/test/input/stelar_candidate_6taxa.tre" "$method" 22 \
    "${WORK}/incomplete-polytomy-${method}.log"
  score_gpu "${ROOT}/test/input/child_only_dedup_incomplete.tre" \
    "${ROOT}/test/input/child_only_candidate_7taxa.tre" "$method" 40 \
    "${WORK}/child-only-${method}.log"
  score_gpu "${ROOT}/all_gt_bs_rooted_37.tre" "${ROOT}/true_37.tre" \
    "$method" 1390544 "${WORK}/large-${method}.log"
done

# The optimization-specific incomplete-tree fixture must also stay exact in both
# wide accumulator kernels; forcing the boundary selector makes these paths run
# even though the fixture itself is small.
for numeric in int128 double; do
  for method in I1 I2 I3 I4; do
    STELARX_WEIGHT_FORCE_DOUBLE=1 score_gpu \
      "${ROOT}/test/input/child_only_dedup_incomplete.tre" \
      "${ROOT}/test/input/child_only_candidate_7taxa.tre" "$method" 40 \
      "${WORK}/child-only-${method}-${numeric}.log" \
      --large-n-score-type "$numeric"
  done
done

# Every method must retain exact results under the three batching controls.
for method in I1 I2 I3 I4; do
  score_gpu "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
    "${ROOT}/test/input/stelar_candidate_6taxa.tre" "$method" 22 \
    "${WORK}/batch-size-${method}.log" --gpu-batch-size 1
  score_gpu "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
    "${ROOT}/test/input/stelar_candidate_6taxa.tre" "$method" 22 \
    "${WORK}/batch-count-${method}.log" --gpu-batches 2
  score_gpu "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
    "${ROOT}/test/input/stelar_candidate_6taxa.tre" "$method" 22 \
    "${WORK}/batch-off-${method}.log" --no-gpu-batch
done

for preset in S1 S2 S3; do
  "${JAVA[@]}" --gpu-strict -q -i "${ROOT}/test/input/test_incomplete.tre" \
    --search-space "$preset" --intersection-method I2 \
    -o "${WORK}/${preset}.tre" >"${WORK}/${preset}.log" 2>&1
  test -s "${WORK}/${preset}.tre"
  grep -q "Phase 6  Weight calculation.*\[GPU\]" "${WORK}/${preset}.log"
  grep -q "Triplet score.*122" "${WORK}/${preset}.log"
  if [[ "$preset" != S1 ]]; then
    grep -q "Phase 1b Auto-complete.*\[GPU\]" "${WORK}/${preset}.log"
    grep -q "Phase 5b Cross-tree transitions.*\[GPU\]" "${WORK}/${preset}.log"
  fi
done

# End-to-end native-polytomy inference: CPU and strict GPU must select the same
# topology and objective score while S2 also exercises completion/similarity/DP.
"${JAVA[@]}" --cpu -q -i "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
  --search-space S2 --intersection-method I2 --keep-polytomy-during-inference \
  -o "${WORK}/poly-cpu.tre" >"${WORK}/poly-cpu.log" 2>&1
"${JAVA[@]}" --gpu-strict -q \
  -i "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
  --search-space S2 --intersection-method I2 --keep-polytomy-during-inference \
  -o "${WORK}/poly-gpu.tre" >"${WORK}/poly-gpu.log" 2>&1
cmp "${WORK}/poly-cpu.tre" "${WORK}/poly-gpu.tre"
grep -q "Phase 6  Weight calculation.*\[GPU\]" "${WORK}/poly-gpu.log"
grep -q "Triplet score.*27" "${WORK}/poly-gpu.log"

# Seeded independent oracles exercise many additional binary, incomplete, and
# internally polytomous layouts, including forced DOUBLE and INT128 kernels.
python3 "${ROOT}/test/test_stelarx_differential.py" \
  --gpu --cases "${STELARX_GPU_RANDOM_CASES:-6}" --no-build
python3 "${ROOT}/test/test_stelarx_inference.py" --gpu --no-build

echo "STELAR-X strict CUDA suite: PASS"
