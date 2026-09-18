#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-comprehensive.XXXXXX")"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

GPU_MODE=auto
QUICK=0
PACKAGING=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --require-gpu) GPU_MODE=require; shift ;;
    --cpu-only) GPU_MODE=off; shift ;;
    --quick) QUICK=1; PACKAGING=0; shift ;;
    --skip-packaging) PACKAGING=0; shift ;;
    -h|--help)
      echo "Usage: $0 [--require-gpu|--cpu-only] [--quick] [--skip-packaging]"
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

started=$SECONDS
echo "=== STELAR-X comprehensive validation ==="
echo "Building Java..."
"${ROOT}/build.sh" >/dev/null

TEST_CLASSES="${WORK}/classes"
mkdir -p "$TEST_CLASSES"
javac -cp "${ROOT}/build" -d "$TEST_CLASSES" \
  "${ROOT}/test/CliPresetsTest.java" \
  "${ROOT}/test/PackedPreflightTest.java" \
  "${ROOT}/test/PackedSimilarityParityTest.java" \
  "${ROOT}/test/SimilarityArgminTest.java" \
  "${ROOT}/test/ThreadingFailureTest.java" \
  "${ROOT}/test/WideSimilarityBoundaryTest.java" \
  "${ROOT}/test/stelarx/FatalReporterTest.java" \
  "${ROOT}/test/stelarx/cluster/ResidualLookupTest.java" \
  "${ROOT}/test/stelarx/completion/PackedMatrixBoundaryTest.java" \
  "${ROOT}/test/stelarx/completion/RootedPolytomyLifecycleTest.java" \
  "${ROOT}/test/stelarx/partition/ChildOnlyDeduplicationTest.java" \
  "${ROOT}/test/stelarx/tree/QuotedTaxonNormalizationTest.java" \
  "${ROOT}/test/stelarx/util/Int128Test.java" \
  "${ROOT}/test/stelarx/weight/WeightModeBoundaryTest.java"

CP="${ROOT}/build:${TEST_CLASSES}"
echo "Running arithmetic, dispatch, threading, and matrix unit tests..."
java -cp "$CP" stelarx.CliPresetsTest
java -cp "$CP" stelarx.util.Int128Test
java -cp "$CP" stelarx.weight.WeightModeBoundaryTest
STELARX_WEIGHT_FORCE_DOUBLE=1 java -cp "$CP" stelarx.weight.WeightModeBoundaryTest double
STELARX_WEIGHT_FORCE_LONG=1 java -cp "$CP" stelarx.weight.WeightModeBoundaryTest long
java -cp "$CP" stelarx.completion.PackedMatrixBoundaryTest
java -cp "$CP" stelarx.FatalReporterTest "${WORK}/fatal-reporter"
java -cp "$CP" stelarx.cluster.ResidualLookupTest \
  "${ROOT}/test/input/test_incomplete.tre"
java -Xmx1g -cp "$CP" PackedPreflightTest
java -cp "$CP" ThreadingFailureTest
if [[ $QUICK -eq 0 ]]; then
  java -Xmx4g -cp "$CP" WideSimilarityBoundaryTest
fi

echo "Running rooted parser/polytomy/completion lifecycle tests..."
java -cp "$CP" stelarx.completion.RootedPolytomyLifecycleTest \
  "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre"
java -cp "$CP" stelarx.partition.ChildOnlyDeduplicationTest \
  "${ROOT}/test/input/child_only_dedup_incomplete.tre"
java -cp "$CP" stelarx.tree.QuotedTaxonNormalizationTest \
  "${ROOT}/test/input/quoted_taxa_equivalent.tre" \
  "${ROOT}/test/input/quoted_taxa_escaped.tre"
python3 "${ROOT}/test/test_taxon_label_normalization.py"
java -cp "$CP" PackedSimilarityParityTest "${ROOT}/test/input/test_incomplete.tre"
java -cp "$CP" SimilarityArgminTest \
  "${ROOT}/test/input/test_incomplete.tre" \
  "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre"

if [[ $QUICK -eq 1 ]]; then
  RANDOM_CASES=3
  MATRIX_SEEDS=(1 2)
else
  RANDOM_CASES=12
  MATRIX_SEEDS=(1 2 3 4 5 6)
fi

echo "Running independent randomized CPU oracles..."
python3 "${ROOT}/test/test_stelarx_differential.py" \
  --cases "$RANDOM_CASES" --no-build
python3 "${ROOT}/test/test_stelarx_inference.py" --no-build
python3 "${ROOT}/test/test_similarity_matrix.py" --stelarx-root "$ROOT" \
  --mode cpu --seeds "${MATRIX_SEEDS[@]}"
python3 "${ROOT}/test/test_distance_matrix.py" --stelarx-root "$ROOT" \
  --mode cpu --seeds "${MATRIX_SEEDS[@]}"
python3 "${ROOT}/test/test_upgma.py" --stelarx-root "$ROOT" \
  --seeds "${MATRIX_SEEDS[@]}"

echo "Running phase verifiers and S1/S2/S3 smoke tests..."
STELARX_SKIP_BUILD=1 "${ROOT}/test/run_stelarx_tests.sh"

JAVA=(java -Dstelarx.crashDir="${WORK}/expected-failure-crash_logs" \
  -Djava.library.path="${ROOT}/native" -cp "${ROOT}/build" stelarx.Main)
expect_failure() {
  local label="$1"
  shift
  if "$@" >"${WORK}/reject-${label}.log" 2>&1; then
    echo "Expected failure was accepted: ${label}" >&2
    exit 1
  fi
}

echo "Running CLI/parser rejection and decorated-Newick tests..."
printf '(A,B,C,D);\n' >"${WORK}/root-polytomy.tre"
printf '(A);\n' >"${WORK}/root-unary.tre"
printf '((A,A),(B,C));\n' >"${WORK}/duplicate.tre"
printf '((A,B),(C,D);\n' >"${WORK}/malformed.tre"
: >"${WORK}/empty.tre"
printf '((A,B),(C,D));\n' >"${WORK}/valid-genes.tre"
printf '(A,(B,C));\n' >"${WORK}/species-missing.tre"
printf '((A,B),(C,X));\n' >"${WORK}/species-unknown.tre"
printf '((A,B),(C,D));\n((A,C),(B,D));\n' >"${WORK}/species-multiple.tre"
printf '(A,B,C,D);\n' >"${WORK}/species-root-polytomy.tre"

for invalid in root-polytomy root-unary duplicate malformed empty; do
  expect_failure "$invalid" "${JAVA[@]}" --cpu -q -i "${WORK}/${invalid}.tre"
done
for invalid in species-missing species-unknown species-multiple species-root-polytomy; do
  expect_failure "$invalid" "${JAVA[@]}" --cpu -q -i "${WORK}/valid-genes.tre" \
    --score-species-tree "${WORK}/${invalid}.tre"
done
expect_failure unrooted-mode "${JAVA[@]}" --cpu -q --unrooted -i "${WORK}/valid-genes.tre"
expect_failure bad-method "${JAVA[@]}" --cpu -q -i "${WORK}/valid-genes.tre" --im I5
expect_failure bad-preset "${JAVA[@]}" --cpu -q -i "${WORK}/valid-genes.tre" --search-space S4
expect_failure bad-numeric "${JAVA[@]}" --cpu -q -i "${WORK}/valid-genes.tre" \
  --large-n-score-type decimal
expect_failure input-output-collision "${JAVA[@]}" --cpu -q \
  -i "${WORK}/valid-genes.tre" -o "${WORK}/valid-genes.tre"
printf '((A,B),(C,D));\n' >"${WORK}/protected-species.tre"
cp "${WORK}/protected-species.tre" "${WORK}/protected-species.expected"
expect_failure species-output-collision "${JAVA[@]}" --cpu -q \
  -i "${WORK}/valid-genes.tre" \
  --score-species-tree "${WORK}/protected-species.tre" \
  -o "${WORK}/protected-species.tre"
cmp "${WORK}/protected-species.expected" "${WORK}/protected-species.tre"
expect_failure species-log-collision "${JAVA[@]}" --cpu -q \
  -i "${WORK}/valid-genes.tre" \
  --score-species-tree "${WORK}/protected-species.tre" \
  --log-file "${WORK}/protected-species.tre"
cmp "${WORK}/protected-species.expected" "${WORK}/protected-species.tre"
printf '((A,B),(C,D));\n' >"${WORK}/protected-input.tre"
cp "${WORK}/protected-input.tre" "${WORK}/protected-input.expected"
expect_failure input-dump-collision "${JAVA[@]}" --cpu -q \
  -i "${WORK}/protected-input.tre" --dump-clusters "${WORK}/protected-input.tre"
cmp "${WORK}/protected-input.expected" "${WORK}/protected-input.tre"

printf '((A:0.1,B:2.0)95:0.3,(C:0.4,D:0.5)88:0.6);\n((A,C),(B,D));\n' \
  >"${WORK}/decorated-genes.tre"
printf '((A:3.0,B:4.0)77:1.0,(C:2.0,D:8.0)66:1.5);\n' \
  >"${WORK}/decorated-species.tre"
for method in I1 I2 I3 I4; do
  decorated_output="$("${JAVA[@]}" --cpu -q -i "${WORK}/decorated-genes.tre" \
    --score-species-tree "${WORK}/decorated-species.tre" --im "$method" 2>&1)"
  [[ "$(sed -n 's/^TRIPLET_SCORE: //p' <<<"$decorated_output" | tail -1)" == 4 ]]
done

NO_COLOR=1 "${ROOT}/stelarx" --no-build --version >"${WORK}/version.log" 2>&1
grep -q "STELAR-X  v" "${WORK}/version.log"
NO_COLOR=1 "${ROOT}/stelarx" --no-build --help >"${WORK}/help.log" 2>&1
grep -q -- "--intersection-method, --im" "${WORK}/help.log"

if [[ $PACKAGING -eq 1 ]]; then
  echo "Building and smoke-testing a self-contained CPU package..."
  "${ROOT}/build_portable.sh" --without-cuda --no-archive \
    --output-dir "${WORK}/dist" >"${WORK}/portable.log" 2>&1
  grep -q "Portable application ready:" "${WORK}/portable.log"
fi

gpu_available=0
if [[ "$GPU_MODE" != off ]]; then
  if [[ ! -f "${ROOT}/native/libstelarx_weight.so" ]]; then
    echo "Building native CUDA libraries..."
    "${ROOT}/build_native.sh" >/dev/null
  fi
  if "${JAVA[@]}" --gpu-strict --diagnose >"${WORK}/gpu-probe.log" 2>&1; then
    gpu_available=1
  elif [[ "$GPU_MODE" == require ]]; then
    echo "CUDA was required but strict diagnostics failed:" >&2
    tail -20 "${WORK}/gpu-probe.log" >&2
    exit 1
  else
    echo "CUDA unavailable; strict hardware layer skipped (use --require-gpu to make this fatal)."
  fi
fi

if [[ $gpu_available -eq 1 ]]; then
  echo "Running strict CUDA scoring, batching, numeric, inference, and parity tests..."
  if [[ $QUICK -eq 1 ]]; then
    STELARX_SKIP_BUILD=1 STELARX_GPU_RANDOM_CASES=2 \
      "${ROOT}/test/run_stelarx_gpu_tests.sh"
  else
    STELARX_SKIP_BUILD=1 STELARX_GPU_RANDOM_CASES=6 \
      "${ROOT}/test/run_stelarx_gpu_tests.sh"
  fi
  python3 "${ROOT}/test/test_similarity_matrix.py" --stelarx-root "$ROOT" \
    --mode gpu --seeds "${MATRIX_SEEDS[@]}"
  python3 "${ROOT}/test/test_distance_matrix.py" --stelarx-root "$ROOT" \
    --mode gpu --seeds "${MATRIX_SEEDS[@]}"
fi

echo "STELAR-X comprehensive suite: PASS ($((SECONDS - started)) seconds)"
