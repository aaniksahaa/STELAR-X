#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

# Source mode loads only the uppercase already-completed list and its predicate;
# run-bulk-simulated.sh must not start any simulation or inference work.
source "${ROOT}/run-bulk-simulated.sh"

# The committed list ships empty: every replicate of every configuration runs.
[[ ${#ALREADY_COMPLETED_SIMULATED_CONFIGS[@]} -eq 0 ]] || \
  fail "ALREADY_COMPLETED_SIMULATED_CONFIGS should be empty in the committed script"

if IS_SIMULATED_CONFIG_ALREADY_COMPLETED 75000 1000 0.000001 100000 200000 R5; then
  fail "a replicate was reported as already completed with an empty list"
fi

# Add a test-only dummy tuple and drive the same branch used by the production
# replicate loops. Exactly the dummy R2 run must be skipped as already done;
# its neighbours, other gene-tree counts and other SB settings still run.
DUMMY_CONFIG="42,7,0.125,11,22,R2"
ALREADY_COMPLETED_SIMULATED_CONFIGS+=("$DUMMY_CONFIG")

if IS_SIMULATED_CONFIG_ALREADY_COMPLETED 42 8 0.125 11 22 R2; then
  fail "an already-completed entry leaked into another gene-tree count"
fi
if IS_SIMULATED_CONFIG_ALREADY_COMPLETED 42 7 0.250 11 22 R2; then
  fail "an already-completed entry leaked into another SB setting"
fi

EXECUTED_REPLICATES=()
SKIPPED_REPLICATES=()
for REPLICATE in R1 R2 R3; do
  if IS_SIMULATED_CONFIG_ALREADY_COMPLETED 42 7 0.125 11 22 "$REPLICATE"; then
    SKIPPED_REPLICATES+=("$REPLICATE")
    continue
  fi
  EXECUTED_REPLICATES+=("$REPLICATE")
done

[[ "${SKIPPED_REPLICATES[*]}" == "R2" ]] || \
  fail "dummy already-completed entry did not skip exactly R2: ${SKIPPED_REPLICATES[*]}"
[[ "${EXECUTED_REPLICATES[*]}" == "R1 R3" ]] || \
  fail "dummy already-completed entry suppressed a neighbouring run: ${EXECUTED_REPLICATES[*]}"

echo "PASS: bulk-simulated already-completed replicate list (empty by default, dummy skip)"
