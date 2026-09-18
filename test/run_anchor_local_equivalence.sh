#!/usr/bin/env bash
# Compare local-mode rooted-triplet scores with and without anchor-free X.
#
# Unlike the general test runner, this keeps every fixture in true LOCAL mode,
# including polytomous cases. It is a regression check for leaf-induced Type-2
# rotations and anchored-outgroup search-space reduction.
#
# Usage:
#   bash test/run_anchor_local_equivalence.sh [TC_FILTER] [--cpu|--gpu] [--no-build]
# Optional: WEIGHT_METHOD=prefix-sum|smaller-side-traversal|bitset|simple-tree-walk

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT_DIR="$SCRIPT_DIR/input"

FILTER="tc[0-9]*"
COMPUTE_MODE="--cpu"
SKIP_BUILD=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpu)      COMPUTE_MODE="--cpu"; shift ;;
        --gpu)      COMPUTE_MODE="--gpu"; shift ;;
        --no-build) SKIP_BUILD=1; shift ;;
        *)          FILTER="$1"; shift ;;
    esac
done

if [[ $SKIP_BUILD -eq 0 ]]; then
    bash "$ROOT_DIR/build.sh" >/dev/null
    if [[ "$COMPUTE_MODE" == "--gpu" ]]; then
        bash "$ROOT_DIR/build_native.sh" >/dev/null
    fi
fi

EXTRA_OPTS=()
if [[ -n "${WEIGHT_METHOD:-}" ]]; then
    EXTRA_OPTS+=(--weight-intersection-method "$WEIGHT_METHOD")
fi
BASE_CMD=(java -Djava.library.path="$ROOT_DIR/native" -cp "$ROOT_DIR/build"
          stelarx.Main "$COMPUTE_MODE" --search-mode local "${EXTRA_OPTS[@]}")

score_of() {
    local input="$1"
    shift
    local output score
    output="$("${BASE_CMD[@]}" "$@" -i "$input" 2>&1)"
    score="$(sed -n 's/.*Final triplet score = \([0-9][0-9]*\).*/\1/p' <<<"$output" | tail -1)"
    if [[ -z "$score" ]]; then
        echo "ERROR: no score for $input" >&2
        echo "$output" | tail -8 >&2
        return 1
    fi
    printf '%s' "$score"
}

mapfile -t inputs < <(
    find "$INPUT_DIR" -maxdepth 1 -name "${FILTER}_*.tre" ! -name "*_true.tre" | sort
)

if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "No test inputs matched '$FILTER' in $INPUT_DIR" >&2
    exit 1
fi

pass=0
fail=0
printf '\n=== Anchored local equivalence (%s) ===\n\n' "$COMPUTE_MODE"
for input in "${inputs[@]}"; do
    name="$(basename "$input" .tre)"
    base="$(score_of "$input" --no-anchor-outgroup)"
    anchored="$(score_of "$input")"
    if [[ "$base" == "$anchored" ]]; then
        printf '  %-43s PASS  score=%s\n' "$name" "$base"
        ((pass++)) || true
    else
        printf '  %-43s FAIL  local=%s anchored=%s\n' "$name" "$base" "$anchored"
        ((fail++)) || true
    fi
done

printf '\nResults: %d passed, %d failed\n\n' "$pass" "$fail"
[[ $fail -eq 0 ]]
