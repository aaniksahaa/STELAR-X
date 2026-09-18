#!/usr/bin/env bash
# run_tests.sh — Run all STELAR-X test cases.
#
# For each test input:
#   1. Run verify_weights.py (independent Python reimplementation) → expected score
#   2. Run STELAR-X binary with the configured flags → actual score
#   3. PASS only if both succeed and scores match
#   4. If a tc*_true.tre file exists, also compute RF distance on the STELAR-X output
#
# Usage:  bash test/run_tests.sh [TC_FILTER] [--search-mode local|full] [--cpu|--gpu]
#   TC_FILTER (optional): pattern to match test names, e.g. "tc1" or "tc1[23]"
#
# Exit code: 0 if all tests pass, 1 if any fail.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT_DIR="${SCRIPT_DIR}/input"
VERIFIER="${SCRIPT_DIR}/verify_weights.py"

BUILD_DIR="${ROOT_DIR}/build"
NATIVE_DIR="${ROOT_DIR}/native"

FILTER="tc[0-9]*"
COMPUTE_MODE="--gpu"
SEARCH_MODE="local"
SKIP_BUILD=0

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpu)          COMPUTE_MODE="--cpu"; shift ;;
        --gpu)          COMPUTE_MODE="--gpu"; shift ;;
        --search-mode)  SEARCH_MODE="$2"; shift 2 ;;
        --no-build)     SKIP_BUILD=1; shift ;;
        *)              FILTER="$1"; shift ;;
    esac
done

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

pass=0
fail=0

# Optional: pass extra STELAR-X flags via WEIGHT_METHOD env, e.g.
#   WEIGHT_METHOD=bitset bash test/run_tests.sh --gpu
EXTRA_OPTS=()
[[ -n "${WEIGHT_METHOD:-}" ]] && EXTRA_OPTS+=(--weight-intersection-method "$WEIGHT_METHOD")
# EXTRA passes arbitrary flags through, e.g. EXTRA="--no-anchor-outgroup"
[[ -n "${EXTRA:-}" ]] && EXTRA_OPTS+=($EXTRA)

STELARX_CMD=(java -Djava.library.path="$NATIVE_DIR" -cp "$BUILD_DIR" stelarx.Main
             $COMPUTE_MODE --search-mode "$SEARCH_MODE" "${EXTRA_OPTS[@]}")

run_tc () {
    local input="$1"
    local name
    name="$(basename "$input" .tre)"
    local tc_id="${name%%_*}"

    [[ "$name" == *_true ]] && return 0

    printf "  %-45s" "$name"

    # ── Step 1: Python verifier → expected score ──────────────────────────────
    local verifier_out expected_score
    if ! verifier_out="$(python3 "$VERIFIER" "$input" 2>&1)"; then
        printf "${RED}FAIL${NC}  verifier error\n"
        echo "$verifier_out" | tail -3 | sed 's/^/    /'
        ((fail++)) || true
        return 0
    fi
    expected_score=$(echo "$verifier_out" | grep -oP 'quartet score = \K[0-9]+' || true)
    if [[ -z "$expected_score" ]]; then
        printf "${RED}FAIL${NC}  verifier gave no score\n"
        ((fail++)) || true
        return 0
    fi

    # ── Step 2: STELAR-X binary → actual score ────────────────────────────────
    # Polytomous fixtures need FULL search to match the full-DP oracle (local search
    # finds a valid but lower optimum on the richer polytomous candidate set).
    local cmd=("${STELARX_CMD[@]}")
    if [[ "$name" == *polytomy* && "$SEARCH_MODE" != "full" ]]; then
        cmd=(java -Djava.library.path="$NATIVE_DIR" -cp "$BUILD_DIR" stelarx.Main
             $COMPUTE_MODE --search-mode full "${EXTRA_OPTS[@]}")
    fi
    # The independent oracle scores unresolved quartets. Keep native polytomies
    # for these fixtures; default CLI refinement is covered separately.
    if [[ "$name" == *polytomy* ]]; then
        cmd+=(--keep-polytomy-during-inference)
    fi
    local stelarx_out actual_score stelarx_tree
    if ! stelarx_out="$("${cmd[@]}" -i "$input" 2>&1)"; then
        printf "${RED}FAIL${NC}  STELAR-X crashed\n"
        echo "$stelarx_out" | tail -5 | sed 's/^/    /'
        ((fail++)) || true
        return 0
    fi
    actual_score=$(echo "$stelarx_out" | grep -oP 'Final quartet score = \K[0-9]+' || true)
    stelarx_tree=$(echo "$stelarx_out" | grep -v '^\[' | grep -v '^[[:space:]]' | grep ';' | tail -1 || true)

    if [[ -z "$actual_score" ]]; then
        printf "${RED}FAIL${NC}  STELAR-X gave no score\n"
        echo "$stelarx_out" | tail -5 | sed 's/^/    /'
        ((fail++)) || true
        return 0
    fi

    # ── Step 3: Compare scores ────────────────────────────────────────────────
    if [[ "$actual_score" != "$expected_score" ]]; then
        printf "${RED}FAIL${NC}  score mismatch: expected=${expected_score} got=${actual_score}\n"
        ((fail++)) || true
        return 0
    fi

    # ── Step 4: RF distance (optional) ───────────────────────────────────────
    local rf_str=""
    local true_file="${INPUT_DIR}/${tc_id}_true.tre"
    if [[ -f "$true_file" && -n "$stelarx_tree" ]]; then
        local true_newick
        true_newick="$(cat "$true_file")"
        local rf_out
        rf_out="$(python3 "$VERIFIER" "$input" --compare "$true_newick" 2>&1)" || true
        local rf
        rf=$(echo "$rf_out" | grep -oP 'RF = \K[0-9]+' || true)
        [[ -n "$rf" ]] && rf_str="RF=${rf}"
    fi

    printf "${GREEN}PASS${NC}  score=%-10s %s\n" "$actual_score" "$rf_str"
    ((pass++)) || true
}

echo -e "\n${BOLD}=== STELAR-X Test Suite ===${NC}"
printf "  mode: %s  search: %s\n\n" "$COMPUTE_MODE" "$SEARCH_MODE"

# ── Build ─────────────────────────────────────────────────────────────────────
if [[ $SKIP_BUILD -eq 0 ]]; then
    echo "  Building Java..."
    bash "$ROOT_DIR/build.sh" > /dev/null \
        && echo "  Build Java   OK" \
        || { echo -e "  ${RED}Build Java   FAILED${NC}"; exit 1; }
    echo "  Building native (CUDA)..."
    bash "$ROOT_DIR/build_native.sh" > /dev/null \
        && echo "  Build native OK" \
        || { echo -e "  ${RED}Build native FAILED${NC}"; exit 1; }
    echo
fi

mapfile -t inputs < <(
    find "$INPUT_DIR" -maxdepth 1 -name "${FILTER}_*.tre" ! -name "*_true.tre" | sort
)

if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "No test inputs matched filter '${FILTER}' in ${INPUT_DIR}/"
    exit 1
fi

for f in "${inputs[@]}"; do
    run_tc "$f"
done

echo
echo -e "  Results: ${GREEN}${pass} passed${NC}  ${RED}${fail} failed${NC}"
echo

[[ $fail -eq 0 ]]
