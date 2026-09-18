#!/bin/bash
# test_rf.sh — Run STELAR-X on a dataset (Mode 1 and Mode 2) and report RF vs true tree.
#
# Usage:
#   bash test_rf.sh 37          # test on 37-taxa dataset
#   bash test_rf.sh 48          # test on 48-taxa dataset
#   bash test_rf.sh 200         # test on 200-taxa dataset
#   bash test_rf.sh all         # test all three datasets in sequence

set -e
cd "$(dirname "$0")"

TAXA="${1:-37}"

run_one() {
    local N="$1"
    local INPUT="all_gt_bs_rooted_${N}.tre"
    local TRUE_TREE="true_${N}.tre"
    local OUT_M1="/tmp/stelarx_${N}_mode1.newick"
    local OUT_M2="/tmp/stelarx_${N}_mode2.newick"

    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "  STELAR-X  |  ${N}-taxa  |  GPU"
    echo "════════════════════════════════════════════════════════"

    if [ ! -f "$INPUT" ]; then
        echo "  [SKIP] Input file not found: $INPUT"
        return
    fi

    # ── Mode 1 (tree-local DP) ────────────────────────────────────────────────
    echo ""
    echo "  ── Mode 1 (local, tree-local DP) ──"
    START=$(date +%s%3N)
    java -Djava.library.path=native -cp build stelarx.Main \
        -i "$INPUT" -o "$OUT_M1" \
        --gpu -vv 2>&1
    END=$(date +%s%3N)
    TIME_M1=$(( END - START ))

    if [ -f "$TRUE_TREE" ]; then
        RF_M1=$(python3 rf.py "$TRUE_TREE" "$OUT_M1" 2>/dev/null | grep "Robinson-Foulds" | awk '{print $3}')
        SIM_M1=$(python3 rf.py "$TRUE_TREE" "$OUT_M1" 2>/dev/null | grep "similarity" | awk '{print $3}')
        echo "  RF distance : $RF_M1   |   Similarity: $SIM_M1"
    else
        echo "  [SKIP] True tree not found: $TRUE_TREE"
    fi
    echo "  Wall time   : ${TIME_M1} ms"

    # ── Mode 2 (cross-tree DP, full search) ───────────────────────────────────
    echo ""
    echo "  ── Mode 2 (full, cross-tree DP) ──"
    START=$(date +%s%3N)
    java -Djava.library.path=native -cp build stelarx.Main \
        -i "$INPUT" -o "$OUT_M2" \
        --gpu --search-mode full -vv 2>&1
    END=$(date +%s%3N)
    TIME_M2=$(( END - START ))

    if [ -f "$TRUE_TREE" ]; then
        RF_M2=$(python3 rf.py "$TRUE_TREE" "$OUT_M2" 2>/dev/null | grep "Robinson-Foulds" | awk '{print $3}')
        SIM_M2=$(python3 rf.py "$TRUE_TREE" "$OUT_M2" 2>/dev/null | grep "similarity" | awk '{print $3}')
        echo "  RF distance : $RF_M2   |   Similarity: $SIM_M2"
    fi
    echo "  Wall time   : ${TIME_M2} ms"

    # ── Summary ───────────────────────────────────────────────────────────────
    echo ""
    echo "  ┌─────────────────────────────────────────┐"
    printf  "  │  %-8s  RF        Similarity  Time     │\n" "${N}-taxa"
    echo    "  ├─────────────────────────────────────────┤"
    printf  "  │  Mode 1   %-8s  %-10s  %6s ms │\n" "$RF_M1" "$SIM_M1" "$TIME_M1"
    printf  "  │  Mode 2   %-8s  %-10s  %6s ms │\n" "$RF_M2" "$SIM_M2" "$TIME_M2"
    echo    "  └─────────────────────────────────────────┘"
}

if [ "$TAXA" = "all" ]; then
    run_one 37
    run_one 48
    run_one 200
else
    run_one "$TAXA"
fi

echo ""
