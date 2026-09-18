package stelarx.greedy;

/**
 * One polytomy that needs resolution.
 *
 * Carries the descriptor of the polytomy node within its snapshot consensus
 * tree plus the group ranges (children + optional "rest") that Step A and
 * Step B will operate on.
 *
 *   thresholdIndex — which T[ti] snapshot the polytomy lives in
 *   tree           — the snapshot tree
 *   node           — the polytomy node (children.size() ≥ 3)
 *   groupLos/His   — per-group [lo, hi) ranges into A_cons.  Length is
 *                    d for a virtual-root polytomy (children tile [0, n)),
 *                    d+1 for a non-root polytomy (the last range is the
 *                    "rest" = A_cons \ [v.lo, v.hi)).  The rest range is
 *                    represented as a single half-open interval ONLY when
 *                    v.lo == 0 OR v.hi == n; otherwise rest is two disjoint
 *                    sub-ranges [0, v.lo) and [v.hi, n).  In that case
 *                    groupLos[d] = 0, groupHis[d] = v.lo, and a second pair
 *                    is stored in {@code restLos2/restHis2}.
 *   degree         — d = number of children (excludes the implicit rest)
 */
public final class PolytomyTask {
    public final int thresholdIndex;
    public final ConsensusTree tree;
    public final ConsensusTree.SNode node;
    public final int degree;
    public final int numGroups;          // d (root polytomy) or d+1 (non-root)
    public final int[] groupLos;         // length numGroups
    public final int[] groupHis;         // length numGroups
    /** True iff the rest group is non-contiguous (two sub-ranges). */
    public final boolean restSplit;
    /** Second half of rest range when {@code restSplit} is true; else -1. */
    public final int restLos2;
    public final int restHis2;

    PolytomyTask(int thresholdIndex, ConsensusTree tree, ConsensusTree.SNode node,
                 int degree, int numGroups,
                 int[] groupLos, int[] groupHis,
                 boolean restSplit, int restLos2, int restHis2) {
        this.thresholdIndex = thresholdIndex;
        this.tree           = tree;
        this.node           = node;
        this.degree         = degree;
        this.numGroups      = numGroups;
        this.groupLos       = groupLos;
        this.groupHis       = groupHis;
        this.restSplit      = restSplit;
        this.restLos2       = restLos2;
        this.restHis2       = restHis2;
    }

    /** Estimated work for load balancing — Step A is O(|v|²) on the polytomy size. */
    public long estimatedCost() {
        int vSize = node.rangeSize();
        return (long) vSize * vSize + (long) degree * degree * (long) Math.max(1, (int) Math.log(degree + 1));
    }
}
