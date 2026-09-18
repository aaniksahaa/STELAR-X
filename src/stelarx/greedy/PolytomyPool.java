package stelarx.greedy;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Collects every polytomy across all 7 snapshot trees, applies the
 * legacy size-limit budget {@code N = 50 + n·25}, and exposes the
 * accepted polytomies as a single work pool (per the doc §8.2 + §10.1).
 *
 * Size-limit rule (matches ASTRAL-MP's loop at WQDataCollection.java
 * lines 1115-1144 exactly):
 *   1. Collect ALL polytomy child-counts across all 7 snapshots.
 *   2. Sort ascending.
 *   3. Walk in order; accumulate Σ degree².  Each iteration starts by
 *      checking accumulated &lt; N — if true, ADD this degree² and advance.
 *      Stop when the next iteration's check would fail.
 *   4. {@code sizeLimit = last degree included}; or 3 if none.
 *   5. A polytomy of degree d is processed iff d ≤ sizeLimit.
 *
 * NOTE on ASTRAL-MP parity: ASTRAL-MP's loop (WQDataCollection.java:1154) actually
 * processes EVERY polytomy; its polytomySizeLimit only disables the *quadratic*
 * NN-balls (line 1247) for over-limit polytomies — the linear path still runs.
 * By default we conservatively drop d > sizeLimit (smaller X, cheaper).  Setting
 * --stepb-process-large-polytomies lifts that, processing every polytomy via the
 * linear path (quadratic still skipped for the large ones), matching ASTRAL-MP.
 */
public final class PolytomyPool {

    public static final int BUDGET_MIN  = 50;       // GREEDY_ADDITION_MAX_POLYTOMY_MIN
    public static final int BUDGET_MULT = 25;       // GREEDY_ADDITION_MAX_POLYTOMY_MULT

    public final List<PolytomyTask> tasks;          // accepted polytomies (degree ≤ sizeLimit)
    public final int sizeLimit;
    public final int numTotal;                      // total polytomies (before size filter)
    public final int numAccepted;
    public final int numSkipped;
    public final int budgetN;
    public final long sumSquaresAccumulated;
    public final int[] degreeHistogram;             // index = degree; capped at min(maxDeg+1, n+1)

    private PolytomyPool(List<PolytomyTask> tasks, int sizeLimit, int numTotal,
                         int numAccepted, int numSkipped, int budgetN,
                         long sumSq, int[] hist) {
        this.tasks                 = tasks;
        this.sizeLimit             = sizeLimit;
        this.numTotal              = numTotal;
        this.numAccepted           = numAccepted;
        this.numSkipped            = numSkipped;
        this.budgetN               = budgetN;
        this.sumSquaresAccumulated = sumSq;
        this.degreeHistogram       = hist;
    }

    /**
     * @param snapshots          7 snapshot trees (some may be null in unusual cases)
     * @param numTaxa            n
     * @param explicitLimit      override; ≤ 0 means "compute from budget"
     */
    public static PolytomyPool build(ConsensusTree[] snapshots, int numTaxa,
                                      int explicitLimit) {
        // ── Pass 1: collect all polytomy descriptors (degree only) ──
        List<int[]> degByThreshold = new ArrayList<>();   // parallel: per-snapshot degree lists
        List<Integer> allDegrees = new ArrayList<>();
        int maxDeg = 0;
        for (int ti = 0; ti < snapshots.length; ti++) {
            ConsensusTree ct = snapshots[ti];
            if (ct == null) { degByThreshold.add(new int[0]); continue; }
            List<Integer> here = new ArrayList<>();
            ct.forEachInternalNode(node -> {
                // ASTRAL-MP's loop iterates EVERY internal node, root
                // included.  A polytomy at the virtual root is just the
                // top-level polytomy of the unrooted tree.
                int d = node.children.size();
                if (d > 2) here.add(d);
            });
            int[] arr = new int[here.size()];
            for (int j = 0; j < arr.length; j++) arr[j] = here.get(j);
            degByThreshold.add(arr);
            for (int d : arr) { allDegrees.add(d); if (d > maxDeg) maxDeg = d; }
        }
        int numTotal = allDegrees.size();

        // ── Compute size limit ──
        int sizeLimit;
        int budgetN = BUDGET_MIN + numTaxa * BUDGET_MULT;
        long sumSq = 0L;
        if (explicitLimit > 0) {
            sizeLimit = explicitLimit;
        } else if (allDegrees.isEmpty()) {
            sizeLimit = 3;   // fully binary; no polytomies anyway
        } else {
            List<Integer> sorted = new ArrayList<>(allDegrees);
            Collections.sort(sorted);
            int i = 0;
            while (i < sorted.size() && sumSq < budgetN) {
                int d = sorted.get(i);
                sumSq += (long) d * d;
                i++;
            }
            sizeLimit = (i > 0) ? sorted.get(i - 1) : 3;
        }

        // ── Pass 2: build PolytomyTask for accepted polytomies ──
        // "Lift the bar": when --stepb-process-large-polytomies is set, accept
        // EVERY polytomy regardless of degree (ASTRAL-MP never drops a polytomy;
        // its size limit only disables the quadratic NN-balls).  Default: drop
        // d > sizeLimit, exactly as before.
        final boolean processAll =
            stelarx.Config.getInstance().isStepBProcessLargePolytomies();
        final int finalSizeLimit = sizeLimit;
        List<PolytomyTask> accepted = new ArrayList<>();
        int[] hist = new int[Math.max(8, maxDeg + 1)];
        int numAccepted = 0, numSkipped = 0;

        for (int ti = 0; ti < snapshots.length; ti++) {
            ConsensusTree ct = snapshots[ti];
            if (ct == null) continue;
            ct.forEachInternalNode(node -> {
                int d = node.children.size();
                if (d <= 2) return;
                if (d < hist.length) hist[d]++;
            });
        }

        for (int ti = 0; ti < snapshots.length; ti++) {
            ConsensusTree ct = snapshots[ti];
            if (ct == null) continue;
            final int tii = ti;
            List<PolytomyTask> localAccepted = new ArrayList<>();
            List<Integer> localSkipped = new ArrayList<>();
            ct.forEachInternalNode(node -> {
                int d = node.children.size();
                if (d <= 2) return;
                if (!processAll && d > finalSizeLimit) {
                    localSkipped.add(d);
                    return;
                }
                localAccepted.add(makeTask(tii, ct, node, numTaxa));
            });
            accepted.addAll(localAccepted);
            numAccepted += localAccepted.size();
            numSkipped += localSkipped.size();
        }

        return new PolytomyPool(accepted, sizeLimit, numTotal,
                                 numAccepted, numSkipped, budgetN, sumSq, hist);
    }

    /** Build the {@link PolytomyTask} for one polytomy node. */
    private static PolytomyTask makeTask(int ti, ConsensusTree ct,
                                          ConsensusTree.SNode node, int numTaxa) {
        int d = node.children.size();
        int vLo = node.rangeLo(), vHi = node.rangeHi();
        boolean isRootPolytomy = (vLo == 0 && vHi == numTaxa);

        int groups = isRootPolytomy ? d : d + 1;
        int[] los = new int[groups];
        int[] his = new int[groups];
        for (int i = 0; i < d; i++) {
            ConsensusTree.SNode c = node.children.get(i);
            los[i] = c.rangeLo();
            his[i] = c.rangeHi();
        }

        boolean restSplit = false;
        int restLos2 = -1, restHis2 = -1;
        if (!isRootPolytomy) {
            // Rest = [0, vLo) ∪ [vHi, n).
            // If one of the two halves is empty, the rest is a single contiguous
            // range; otherwise we mark restSplit and store the second half
            // separately so callers can use combineDisjointSigma{1,2}.
            if (vLo == 0) {                      // single range [vHi, n)
                los[d] = vHi; his[d] = numTaxa;
            } else if (vHi == numTaxa) {         // single range [0, vLo)
                los[d] = 0;   his[d] = vLo;
            } else {                             // split
                los[d] = 0;   his[d] = vLo;
                restSplit = true;
                restLos2 = vHi;
                restHis2 = numTaxa;
            }
        }

        return new PolytomyTask(ti, ct, node, d, groups,
                                 los, his, restSplit, restLos2, restHis2);
    }

    /** Sort tasks by descending estimated cost — longest-first work-stealing (§10.3). */
    public List<PolytomyTask> tasksLongestFirst() {
        List<PolytomyTask> sorted = new ArrayList<>(tasks);
        sorted.sort((a, b) -> Long.compare(b.estimatedCost(), a.estimatedCost()));
        return sorted;
    }
}
