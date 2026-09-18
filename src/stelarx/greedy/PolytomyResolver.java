package stelarx.greedy;

import stelarx.Config;
import stelarx.cluster.ClusterHash;
import stelarx.completion.EulerTourBuilder;
import stelarx.completion.SimilarityMatrix;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.Threading;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Polytomy resolution: per-polytomy Step A (UPGMA on the group similarity
 * matrix → bipartitions) and Step B (sampleAndResolve, in a follow-up).
 *
 * Step A — design §8.3:
 *   1. Build a g×g group-similarity matrix where entry (i,j) =
 *      average sim[x][y] over x ∈ group_i, y ∈ group_j (id-indexed lookups).
 *      Self-similarity (i,i) is left at 0.  The g×g fill dominates (O(|v|²)).
 *   2. Run UPGMA on that matrix.  Every non-root internal dendrogram node
 *      defines a bipartition (groups in its subtree | the rest).
 *   3. For each emission, compute the side's multi-range (union of selected
 *      groups' consensus-tree ranges, plus the rest group's split sub-ranges
 *      when the rest is selected), pick the smaller side, compute the
 *      double-hash signature via the consensus-tree prefix arrays, and add
 *      to the emission buffer (dedup by signature).
 */
public final class PolytomyResolver {

    private PolytomyResolver() {}

    /**
     * Drive Step A + Step B for every polytomy in the pool, in parallel.
     *
     * Dispatch is LPT-first: tasks are sorted by {@link PolytomyTask#estimatedCost}
     * descending and submitted individually to the {@link Threading} executor.
     * The fixed-thread-pool's internal queue picks them up dynamically, so the
     * costly polytomies start first and short ones fill the gaps — effectively
     * "work-stealing-lite" via the standard executor without any custom queue.
     *
     * Per-task RNG is seeded deterministically from the threshold and the
     * polytomy's order-independent taxon-set signature, so the emission set is
     * reproducible across runs and threading configurations.
     *
     * NO nested {@code Threading.processRangeParallel}: Step A's UPGMA and
     * Step B's resolveByDistance use {@link MiniUPGMA}, which is sequential.
     * This is what makes the outer parallelism safe — otherwise blocked
     * polytomy workers would prevent inner sub-tasks from ever starting.
     *
     * Each task writes to a private {@link EmissionBuffer}.  This is essential:
     * Step B decides how many adaptive rounds to run from the number of signatures
     * added by that task.  Measuring a shared concurrent buffer makes that decision
     * depend on which unrelated task happens to publish first.  Each completed task
     * merges immediately into the shared buffer; deterministic duplicate selection
     * there keeps the final set and provenance schedule-independent without retaining
     * every completed task map behind the slowest LPT task.
     */
    public static int runAllParallel(List<PolytomyTask> tasks,
                                      List<Tree> geneTrees, SimilarityMatrix sim,
                                      EmissionBuffer buffer, int numTaxa,
                                      long baseSeed) {
        if (tasks.isEmpty()) return 0;

        // Step B per-gene-tree restriction route: O(d log d) auxiliary tree via
        // Euler-tour LCA (default) vs the O(n) full walk. The Euler structures are
        // immutable and shared read-only across all parallel tasks; build once.
        final boolean fastRestriction = Config.getInstance().isStepBFastRestriction();
        final EulerTourBuilder.TourData[] tours;
        if (fastRestriction) {
            tours = new EulerTourBuilder.TourData[geneTrees.size()];
            for (int i = 0; i < geneTrees.size(); i++)
                tours[i] = EulerTourBuilder.build(geneTrees.get(i), numTaxa);
        } else {
            tours = null;
        }

        List<PolytomyTask> sorted = new ArrayList<>(tasks);
        sorted.sort((a, b) -> Long.compare(b.estimatedCost(), a.estimatedCost()));

        List<Future<Integer>> futures = new ArrayList<>(sorted.size());
        for (PolytomyTask task : sorted) {
            final PolytomyTask t = task;
            int lo = t.node.rangeLo(), hi = t.node.rangeHi();
            long signatureSeed = t.tree.sigma1(0, lo, hi)
                               ^ Long.rotateLeft(t.tree.sigma2(0, lo, hi), 29);
            long seed = baseSeed ^ ((long) t.thresholdIndex << 32)
                       ^ signatureSeed ^ ((long) t.node.rangeSize() << 1);
            futures.add(Threading.submit(() -> {
                Random rng = new Random(seed);
                EmissionBuffer local = new EmissionBuffer();
                if (sim != null) stepA(t, sim, local, numTaxa);
                stepB(t, geneTrees, tours, local, numTaxa, rng, sim);
                int added = 0;
                for (EmittedBipartition emission : local.all()) {
                    if (buffer.add(emission)) added++;
                }
                return added;
            }));
        }

        int total = 0;
        for (Future<Integer> future : futures) {
            try {
                total += future.get();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new RuntimeException(e);
            } catch (ExecutionException e) {
                throw new RuntimeException(e.getCause());
            }
        }
        return total;
    }

    /** Stable insertion-order accessor for the polytomy list — used by callers
     *  that want a fully deterministic single-thread pass for diffing. */
    public static List<PolytomyTask> sortLPT(List<PolytomyTask> tasks) {
        List<PolytomyTask> s = new ArrayList<>(tasks);
        Collections.sort(s, (a, b) -> Long.compare(b.estimatedCost(), a.estimatedCost()));
        return s;
    }

    /**
     * Run Step A for one polytomy.  Returns the number of new bipartitions
     * added to the buffer (signatures not previously seen).
     */
    public static int stepA(PolytomyTask task, SimilarityMatrix sim,
                            EmissionBuffer buffer, int numTaxa) {
        int g = task.numGroups;
        if (g < 3) return 0;  // need ≥ 3 groups for a non-trivial bipartition

        ConsensusTree ct = task.tree;
        int[] aCons = ct.aCons();

        // ── (1) Build g×g group similarity matrix ───────────────────────────
        double[] groupSim = new double[g * g];
        for (int i = 0; i < g; i++) {
            int iLo = task.groupLos[i], iHi = task.groupHis[i];
            int iLo2 = -1, iHi2 = -1;
            if (i == g - 1 && task.numGroups > task.degree && task.restSplit) {
                iLo2 = task.restLos2; iHi2 = task.restHis2;
            }
            int iSize = (iHi - iLo) + (iLo2 < 0 ? 0 : (iHi2 - iLo2));

            for (int j = i + 1; j < g; j++) {
                int jLo = task.groupLos[j], jHi = task.groupHis[j];
                int jLo2 = -1, jHi2 = -1;
                if (j == g - 1 && task.numGroups > task.degree && task.restSplit) {
                    jLo2 = task.restLos2; jHi2 = task.restHis2;
                }
                int jSize = (jHi - jLo) + (jLo2 < 0 ? 0 : (jHi2 - jLo2));

                double s = avgSim(sim, aCons, iLo, iHi, iLo2, iHi2,
                                  jLo, jHi, jLo2, jHi2);
                // Average over (iSize × jSize) pairs
                if (iSize > 0 && jSize > 0) s /= ((double) iSize * (double) jSize);
                groupSim[i * g + j] = s;
                groupSim[j * g + i] = s;
            }
        }

        // ── (2) Run UPGMA on the g×g matrix ─────────────────────────────────
        // Large polytomy (g > 31) under the lift-the-bar flag → exact O(d²) NN-chain
        // (build()'s O(g³) closest-pair scan would dominate). g ≤ 31 keeps the
        // proven O(g³) path so existing emissions are byte-identical.
        Tree dendrogram = (g > 31 && Config.getInstance().isStepBProcessLargePolytomies())
            ? MiniUPGMA.buildFast(groupSim, g, /*treeIndex*/0)
            : MiniUPGMA.build(groupSim, g, /*treeIndex*/0);

        // ── (3) Walk dendrogram; emit each non-root internal node ──────────
        int[] addedNew = {0};
        walkDendrogram(dendrogram.root, dendrogram.postorderArray,
                       task, buffer, numTaxa, addedNew);
        return addedNew[0];
    }

    // ── Group similarity: sum over pairs (without averaging — caller divides) ──

    private static double avgSim(SimilarityMatrix sim, int[] aCons,
                                  int iLo, int iHi, int iLo2, int iHi2,
                                  int jLo, int jHi, int jLo2, int jHi2) {
        double s = 0;
        s += pairSum(sim, aCons, iLo, iHi, jLo, jHi);
        if (jLo2 >= 0) s += pairSum(sim, aCons, iLo, iHi, jLo2, jHi2);
        if (iLo2 >= 0) {
            s += pairSum(sim, aCons, iLo2, iHi2, jLo, jHi);
            if (jLo2 >= 0) s += pairSum(sim, aCons, iLo2, iHi2, jLo2, jHi2);
        }
        return s;
    }

    private static double pairSum(SimilarityMatrix sim, int[] aCons,
                                   int iLo, int iHi, int jLo, int jHi) {
        double s = 0;
        for (int ai = iLo; ai < iHi; ai++) {
            int x = aCons[ai];
            for (int aj = jLo; aj < jHi; aj++) {
                int y = aCons[aj];
                s += sim.getSim(x, y);
            }
        }
        return s;
    }

    // ── Dendrogram walk: every non-root internal node = one bipartition ──

    private static void walkDendrogram(TreeNode dn, int[] postArr,
                                        PolytomyTask task, EmissionBuffer buffer,
                                        int numTaxa, int[] addedNew) {
        if (dn.isLeaf()) return;
        walkDendrogram(dn.left,  postArr, task, buffer, numTaxa, addedNew);
        walkDendrogram(dn.right, postArr, task, buffer, numTaxa, addedNew);
        if (dn.isRoot()) return;

        int nGroups = task.numGroups;
        boolean[] selected = new boolean[nGroups];
        int selectedSize = 0;
        for (int p = dn.rangeStart; p < dn.rangeEnd; p++) {
            int gi = postArr[p];
            if (!selected[gi]) {
                selected[gi] = true;
                selectedSize += groupSize(task, gi);
            }
        }
        int complementSize = numTaxa - selectedSize;
        if (selectedSize <= 1 || complementSize <= 1) return;        // trivial
        if (selectedSize == numTaxa) return;                          // whole tree

        // Pick smaller side; build its multi-range
        MultiRange canonical;
        int canonicalSize;
        if (selectedSize <= complementSize) {
            canonical = buildSideMultiRange(task, selected, /*inverted=*/false);
            canonicalSize = selectedSize;
        } else {
            canonical = buildSideMultiRange(task, selected, /*inverted=*/true);
            canonicalSize = complementSize;
        }

        // Compute double-hash signature via consensus prefix scan
        int m = task.tree.numSeeds();
        long[] sums = new long[m];
        long[] xors = new long[m];
        for (int s = 0; s < m; s++) {
            sums[s] = task.tree.combineDisjointSigma1(s, canonical.los, canonical.his);
            xors[s] = task.tree.combineDisjointSigma2(s, canonical.los, canonical.his);
        }
        ClusterHash sig = new ClusterHash(sums, xors, canonicalSize, m);

        if (buffer.add(new EmittedBipartition(
                sig, canonical, canonicalSize, 'A', task.thresholdIndex))) {
            addedNew[0]++;
        }
    }

    /** Number of taxa in group {@code gi}, accounting for the split-rest case. */
    private static int groupSize(PolytomyTask task, int gi) {
        int sz = task.groupHis[gi] - task.groupLos[gi];
        boolean isRest = (gi == task.numGroups - 1) && (task.numGroups > task.degree);
        if (isRest && task.restSplit) sz += task.restHis2 - task.restLos2;
        return sz;
    }

    // ── §8.4 Step B: sampleAndResolve with d-rep restriction ─────────────

    /** Legacy tuning constants — matches WQDataCollection.java lines 61-67. */
    public static final int STEPB_DEFAULT_RUNS      = 10;
    public static final int STEPB_MAX               = 100;
    public static final int STEPB_IMPROVEMENT_REWARD = 2;
    public static final int STEPB_MIN_FREQ          = 5;
    public static final double STEPB_MIN_RATIO      = 0.01;
    /** Quadratic NN-ball emission fires only for thresholdIndex < this (ASTRAL-MP
     *  GREEDY_DIST_ADDITTION_LAST_THRESHOLD_INDX = 3): the loosest 3 thresholds. */
    public static final int STEPB_QUADRATIC_MAX_THRESHOLD_INDEX = 3;

    /**
     * Step B — sampleAndResolve.  Runs {@link #STEPB_DEFAULT_RUNS} base rounds
     * and up to {@link #STEPB_MAX} adaptive bonus rounds (each productive round
     * adds {@link #STEPB_IMPROVEMENT_REWARD} more), per the legacy adaptive
     * scheme.  As in ASTRAL-MP, a round is productive only when it adds a new
     * accepted cluster whose support is greater than {@link #STEPB_MIN_FREQ}
     * and at least {@link #STEPB_MIN_RATIO} of the gene trees.
     *
     * Each round:
     *   1. Pick one random representative taxon per group.
     *   2. For each gene tree, walk it postorder; propagate a per-node int
     *      bitmap of "present reps in subtree".  When a node's popcount is in
     *      [2, presentCount - 1], it defines an induced split — map the rep
     *      bitmap back to group indices, build the multi-range, compute the
     *      smaller-side signature, emit.
     *
     * Per-tree walk is O(n) — fine for typical inputs.  An O(d log n) variant
     * via marked-ancestor walks or precomputed LCA is a future optimization.
     *
     * Limitation: assumes {@code numGroups ≤ 31} so the rep-membership bitmap
     * fits in an int.  The legacy size-limit budget already caps polytomy
     * degree near √(50 + 25n) ≪ 31 for typical n; we assert anyway.
     */
    public static int stepB(PolytomyTask task, List<Tree> geneTrees,
                            EulerTourBuilder.TourData[] tours,
                            EmissionBuffer buffer, int numTaxa, Random rng,
                            SimilarityMatrix sim) {
        int d = task.numGroups;
        if (d < 4) return 0;                // need ≥ 4 groups for non-trivial induced split
        // d > 31 needs the long[]-bitmap path; only enabled by the lift-the-bar flag.
        if (d > 31 && !Config.getInstance().isStepBProcessLargePolytomies()) return 0;

        // One reusable d×d induced-similarity buffer per task for resolveByDistance,
        // shared across all rounds (avoids re-allocating a d²-double matrix every round —
        // the remaining churn after hash-counts).  null for the int path (d≤31, tiny).
        double[] simBuf = (d > 31) ? new double[d * d] : null;

        int totalNewSignatures = 0;
        int adaptBonus = 0;
        int j = 0;
        while (j < STEPB_DEFAULT_RUNS + adaptBonus) {
            int beforeSize = buffer.size();
            boolean productive = stepBRound(
                task, geneTrees, tours, buffer, numTaxa, rng, sim, j, simBuf);
            int newThisRound = buffer.size() - beforeSize;
            totalNewSignatures += newThisRound;
            if (productive && adaptBonus < STEPB_MAX) {
                adaptBonus += STEPB_IMPROVEMENT_REWARD;
            }
            j++;
        }
        return totalNewSignatures;
    }

    /**
     * Run one Step B round, matching the legacy {@code sampleAndResolve} →
     * {@code resolveLinearly} flow:
     *   1. Pick one rep per group.
     *   2. For each gene tree, walk postorder and emit each non-root internal
     *      node's rep-bitmap if it passes the size filter (cnt ∈ [2, d-2]).
     *      This mirrors {@code Utils.getBitsets}.
     *   3. Aggregate bitmaps into a frequency map; each bipartition is counted
     *      under whichever side was encountered first (complementary dedupe
     *      against the existing key — same shape as
     *      {@code returnBitSetCounts}).
     *   4. Sort by frequency descending (with deterministic tie-break).
     *   5. Run a mini-greedy laminar build on the d reps: accept a bipartition
     *      iff it is pairwise nested-or-disjoint with every previously
     *      accepted bipartition.  ASTRAL-MP's exact buildTreeFromClusters
     *      check additionally requires the new cluster to MOVE ≥ 2 of the
     *      LCA's children — but exact duplicates (single-child matches) are
     *      already caught by the global signature dedup in {@link EmissionBuffer},
     *      so pairwise laminar gives an equivalent emission set here.
     *   6. Each accepted bipartition → emit full-taxa bipartition via
     *      smaller-side selection, hashed via the consensus prefix-scan.
     */
    private static boolean stepBRound(PolytomyTask task, List<Tree> geneTrees,
                                      EulerTourBuilder.TourData[] tours,
                                      EmissionBuffer buffer, int numTaxa, Random rng,
                                      SimilarityMatrix sim, int roundIndex, double[] simBuf) {
        int d = task.numGroups;
        // Large polytomy: rep-subsets no longer fit in an int → long[] path.
        if (d > 31) {
            return stepBRoundLong(
                task, geneTrees, tours, buffer, numTaxa, rng, sim, roundIndex, simBuf);
        }
        int[] aCons = task.tree.aCons();
        int allBits = (1 << d) - 1;

        int[] reps = new int[d];
        for (int gi = 0; gi < d; gi++) {
            int firstLo = task.groupLos[gi], firstHi = task.groupHis[gi];
            int firstSize = firstHi - firstLo;
            int totalSize = firstSize;
            boolean isRest = (gi == d - 1) && (task.numGroups > task.degree);
            int restExtraSize = 0;
            if (isRest && task.restSplit) {
                restExtraSize = task.restHis2 - task.restLos2;
                totalSize += restExtraSize;
            }
            if (totalSize <= 0) { reps[gi] = -1; continue; }
            int idx = rng.nextInt(totalSize);
            int pos = (idx < firstSize)
                ? firstLo + idx
                : task.restLos2 + (idx - firstSize);
            reps[gi] = aCons[pos];
        }

        // ── Step (2)+(3): collect induced bipartition counts ────────────────
        java.util.HashMap<Integer, Integer> counts = new java.util.HashMap<>();
        java.util.ArrayList<Integer> perTreeBitmaps = new java.util.ArrayList<>(8);
        for (int ti = 0; ti < geneTrees.size(); ti++) {
            Tree gt = geneTrees.get(ti);
            if (tours != null) collectGeneTreeBitmapsFast(gt, tours[ti], reps, d, perTreeBitmaps);
            else               collectGeneTreeBitmaps(gt, task, reps, d, perTreeBitmaps);
            for (int bm : perTreeBitmaps) {
                Integer cur = counts.get(bm);
                if (cur != null) {
                    counts.put(bm, cur + 1);
                    continue;
                }
                int comp = allBits ^ bm;
                Integer compCur = counts.get(comp);
                if (compCur != null) {
                    counts.put(comp, compCur + 1);
                    continue;
                }
                counts.put(bm, 1);
            }
        }
        if (counts.isEmpty()) return false;

        // ── Step (4): sort by freq desc; deterministic tie-break by bitmap ──
        List<int[]> sorted = new ArrayList<>(counts.size());
        for (var e : counts.entrySet()) sorted.add(new int[]{e.getKey(), e.getValue()});
        sorted.sort((a, b) -> {
            int c = Integer.compare(b[1], a[1]);
            return (c != 0) ? c : Integer.compare(a[0], b[0]);
        });

        // ── Step (5): mini-greedy laminar build with buildTreeFromClusters
        //              semantics (LCA + ≥ 2 children moved).
        MiniGreedyBuilder mg = new MiniGreedyBuilder(d);
        boolean anyAccepted = false;
        boolean productive = false;
        for (int[] entry : sorted) {
            if (mg.tryInsert(entry[0])) {
                anyAccepted = true;
                boolean added = emitInducedSplit(entry[0], task, numTaxa, buffer);
                if (added && isHighSupport(entry[1], geneTrees.size())) productive = true;
            }
        }

        // ── Step (6b) D2: random resolution of leftover multifurcations (opt-in).
        //   ASTRAL-MP resolveLinearly runs this only when the round accepted ≥1 cluster.
        if (Config.getInstance().isStepBRandomLeftoverResolution() && anyAccepted) {
            mg.resolveLeftoverPolytomiesRandomly(rng,
                bm -> emitInducedSplit(bm, task, numTaxa, buffer));
        }

        // ── Step (7): resolveByDistance — UPGMA on the d×d induced similarity
        //              matrix (per-round, on the sampled reps).  Each non-root
        //              internal node of the dendrogram → one emission, mapped
        //              back via group ranges.  This matches ASTRAL-MP's
        //              {@code resolveByDistance} call from sampleAndResolve.
        if (sim != null) {
            stepBResolveByDistance(task, reps, sim, numTaxa, buffer, roundIndex);
        }
        return productive;
    }

    private static boolean isHighSupport(int frequency, int numGeneTrees) {
        return frequency > STEPB_MIN_FREQ
            && (double) frequency / Math.max(1, numGeneTrees) >= STEPB_MIN_RATIO;
    }

    /**
     * resolveByDistance — on the d×d induced similarity matrix over this round's
     * sampled representatives (one taxon per arm; {@code inducedSim[i][j] =
     * sim(rep_i, rep_j)}, the same submatrix ASTRAL-MP's {@code getInducedMatrix}
     * builds).  Two families of emissions, mirroring ASTRAL-MP:
     *   (a) the induced UPGMA tree's bipartitions ({@code inferTreeBitsets}); and
     *   (b) the "quadratic" nearest-neighbour balls ({@code getQuadraticBitsets}) —
     *       gated exactly like ASTRAL-MP: only for the loosest thresholds
     *       ({@code thresholdIndex < STEPB_QUADRATIC_MAX_THRESHOLD_INDEX}) and the
     *       non-bonus rounds ({@code roundIndex < STEPB_DEFAULT_RUNS}).
     */
    private static void stepBResolveByDistance(PolytomyTask task, int[] reps,
                                                 SimilarityMatrix sim,
                                                 int numTaxa, EmissionBuffer buffer,
                                                 int roundIndex) {
        int d = task.numGroups;
        double[] inducedSim = new double[d * d];
        for (int i = 0; i < d; i++) {
            int ri = reps[i];
            if (ri < 0) continue;
            for (int j = i + 1; j < d; j++) {
                int rj = reps[j];
                if (rj < 0) continue;
                double s = sim.getSim(ri, rj);
                inducedSim[i * d + j] = s;
                inducedSim[j * d + i] = s;
            }
        }
        // (a) induced UPGMA tree
        Tree dendro = MiniUPGMA.build(inducedSim, d, /*treeIndex*/0);
        walkDendroAsRepBitmap(dendro.root, dendro.postorderArray, task, numTaxa, buffer);

        // (b) quadratic nearest-neighbour balls (ASTRAL-MP getQuadraticBitsets):
        // opt-in (enlarges X), and gated to the loosest thresholds + non-bonus rounds.
        if (Config.getInstance().isStepBQuadraticNnBalls()
                && task.thresholdIndex < STEPB_QUADRATIC_MAX_THRESHOLD_INDEX
                && roundIndex < STEPB_DEFAULT_RUNS) {
            emitQuadraticBalls(task, reps, inducedSim, d, numTaxa, buffer);
        }
    }

    /**
     * ASTRAL-MP {@code getQuadraticBitsets} on the induced rep matrix: for each
     * arm {@code i}, emit the nested "balls" {i}, {i+nearest}, {i+2 nearest}, …,
     * where neighbours are ordered by descending similarity to {@code i} (ties by
     * arm index — matching ASTRAL-MP's {@code -Float.compare} + index tie-break).
     * Each ball (a subset of arms) is emitted via {@link #emitInducedSplit}, which
     * expands it to the arm union (a multi-range cluster), size-filters the trivial
     * balls, and dedups by signature. O(d² log d) per call; d ≤ 31.
     */
    private static void emitQuadraticBalls(PolytomyTask task, int[] reps,
                                            double[] inducedSim, int d,
                                            int numTaxa, EmissionBuffer buffer) {
        // valid arm indices (skip absent arms)
        Integer[] order = new Integer[d];
        for (int i = 0; i < d; i++) {
            if (reps[i] < 0) continue;
            int rowBase = i * d;
            int cnt = 0;
            for (int j = 0; j < d; j++) if (j != i && reps[j] >= 0) order[cnt++] = j;
            final int fi = i;
            // nearest first: higher induced similarity to i, tie-break smaller index
            java.util.Arrays.sort(order, 0, cnt, (a, b) -> {
                int c = Double.compare(inducedSim[rowBase + b], inducedSim[rowBase + a]);
                return (c != 0) ? c : Integer.compare(a, b);
            });
            int bm = 1 << fi;                       // ball rooted at arm i (i first)
            emitInducedSplit(bm, task, numTaxa, buffer);
            for (int k = 0; k < cnt; k++) {
                bm |= (1 << order[k]);
                emitInducedSplit(bm, task, numTaxa, buffer);
            }
        }
    }

    /** Walk dendrogram; for each non-root internal node, emit the rep-bitmap as a bipartition. */
    private static void walkDendroAsRepBitmap(TreeNode dn, int[] postArr,
                                               PolytomyTask task, int numTaxa,
                                               EmissionBuffer buffer) {
        if (dn.isLeaf()) return;
        walkDendroAsRepBitmap(dn.left,  postArr, task, numTaxa, buffer);
        walkDendroAsRepBitmap(dn.right, postArr, task, numTaxa, buffer);
        if (dn.isRoot()) return;
        int bm = 0;
        for (int p = dn.rangeStart; p < dn.rangeEnd; p++) bm |= (1 << postArr[p]);
        int sz = Integer.bitCount(bm);
        int d = task.numGroups;
        if (sz < 2 || sz > d - 1) return;
        emitInducedSplit(bm, task, numTaxa, buffer);
    }

    /** Postorder walk that fills {@code out} with rep-bitmaps for each qualifying
     *  non-root internal node (mirrors {@code Utils.getBitsets}).  O(n) reference route.
     *  Package-visible for the equivalence harness. */
    static void collectGeneTreeBitmaps(Tree gt, PolytomyTask task, int[] reps,
                                                int d, java.util.ArrayList<Integer> out) {
        out.clear();
        int[] repAtPos = new int[gt.leafCount];
        Arrays.fill(repAtPos, -1);
        int presentCount = 0;
        for (int gi = 0; gi < d; gi++) {
            int r = reps[gi];
            if (r < 0) continue;
            int p = gt.positionMap[r];
            if (p < 0) continue;
            repAtPos[p] = gi;
            presentCount++;
        }
        if (presentCount < 2) return;
        walkCollect(gt.root, /*isRoot=*/true, repAtPos, d, out);
    }

    private static int walkCollect(TreeNode node, boolean isRoot, int[] repAtPos,
                                    int d, java.util.ArrayList<Integer> out) {
        if (node.isLeaf()) {
            int gi = repAtPos[node.rangeStart];
            return (gi >= 0) ? (1 << gi) : 0;
        }
        int bm = 0, legit = 0;
        if (node.isPolytomous()) {              // n-ary gene-tree node: union over all children
            for (TreeNode c : node.children) {
                int cbm = walkCollect(c, false, repAtPos, d, out);
                bm |= cbm;
                if (cbm != 0) legit++;
            }
        } else {
            int leftBM  = walkCollect(node.left,  false, repAtPos, d, out);
            int rightBM = walkCollect(node.right, false, repAtPos, d, out);
            bm = leftBM | rightBM;
            legit = (leftBM != 0 ? 1 : 0) + (rightBM != 0 ? 1 : 0);
        }
        // Skip root (matches the `isRoot` skip in Utils.getBitsets)
        if (isRoot) return bm;
        if (legit < 2) return bm;
        int sz = Integer.bitCount(bm);
        if (sz < 2 || sz >= d - 1) return bm;
        out.add(bm);
        return bm;
    }

    // ── O(d log d) restriction (auxiliary/induced tree via Euler-tour LCA) ──────
    //
    // walkCollect emits a rep-bitmap at exactly the binary MERGE nodes (legit==2),
    // which are precisely the internal nodes of the gene tree restricted to the
    // present reps — each once. So we can enumerate those clades directly from the
    // induced tree without touching all n leaves:
    //   1. Locate the ≤ d present reps and sort by leaf position (DFS order). O(d log d).
    //   2. separator depth between consecutive reps = depth(LCA) via O(1) Euler RMQ.
    //   3. The induced tree is the Cartesian tree on those separator depths; each
    //      internal node spans a contiguous rep interval. Enumerate via a min-split
    //      recursion (O(d²) on the ≤31-rep sequence — trivially fast), applying the
    //      SAME size filter (2 ≤ sz ≤ d-2) and gene-tree-root skip (depth 0) as
    //      walkCollect. Produces the identical per-tree bitmap set.

    /** Package-visible for the equivalence harness. */
    static void collectGeneTreeBitmapsFast(Tree gt, EulerTourBuilder.TourData tour,
                                           int[] reps, int d, java.util.ArrayList<Integer> out) {
        out.clear();
        int[] pos = new int[d], tax = new int[d], grp = new int[d];
        int k = 0;
        for (int gi = 0; gi < d; gi++) {
            int r = reps[gi];
            if (r < 0) continue;
            int p = gt.positionMap[r];
            if (p < 0) continue;
            pos[k] = p; tax[k] = r; grp[k] = gi; k++;
        }
        if (k < 2) return;
        // insertion sort by leaf position (k ≤ 31)
        for (int i = 1; i < k; i++) {
            int pp = pos[i], tt = tax[i], gg = grp[i], j = i - 1;
            while (j >= 0 && pos[j] > pp) { pos[j+1]=pos[j]; tax[j+1]=tax[j]; grp[j+1]=grp[j]; j--; }
            pos[j+1]=pp; tax[j+1]=tt; grp[j+1]=gg;
        }
        int[] sep = new int[k - 1];
        for (int i = 0; i < k - 1; i++) sep[i] = lcaDepth(tour, tax[i], tax[i+1]);
        int[] bit = new int[k];
        for (int i = 0; i < k; i++) bit[i] = 1 << grp[i];
        emitInducedClades(0, k - 1, sep, bit, d, out);
    }

    /** Min-split recursion over the sorted rep interval [lo,hi]; returns its bitmap. */
    private static int emitInducedClades(int lo, int hi, int[] sep, int[] bit, int d,
                                         java.util.ArrayList<Integer> out) {
        if (lo == hi) return bit[lo];
        // The LCA of reps[lo..hi] is the unique minimal-depth separator (binary tree).
        int m = lo, minDepth = sep[lo];
        for (int i = lo + 1; i <= hi - 1; i++) if (sep[i] < minDepth) { minDepth = sep[i]; m = i; }
        int leftBM  = emitInducedClades(lo, m, sep, bit, d, out);
        int rightBM = emitInducedClades(m + 1, hi, sep, bit, d, out);
        int bm = leftBM | rightBM;
        int sz = Integer.bitCount(bm);
        // depth 0 ⇒ this merge is the gene-tree root (walkCollect's isRoot skip).
        if (minDepth != 0 && sz >= 2 && sz <= d - 2) out.add(bm);
        return bm;
    }

    /** depth(LCA(taxonA, taxonB)) via O(1) sparse-table RMQ over the Euler tour.
     *  Queries are leaf-pairs, so the range length stays < tourLen ≤ 1<<log. */
    private static int lcaDepth(EulerTourBuilder.TourData tour, int taxonA, int taxonB) {
        int fa = tour.firstOcc[taxonA], fb = tour.firstOcc[taxonB];
        int lo = Math.min(fa, fb), hi = Math.max(fa, fb);
        int len = hi - lo + 1;
        int k = 31 - Integer.numberOfLeadingZeros(len);     // floor(log2(len))
        int d1 = tour.sparseMin[k][lo];
        int d2 = tour.sparseMin[k][hi - (1 << k) + 1];
        return Math.min(d1, d2);
    }

    /** Convert a rep-bitmap into a group-bipartition and emit (smaller side). */
    private static boolean emitInducedSplit(int repBitmap, PolytomyTask task,
                                             int numTaxa, EmissionBuffer buffer) {
        int d = task.numGroups;
        boolean[] selected = new boolean[d];
        int selectedSize = 0;
        for (int gi = 0; gi < d; gi++) {
            if ((repBitmap & (1 << gi)) != 0) {
                selected[gi] = true;
                selectedSize += groupSize(task, gi);
            }
        }
        int complementSize = numTaxa - selectedSize;
        if (selectedSize <= 1 || complementSize <= 1) return false;
        if (selectedSize == numTaxa) return false;

        MultiRange canonical;
        int canonicalSize;
        if (selectedSize <= complementSize) {
            canonical = buildSideMultiRange(task, selected, /*inverted=*/false);
            canonicalSize = selectedSize;
        } else {
            canonical = buildSideMultiRange(task, selected, /*inverted=*/true);
            canonicalSize = complementSize;
        }

        int m = task.tree.numSeeds();
        long[] sums = new long[m];
        long[] xors = new long[m];
        for (int s = 0; s < m; s++) {
            sums[s] = task.tree.combineDisjointSigma1(s, canonical.los, canonical.his);
            xors[s] = task.tree.combineDisjointSigma2(s, canonical.los, canonical.his);
        }
        ClusterHash sig = new ClusterHash(sums, xors, canonicalSize, m);
        return buffer.add(new EmittedBipartition(
            sig, canonical, canonicalSize, 'B', task.thresholdIndex));
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  LARGE-POLYTOMY (d > 31) long[]-bitmap Step B path
    //  ----------------------------------------------------------------------
    //  Mirrors the int path above method-for-method, with rep-subsets stored as
    //  long[W] (W = ceil(d/64)) instead of a single int.  Reached only when
    //  --stepb-process-large-polytomies is set (so the int path / existing
    //  behaviour for d ≤ 31 is untouched).  Quadratic NN-balls are intentionally
    //  NOT run here — matching ASTRAL-MP, which disables the quadratic family for
    //  polytomies above its size limit while still running the linear path.
    // ═══════════════════════════════════════════════════════════════════════

    /** Hashable / orderable wrapper over a long[] rep-subset bitmap. */
    static final class BitKey {
        final long[] w;
        private final int h;
        BitKey(long[] w) {
            this.w = w;
            int hh = 1;
            for (long x : w) hh = 31 * hh + (int) (x ^ (x >>> 32));
            this.h = hh;
        }
        @Override public int hashCode() { return h; }
        @Override public boolean equals(Object o) {
            return (o instanceof BitKey) && java.util.Arrays.equals(w, ((BitKey) o).w);
        }
        /** Ascending unsigned-magnitude order (deterministic tie-break). */
        static int compare(BitKey a, BitKey b) {
            for (int k = a.w.length - 1; k >= 0; k--) {
                int c = Long.compareUnsigned(a.w[k], b.w[k]);
                if (c != 0) return c;
            }
            return 0;
        }
    }

    // ── Memory-lean hash-based counting for the large (d>31) fast path ──────────
    //
    // Instead of storing a long[⌈d/64⌉] bitmap per distinct induced clade (the
    // dominant RAM term for big polytomies), we identify each clade by a 128-bit
    // (sum,xor) hash of its arm-group set — the same hashing principle STELAR-X uses
    // for every cluster (ClusterHash); dedup reliability is the identical 2⁻¹²⁸ level.
    // We keep only frequency + a compact provenance (treeIdx, [lo,hi] interval of the
    // position-sorted reps that produced it) and a small per-tree grp[] cache, then
    // materialize the bitmap just-in-time, one at a time, for the mini-greedy.
    // (Tie-break among equal-frequency candidates becomes hash-order instead of
    // bitmap-magnitude — a different but equal-quality emission set, within the same
    // run-to-run nondeterminism Step B already has.  d≤31 int path is untouched.)

    /** 128-bit set signature (sum,xor of per-group hashes) — the counts-map key. */
    static final class CladeKey {
        final long s, x;
        CladeKey(long s, long x) { this.s = s; this.x = x; }
        @Override public int hashCode() {
            long h = s * 0x9E3779B97F4A7C15L ^ x; return (int) (h ^ (h >>> 32));
        }
        @Override public boolean equals(Object o) {
            return (o instanceof CladeKey k) && k.s == s && k.x == x;
        }
        static int compare(CladeKey a, CladeKey b) {
            int c = Long.compareUnsigned(a.s, b.s);
            return (c != 0) ? c : Long.compareUnsigned(a.x, b.x);
        }
    }
    /** Frequency + provenance to re-materialize the clade's rep bitmap on demand. */
    static final class CladeEntry {
        int freq; final int ti, lo, hi;
        CladeEntry(int ti, int lo, int hi) { this.ti = ti; this.lo = lo; this.hi = hi; this.freq = 1; }
    }

    /** SplitMix64 finalizer — well-distributed per-group hash. */
    private static long mix64(long z) {
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    /**
     * Fast-path clade enumeration with HASH counting (no per-clade long[] bitmaps).
     * Sorts the present reps, runs the SAME min-split induced-clade recursion as
     * {@link #collectGeneTreeBitmapsFastLong}, but records each qualifying clade by a
     * 128-bit (sum,xor) hash + (ti, interval) provenance into {@code counts}.
     * @return the sorted grp[] (provenance source for this tree), or null if &lt;2 reps.
     */
    private static int[] collectCladesHash(Tree gt, EulerTourBuilder.TourData tour,
            int[] reps, int d, long[] vg, int ti, long allSum, long allXor,
            java.util.HashMap<CladeKey, CladeEntry> counts) {
        int[] pos = new int[d], tax = new int[d], grp = new int[d];
        int k = 0;
        for (int gi = 0; gi < d; gi++) {
            int r = reps[gi];
            if (r < 0) continue;
            int p = gt.positionMap[r];
            if (p < 0) continue;
            pos[k] = p; tax[k] = r; grp[k] = gi; k++;
        }
        if (k < 2) return null;
        for (int i = 1; i < k; i++) {            // insertion sort by leaf position
            int pp = pos[i], tt = tax[i], gg = grp[i], j = i - 1;
            while (j >= 0 && pos[j] > pp) { pos[j+1]=pos[j]; tax[j+1]=tax[j]; grp[j+1]=grp[j]; j--; }
            pos[j+1]=pp; tax[j+1]=tt; grp[j+1]=gg;
        }
        int[] sep = new int[k - 1];
        for (int i = 0; i < k - 1; i++) sep[i] = lcaDepth(tour, tax[i], tax[i+1]);
        int[] grpArr = java.util.Arrays.copyOf(grp, k);
        long[] psum = new long[k + 1], pxor = new long[k + 1];
        for (int i = 0; i < k; i++) { long h = vg[grpArr[i]]; psum[i+1] = psum[i] + h; pxor[i+1] = pxor[i] ^ h; }
        emitCladesHash(0, k - 1, sep, psum, pxor, d, ti, allSum, allXor, counts);
        return grpArr;
    }

    private static void emitCladesHash(int lo, int hi, int[] sep, long[] psum, long[] pxor,
            int d, int ti, long allSum, long allXor, java.util.HashMap<CladeKey, CladeEntry> counts) {
        if (lo == hi) return;                    // singleton — not an emitted clade
        int m = lo, minDepth = sep[lo];
        for (int i = lo + 1; i <= hi - 1; i++) if (sep[i] < minDepth) { minDepth = sep[i]; m = i; }
        emitCladesHash(lo, m, sep, psum, pxor, d, ti, allSum, allXor, counts);
        emitCladesHash(m + 1, hi, sep, psum, pxor, d, ti, allSum, allXor, counts);
        int sz = hi - lo + 1;                    // reps are distinct groups ⇒ popcount = interval length
        if (minDepth != 0 && sz >= 2 && sz <= d - 2) {
            long sum = psum[hi + 1] - psum[lo];
            long xor = pxor[hi + 1] ^ pxor[lo];
            CladeKey key = new CladeKey(sum, xor);
            CladeEntry e = counts.get(key);
            if (e != null) { e.freq++; return; }
            CladeEntry ce = counts.get(new CladeKey(allSum - sum, allXor ^ xor));  // complement
            if (ce != null) { ce.freq++; return; }
            counts.put(key, new CladeEntry(ti, lo, hi));
        }
    }

    /** Re-materialize a clade's rep bitmap into {@code buf} from its provenance. */
    private static void materializeClade(long[] buf, int W, int[] grp, int lo, int hi) {
        java.util.Arrays.fill(buf, 0, W, 0L);
        for (int p = lo; p <= hi; p++) { int g = grp[p]; buf[g >>> 6] |= (1L << (g & 63)); }
    }

    /** Shared large-polytomy round tail: D2 leftover and resolveByDistance. */
    private static void finishRoundLong(MiniGreedyBuilderLong mg, boolean anyAccepted,
            PolytomyTask task, int numTaxa, EmissionBuffer buffer, Random rng,
            SimilarityMatrix sim, int[] reps, int W, double[] simBuf) {
        if (Config.getInstance().isStepBRandomLeftoverResolution() && anyAccepted) {
            mg.resolveLeftoverPolytomiesRandomly(rng,
                bm -> emitInducedSplitLong(bm, task, numTaxa, buffer));
        }
        if (sim != null) stepBResolveByDistanceLong(task, reps, sim, numTaxa, buffer, W, simBuf);
    }

    private static boolean stepBRoundLong(PolytomyTask task, List<Tree> geneTrees,
                                          EulerTourBuilder.TourData[] tours,
                                          EmissionBuffer buffer, int numTaxa, Random rng,
                                          SimilarityMatrix sim, int roundIndex, double[] simBuf) {
        int d = task.numGroups;
        int W = (d + 63) >>> 6;
        int[] aCons = task.tree.aCons();

        // ── (1) one random rep per group (identical logic to the int path) ──
        int[] reps = new int[d];
        for (int gi = 0; gi < d; gi++) {
            int firstLo = task.groupLos[gi], firstHi = task.groupHis[gi];
            int firstSize = firstHi - firstLo;
            int totalSize = firstSize;
            boolean isRest = (gi == d - 1) && (task.numGroups > task.degree);
            if (isRest && task.restSplit) totalSize += task.restHis2 - task.restLos2;
            if (totalSize <= 0) { reps[gi] = -1; continue; }
            int idx = rng.nextInt(totalSize);
            int pos = (idx < firstSize) ? firstLo + idx
                                        : task.restLos2 + (idx - firstSize);
            reps[gi] = aCons[pos];
        }

        // ── (2)–(5) collect induced bipartitions, sort by freq, mini-greedy build ──
        MiniGreedyBuilderLong mg = new MiniGreedyBuilderLong(d);
        boolean anyAccepted = false;
        boolean productive = false;

        if (tours != null) {
            // FAST path: 128-bit hash counting + JIT bitmap materialization (memory-lean).
            long[] vg = new long[d];
            long allSum = 0, allXor = 0;
            for (int g = 0; g < d; g++) {
                long h = mix64((g + 1) * 0x9E3779B97F4A7C15L);
                vg[g] = h; allSum += h; allXor ^= h;
            }
            int nT = geneTrees.size();
            int[][] grpCache = new int[nT][];
            java.util.HashMap<CladeKey, CladeEntry> counts = new java.util.HashMap<>();
            for (int ti = 0; ti < nT; ti++) {
                grpCache[ti] = collectCladesHash(geneTrees.get(ti), tours[ti], reps, d, vg,
                                                 ti, allSum, allXor, counts);
            }
            if (counts.isEmpty()) return false;
            List<java.util.Map.Entry<CladeKey, CladeEntry>> sorted =
                new ArrayList<>(counts.entrySet());
            sorted.sort((a, b) -> {
                int c = Integer.compare(b.getValue().freq, a.getValue().freq);
                return (c != 0) ? c : CladeKey.compare(a.getKey(), b.getKey());
            });
            long[] buf = new long[W];
            for (var e : sorted) {
                CladeEntry ce = e.getValue();
                materializeClade(buf, W, grpCache[ce.ti], ce.lo, ce.hi);
                if (mg.tryInsert(buf)) {
                    anyAccepted = true;
                    boolean added = emitInducedSplitLong(buf, task, numTaxa, buffer);
                    if (added && isHighSupport(ce.freq, geneTrees.size())) productive = true;
                }
            }
        } else {
            // REFERENCE (non-fast) path: exact long[] bitmap counting.
            long[] allBits = new long[W];
            for (int b = 0; b < d; b++) allBits[b >>> 6] |= (1L << (b & 63));
            java.util.HashMap<BitKey, Integer> counts = new java.util.HashMap<>();
            java.util.ArrayList<long[]> perTree = new java.util.ArrayList<>(8);
            for (int ti = 0; ti < geneTrees.size(); ti++) {
                collectGeneTreeBitmapsLong(geneTrees.get(ti), reps, d, W, perTree);
                for (long[] bm : perTree) {
                    BitKey key = new BitKey(bm);
                    Integer cur = counts.get(key);
                    if (cur != null) { counts.put(key, cur + 1); continue; }
                    long[] comp = new long[W];
                    for (int k = 0; k < W; k++) comp[k] = allBits[k] ^ bm[k];
                    BitKey ckey = new BitKey(comp);
                    Integer cc = counts.get(ckey);
                    if (cc != null) { counts.put(ckey, cc + 1); continue; }
                    counts.put(key, 1);
                }
            }
            if (counts.isEmpty()) return false;
            List<java.util.Map.Entry<BitKey, Integer>> sorted =
                new ArrayList<>(counts.entrySet());
            sorted.sort((a, b) -> {
                int c = Integer.compare(b.getValue(), a.getValue());
                return (c != 0) ? c : BitKey.compare(a.getKey(), b.getKey());
            });
            for (var e : sorted) {
                if (mg.tryInsert(e.getKey().w)) {
                    anyAccepted = true;
                    boolean added = emitInducedSplitLong(
                        e.getKey().w, task, numTaxa, buffer);
                    if (added && isHighSupport(e.getValue(), geneTrees.size())) productive = true;
                }
            }
        }

        // ── (6)+(7) emit accepted, D2 leftover, resolveByDistance (shared tail) ──
        finishRoundLong(mg, anyAccepted, task, numTaxa, buffer, rng, sim, reps, W, simBuf);
        return productive;
    }

    private static void stepBResolveByDistanceLong(PolytomyTask task, int[] reps,
                                                   SimilarityMatrix sim, int numTaxa,
                                                   EmissionBuffer buffer, int W, double[] simBuf) {
        int d = task.numGroups;
        // Reuse the per-task buffer (zero-filled here; only i<j and j>i entries are
        // written below, so stale values from a prior round must be cleared first).
        double[] inducedSim = (simBuf != null) ? simBuf : new double[d * d];
        java.util.Arrays.fill(inducedSim, 0, d * d, 0.0);
        for (int i = 0; i < d; i++) {
            int ri = reps[i];
            if (ri < 0) continue;
            for (int j = i + 1; j < d; j++) {
                int rj = reps[j];
                if (rj < 0) continue;
                double s = sim.getSim(ri, rj);
                inducedSim[i * d + j] = s;
                inducedSim[j * d + i] = s;
            }
        }
        Tree dendro = MiniUPGMA.buildFast(inducedSim, d, /*treeIndex*/0);
        walkDendroAsRepBitmapLong(dendro.root, dendro.postorderArray, task, numTaxa, buffer, d, W);
    }

    private static void walkDendroAsRepBitmapLong(TreeNode dn, int[] postArr,
                                                  PolytomyTask task, int numTaxa,
                                                  EmissionBuffer buffer, int d, int W) {
        if (dn.isLeaf()) return;
        walkDendroAsRepBitmapLong(dn.left,  postArr, task, numTaxa, buffer, d, W);
        walkDendroAsRepBitmapLong(dn.right, postArr, task, numTaxa, buffer, d, W);
        if (dn.isRoot()) return;
        long[] bm = new long[W];
        for (int p = dn.rangeStart; p < dn.rangeEnd; p++) {
            int g = postArr[p];
            bm[g >>> 6] |= (1L << (g & 63));
        }
        int sz = popcountL(bm);
        if (sz < 2 || sz > d - 1) return;
        emitInducedSplitLong(bm, task, numTaxa, buffer);
    }

    /** O(d log d) Euler-tour induced-clade enumeration (long[] variant of
     *  {@link #collectGeneTreeBitmapsFast}). */
    static void collectGeneTreeBitmapsFastLong(Tree gt, EulerTourBuilder.TourData tour,
                                               int[] reps, int d, int W,
                                               java.util.ArrayList<long[]> out) {
        out.clear();
        int[] pos = new int[d], tax = new int[d], grp = new int[d];
        int k = 0;
        for (int gi = 0; gi < d; gi++) {
            int r = reps[gi];
            if (r < 0) continue;
            int p = gt.positionMap[r];
            if (p < 0) continue;
            pos[k] = p; tax[k] = r; grp[k] = gi; k++;
        }
        if (k < 2) return;
        for (int i = 1; i < k; i++) {       // insertion sort by leaf position
            int pp = pos[i], tt = tax[i], gg = grp[i], j = i - 1;
            while (j >= 0 && pos[j] > pp) { pos[j+1]=pos[j]; tax[j+1]=tax[j]; grp[j+1]=grp[j]; j--; }
            pos[j+1]=pp; tax[j+1]=tt; grp[j+1]=gg;
        }
        int[] sep = new int[k - 1];
        for (int i = 0; i < k - 1; i++) sep[i] = lcaDepth(tour, tax[i], tax[i+1]);
        int[] grpArr = java.util.Arrays.copyOf(grp, k);
        emitInducedCladesLong(0, k - 1, sep, grpArr, d, W, out);
    }

    private static long[] emitInducedCladesLong(int lo, int hi, int[] sep, int[] grp,
                                                int d, int W, java.util.ArrayList<long[]> out) {
        if (lo == hi) {
            long[] bm = new long[W];
            bm[grp[lo] >>> 6] |= (1L << (grp[lo] & 63));
            return bm;
        }
        int m = lo, minDepth = sep[lo];
        for (int i = lo + 1; i <= hi - 1; i++) if (sep[i] < minDepth) { minDepth = sep[i]; m = i; }
        long[] leftBM  = emitInducedCladesLong(lo, m, sep, grp, d, W, out);
        long[] rightBM = emitInducedCladesLong(m + 1, hi, sep, grp, d, W, out);
        long[] bm = new long[W];
        for (int k = 0; k < W; k++) bm[k] = leftBM[k] | rightBM[k];
        int sz = popcountL(bm);
        if (minDepth != 0 && sz >= 2 && sz <= d - 2) out.add(bm);
        return bm;
    }

    /** O(n) reference walk (long[] variant of {@link #collectGeneTreeBitmaps}). */
    static void collectGeneTreeBitmapsLong(Tree gt, int[] reps, int d, int W,
                                           java.util.ArrayList<long[]> out) {
        out.clear();
        int[] repAtPos = new int[gt.leafCount];
        Arrays.fill(repAtPos, -1);
        int presentCount = 0;
        for (int gi = 0; gi < d; gi++) {
            int r = reps[gi];
            if (r < 0) continue;
            int p = gt.positionMap[r];
            if (p < 0) continue;
            repAtPos[p] = gi;
            presentCount++;
        }
        if (presentCount < 2) return;
        walkCollectLong(gt.root, /*isRoot=*/true, repAtPos, d, W, out);
    }

    private static long[] walkCollectLong(TreeNode node, boolean isRoot, int[] repAtPos,
                                          int d, int W, java.util.ArrayList<long[]> out) {
        if (node.isLeaf()) {
            int gi = repAtPos[node.rangeStart];
            long[] bm = new long[W];
            if (gi >= 0) bm[gi >>> 6] |= (1L << (gi & 63));
            return bm;
        }
        long[] bm = new long[W];
        int legit = 0;
        if (node.isPolytomous()) {               // n-ary gene-tree node
            for (TreeNode c : node.children) {
                long[] cbm = walkCollectLong(c, false, repAtPos, d, W, out);
                for (int k = 0; k < W; k++) bm[k] |= cbm[k];
                if (nonZeroL(cbm)) legit++;
            }
        } else {
            long[] leftBM  = walkCollectLong(node.left,  false, repAtPos, d, W, out);
            long[] rightBM = walkCollectLong(node.right, false, repAtPos, d, W, out);
            for (int k = 0; k < W; k++) bm[k] = leftBM[k] | rightBM[k];
            legit = (nonZeroL(leftBM) ? 1 : 0) + (nonZeroL(rightBM) ? 1 : 0);
        }
        if (isRoot) return bm;
        if (legit < 2) return bm;
        int sz = popcountL(bm);
        if (sz < 2 || sz >= d - 1) return bm;
        out.add(bm);
        return bm;
    }

    /** Convert a long[] rep-bitmap into a group bipartition and emit (smaller side). */
    private static boolean emitInducedSplitLong(long[] repBitmap, PolytomyTask task,
                                                int numTaxa, EmissionBuffer buffer) {
        int d = task.numGroups;
        boolean[] selected = new boolean[d];
        int selectedSize = 0;
        for (int gi = 0; gi < d; gi++) {
            if ((repBitmap[gi >>> 6] & (1L << (gi & 63))) != 0) {
                selected[gi] = true;
                selectedSize += groupSize(task, gi);
            }
        }
        int complementSize = numTaxa - selectedSize;
        if (selectedSize <= 1 || complementSize <= 1) return false;
        if (selectedSize == numTaxa) return false;

        MultiRange canonical;
        int canonicalSize;
        if (selectedSize <= complementSize) {
            canonical = buildSideMultiRange(task, selected, /*inverted=*/false);
            canonicalSize = selectedSize;
        } else {
            canonical = buildSideMultiRange(task, selected, /*inverted=*/true);
            canonicalSize = complementSize;
        }

        int m = task.tree.numSeeds();
        long[] sums = new long[m];
        long[] xors = new long[m];
        for (int s = 0; s < m; s++) {
            sums[s] = task.tree.combineDisjointSigma1(s, canonical.los, canonical.his);
            xors[s] = task.tree.combineDisjointSigma2(s, canonical.los, canonical.his);
        }
        ClusterHash sig = new ClusterHash(sums, xors, canonicalSize, m);
        return buffer.add(new EmittedBipartition(
            sig, canonical, canonicalSize, 'B', task.thresholdIndex));
    }

    private static int popcountL(long[] bm) {
        int c = 0;
        for (long x : bm) c += Long.bitCount(x);
        return c;
    }
    private static boolean nonZeroL(long[] bm) {
        for (long x : bm) if (x != 0) return true;
        return false;
    }

    // ─────────────────────────────────────────────────────────────────────

    /** Build a MultiRange of the chosen groups (or their complement when inverted). */
    private static MultiRange buildSideMultiRange(PolytomyTask task,
                                                   boolean[] selected,
                                                   boolean inverted) {
        // Collect (lo, hi) pairs, splitting rest into two when applicable.
        List<int[]> ranges = new ArrayList<>();
        for (int i = 0; i < task.numGroups; i++) {
            boolean want = inverted ? !selected[i] : selected[i];
            if (!want) continue;
            int lo = task.groupLos[i], hi = task.groupHis[i];
            if (lo < hi) ranges.add(new int[]{lo, hi});
            boolean isRest = (i == task.numGroups - 1) && (task.numGroups > task.degree);
            if (isRest && task.restSplit) {
                ranges.add(new int[]{task.restLos2, task.restHis2});
            }
        }
        int[] los = new int[ranges.size()];
        int[] his = new int[ranges.size()];
        for (int k = 0; k < ranges.size(); k++) {
            los[k] = ranges.get(k)[0];
            his[k] = ranges.get(k)[1];
        }
        return new MultiRange(task.tree, los, his);
    }
}
