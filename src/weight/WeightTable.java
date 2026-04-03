package weight;

import utils.Logging;
import cluster.Cluster;
import cluster.ClusterHash;
import cluster.ClusterTable;
import dp.BipartitionSplit;
import dp.DPTable;
import partition.RangePartition;
import partition.PartitionTable;
import tree.Tree;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Precomputed STELAR-X triplet scores for every candidate bipartition split.
 *
 * For a candidate split (A | B) inducing species-tree tripartition (A | B | C)
 * where C = S \ (A u B), and a gene-tree tripartition (M1 | M2 | M3):
 *
 *   Compute 4 core intersections: a0=|A n M1|, a1=|A n M2|, b0=|B n M1|, b1=|B n M2|
 *
 *   Derive remaining 5 by row/column constraints:
 *     lgA = |A n Lg_GT|   (= sizeA for complete trees)
 *     lgB = |B n Lg_GT|   (= sizeB for complete trees)
 *     a2 = lgA - a0 - a1
 *     b2 = lgB - b0 - b1
 *     c0 = sz1 - a0 - b0
 *     c1 = sz2 - a1 - b1
 *     c2 = sz3 - a2 - b2
 *
 *   twoTripletScore = sum over 6 permutations (i,j,k) of {0,1,2}:
 *     a[i] * b[j] * c[k] * (a[i] + b[j] + c[k] - 3)
 *     (only terms where a[i]+b[j]+c[k]-3 > 0 contribute)
 *
 *   score(A|B) = sum over unique gene-tree tripartitions P:
 *                  P.frequency * twoTripletScore / 2
 *
 * Execution: multi-threaded CPU using Threading.processRangeParallel pattern.
 */
public class WeightTable {

    private final Map<BipartitionSplit, Long> scores = new HashMap<>();
    private final int n;   // total taxa

    private long maxScore;
    private long totalScore;

    // -------------------------------------------------------------------------

    public WeightTable(DPTable dpTable, PartitionTable partTable,
                       ClusterTable clusterTable, List<Tree> trees) {
        long t0 = System.nanoTime();
        this.n = clusterTable.getAllTaxaHash().size;

        // Collect all unique splits from DPTable into an indexed list
        List<BipartitionSplit> splitList = new ArrayList<>();
        for (Map.Entry<ClusterHash, Set<BipartitionSplit>> entry : dpTable.entries()) {
            splitList.addAll(entry.getValue());
        }
        int numSplits = splitList.size();

        Collection<PartitionTable.Entry> partitions = partTable.entries();
        long[] scoreArray = new long[numSplits];

        // Multi-threaded scoring
        int numThreads = Runtime.getRuntime().availableProcessors();
        utils.Threading.startThreading(numThreads);

        int chunkSize = Math.max(1, (numSplits + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (numSplits + chunkSize - 1) / chunkSize);
        java.util.concurrent.CountDownLatch latch = new java.util.concurrent.CountDownLatch(actualThreads);
        AtomicInteger done = new AtomicInteger(0);

        for (int t = 0; t < actualThreads; t++) {
            final int startIdx = t * chunkSize;
            final int endIdx   = Math.min(startIdx + chunkSize, numSplits);
            utils.Threading.execute(() -> {
                try {
                    for (int idx = startIdx; idx < endIdx; idx++) {
                        scoreArray[idx] = computeScore(splitList.get(idx), partitions, clusterTable, trees);
                    }
                } finally {
                    latch.countDown();
                }
            });
        }

        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Weight calculation interrupted", e);
        } finally {
            utils.Threading.shutdown();
        }

        for (int i = 0; i < numSplits; i++) {
            scores.put(splitList.get(i), scoreArray[i]);
            if (scoreArray[i] > maxScore) maxScore = scoreArray[i];
            totalScore += scoreArray[i];
        }

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Weight table: %d splits scored, maxScore=%d, totalScore=%d in %d ms",
            scores.size(), maxScore, totalScore, ms);
    }

    // -------------------------------------------------------------------------
    // CPU scoring
    // -------------------------------------------------------------------------

    private long computeScore(BipartitionSplit split,
                               Collection<PartitionTable.Entry> partitions,
                               ClusterTable clusterTable, List<Tree> trees) {
        ClusterTable.Entry eA = clusterTable.get(split.lo);
        ClusterTable.Entry eB = clusterTable.get(split.hi);
        if (eA == null || eB == null) return 0L;

        Cluster cA = eA.exemplar;
        Cluster cB = eB.exemplar;
        Tree tA = trees.get(cA.treeIndex);
        Tree tB = trees.get(cB.treeIndex);
        int sizeA = cA.size;
        int sizeB = cB.size;
        int sizeC = n - sizeA - sizeB;
        if (sizeC < 0) return 0L;

        long twoScore = 0L;

        for (PartitionTable.Entry pe : partitions) {
            RangePartition p = pe.exemplar;
            Tree tGT = trees.get(p.treeIndex);

            int lo1 = p.leftStart,  hi1 = p.leftEnd;   // M1 range
            int lo2 = p.rightStart, hi2 = p.rightEnd;  // M2 range
            int sz1 = p.size1, sz2 = p.size2, sz3 = p.size3;

            // 4 core intersections
            int a0 = IntersectionCounter.intersect(tGT, lo1, hi1, tA, cA.left, cA.right, cA.complement, sz1);
            int a1 = IntersectionCounter.intersect(tGT, lo2, hi2, tA, cA.left, cA.right, cA.complement, sz2);
            int b0 = IntersectionCounter.intersect(tGT, lo1, hi1, tB, cB.left, cB.right, cB.complement, sz1);
            int b1 = IntersectionCounter.intersect(tGT, lo2, hi2, tB, cB.left, cB.right, cB.complement, sz2);

            // Row sums: for incomplete gene trees |A n Lg_GT| < sizeA; must compute explicitly
            int lgA = tGT.isComplete ? sizeA
                    : IntersectionCounter.intersectWithFullTree(tGT, tA, cA.left, cA.right, cA.complement);
            int lgB = tGT.isComplete ? sizeB
                    : IntersectionCounter.intersectWithFullTree(tGT, tB, cB.left, cB.right, cB.complement);

            // Derive remaining 5
            int a2 = lgA - a0 - a1;
            int b2 = lgB - b0 - b1;
            int c0 = sz1 - a0 - b0;
            int c1 = sz2 - a1 - b1;
            int c2 = sz3 - a2 - b2;

            if (a2 < 0 || b2 < 0 || c0 < 0 || c1 < 0 || c2 < 0) continue;

            long twoTriplet = computeTwoTripletScore(a0, a1, a2, b0, b1, b2, c0, c1, c2);
            twoScore += (long) pe.frequency * twoTriplet;
        }

        return twoScore / 2;
    }

    /**
     * Compute 2 * tripletScore = sum over 6 permutations (i,j,k) of {0,1,2}:
     *   a[i] * b[j] * c[k] * (a[i] + b[j] + c[k] - 3)
     *
     * Only terms where (a[i] + b[j] + c[k] - 3) > 0 contribute.
     * The result is always a non-negative even integer.
     */
    private static long computeTwoTripletScore(int a0, int a1, int a2,
                                                int b0, int b1, int b2,
                                                int c0, int c1, int c2) {
        long[] a = {a0, a1, a2};
        long[] b = {b0, b1, b2};
        long[] c = {c0, c1, c2};

        long sum = 0;
        int[][] perms = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
        for (int[] perm : perms) {
            long ai = a[perm[0]], bj = b[perm[1]], ck = c[perm[2]];
            long s = ai + bj + ck - 3;
            if (s > 0) sum += ai * bj * ck * s;
        }
        return sum;
    }

    // -------------------------------------------------------------------------
    // Queries
    // -------------------------------------------------------------------------

    public long getScore(BipartitionSplit split) {
        return scores.getOrDefault(split, 0L);
    }

    public long getMaxScore()   { return maxScore; }
    public long getTotalScore() { return totalScore; }
    public int  size()          { return scores.size(); }

    public Set<Map.Entry<BipartitionSplit, Long>> entries() { return scores.entrySet(); }
}
