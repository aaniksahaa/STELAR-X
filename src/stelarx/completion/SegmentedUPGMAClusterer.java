package stelarx.completion;

import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.Threading;

import java.util.Arrays;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;

/**
 * Exact large-matrix counterpart of {@link UPGMAClusterer}.
 *
 * The algorithm, active-set ordering, comparisons, and tie behavior intentionally
 * mirror the dense implementation.  Only the mutable upper-triangle backing
 * store and its index arithmetic differ.
 */
final class SegmentedUPGMAClusterer {
    private final int n;
    private final SegmentedDoubleArray mat;
    private final double[] bestSim;
    private final int[] bestJ;
    private final int[] weight;
    private final TreeNode[] clusterRoot;
    private final int[] activeArr;
    private final int[] posInActive;
    private int activeCount;

    private SegmentedUPGMAClusterer(SimilarityMatrix sim) {
        this.n = sim.n;
        this.mat = sim.copyPackedSimilarity();
        this.bestSim = new double[n];
        this.bestJ = new int[n];
        this.weight = new int[n];
        this.clusterRoot = new TreeNode[n];
        this.activeArr = new int[n];
        this.posInActive = new int[n];
        this.activeCount = n;

        for (int i = 0; i < n; i++) {
            activeArr[i] = i;
            posInActive[i] = i;
            weight[i] = 1;
            TreeNode leaf = new TreeNode();
            leaf.taxonId = i;
            clusterRoot[i] = leaf;
        }

        // Same i-then-j scan and strict comparison as the dense implementation.
        for (int i = 0; i < n; i++) {
            double bs = Double.NEGATIVE_INFINITY;
            int bj = -1;
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                double s = get(i, j);
                if (s > bs) { bs = s; bj = j; }
            }
            bestSim[i] = bs;
            bestJ[i] = bj;
        }
    }

    static Tree build(SimilarityMatrix sim, int treeIndex) {
        return new SegmentedUPGMAClusterer(sim).run(treeIndex);
    }

    private Tree run(int treeIndex) {
        int threads = Threading.getNumThreads();
        for (int iter = 0; iter < n - 1; iter++) {
            int I = findBestI(threads);
            int J = bestJ[I];

            TreeNode newNode = new TreeNode();
            newNode.left = clusterRoot[I];
            newNode.right = clusterRoot[J];
            clusterRoot[I].parent = newNode;
            clusterRoot[J].parent = newNode;
            clusterRoot[I] = newNode;

            final double wI = weight[I];
            final double wJ = weight[J];
            weight[I] += weight[J];
            removeActive(J);

            final int snapCount = activeCount;
            final int finalI = I;
            final int finalJ = J;
            ConcurrentLinkedQueue<Integer> staleQueue = new ConcurrentLinkedQueue<>();
            Threading.processRangeParallel(snapCount, ai -> {
                int k = activeArr[ai];
                if (k == finalI) return;

                double oldIK = get(finalI, k);
                double oldJK = get(finalJ, k);
                double newIK = (oldIK * wI + oldJK * wJ) / (wI + wJ);
                set(finalI, k, newIK);

                if (bestJ[k] == finalJ) {
                    staleQueue.add(k);
                } else if (newIK > bestSim[k]) {
                    bestJ[k] = finalI;
                    bestSim[k] = newIK;
                } else if (bestJ[k] == finalI && newIK < bestSim[k]) {
                    staleQueue.add(k);
                }
            });

            Integer[] staleArr = staleQueue.toArray(new Integer[0]);
            if (staleArr.length > 0) {
                Threading.processRangeParallel(staleArr.length, si -> recomputeBest(staleArr[si]));
            }
            recomputeBest(I);
        }
        return buildTree(treeIndex);
    }

    private int findBestI(int threads) {
        final int snap = activeCount;
        int actual = Math.min(threads, snap);
        int chunk = Math.max(1, (snap + actual - 1) / actual);
        int[] localI = new int[actual];
        double[] localS = new double[actual];
        Arrays.fill(localI, -1);
        Arrays.fill(localS, Double.NEGATIVE_INFINITY);

        CountDownLatch latch = new CountDownLatch(actual);
        for (int t = 0; t < actual; t++) {
            int tid = t;
            int lo = t * chunk;
            int hi = Math.min(lo + chunk, snap);
            Threading.submit(() -> {
                int bI = -1;
                double bS = Double.NEGATIVE_INFINITY;
                for (int ai = lo; ai < hi; ai++) {
                    int i = activeArr[ai];
                    if (bestSim[i] > bS) { bS = bestSim[i]; bI = i; }
                }
                localI[tid] = bI;
                localS[tid] = bS;
                latch.countDown();
            });
        }
        try { latch.await(); }
        catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException(e);
        }

        int I = -1;
        double best = Double.NEGATIVE_INFINITY;
        for (int t = 0; t < actual; t++) {
            if (localS[t] > best) { best = localS[t]; I = localI[t]; }
        }
        return I;
    }

    private void removeActive(int j) {
        int pos = posInActive[j];
        int last = activeArr[activeCount - 1];
        activeArr[pos] = last;
        posInActive[last] = pos;
        activeCount--;
    }

    private void recomputeBest(int k) {
        double bs = Double.NEGATIVE_INFINITY;
        int bj = -1;
        for (int ai = 0; ai < activeCount; ai++) {
            int m = activeArr[ai];
            if (m == k) continue;
            double s = get(k, m);
            if (s > bs) { bs = s; bj = m; }
        }
        bestSim[k] = bs;
        bestJ[k] = bj;
    }

    private long index(int i, int j) {
        if (i > j) { int t = i; i = j; j = t; }
        return (long)i * n - (long)i * (i + 1L) / 2L + j;
    }

    private double get(int i, int j) { return mat.get(index(i, j)); }
    private void set(int i, int j, double value) { mat.set(index(i, j), value); }

    private Tree buildTree(int treeIndex) {
        TreeNode root = clusterRoot[activeArr[0]];
        int[] postArr = new int[n];
        int[] cursor = {0};
        assignRanges(root, postArr, cursor);
        int[] posMap = new int[n];
        Arrays.fill(posMap, -1);
        for (int pos = 0; pos < n; pos++) posMap[postArr[pos]] = pos;
        return new Tree(treeIndex, root, postArr, posMap, n, n);
    }

    private void assignRanges(TreeNode node, int[] postArr, int[] cursor) {
        if (node.isLeaf()) {
            node.rangeStart = cursor[0];
            postArr[cursor[0]++] = node.taxonId;
            node.rangeEnd = cursor[0];
        } else {
            node.rangeStart = cursor[0];
            assignRanges(node.left, postArr, cursor);
            assignRanges(node.right, postArr, cursor);
            node.rangeEnd = cursor[0];
        }
    }
}
