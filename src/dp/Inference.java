package dp;

import utils.Logging;
import cluster.Cluster;
import cluster.ClusterHash;
import cluster.ClusterTable;
import taxon.TaxonRegistry;
import tree.Tree;
import weight.WeightTable;

import java.util.*;

/**
 * Memoized top-down inference DP.
 *
 * dp[cluster] = max over all splits (A|B) in dpTable:
 *                   score(A|B) + dp[A] + dp[B]
 * dp[singleton] = 0
 *
 * Runs in O(|DPTable| * max_splits_per_cluster) time with memoization.
 * Reconstructs the optimal species tree as a Newick string.
 */
public class Inference {

    private final Map<ClusterHash, Long>             dpMemo     = new HashMap<>();
    private final Map<ClusterHash, BipartitionSplit> bestSplits = new HashMap<>();

    // -------------------------------------------------------------------------

    /**
     * Run the inference DP and return the optimal species tree as a Newick string.
     *
     * @param dpTable       DP search space (transitions)
     * @param weightTable   precomputed split scores
     * @param clusterTable  for taxon-name lookups
     * @param trees         gene trees (for exemplar access)
     * @param registry      taxon ID to name mapping
     * @return Newick string ending with ";"
     */
    public String run(DPTable dpTable, WeightTable weightTable,
                      ClusterTable clusterTable, List<Tree> trees,
                      TaxonRegistry registry) {
        long t0 = System.nanoTime();

        ClusterHash root = dpTable.getRootHash();
        long totalScore = solve(root, dpTable, weightTable);

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Inference DP: optimal triplet score = %d  (%d ms)", totalScore, ms);
        System.out.println("Inference DP: optimal triplet score = " + totalScore + "  (" + ms + " ms)");

        String newick = buildNewick(root, dpTable, clusterTable, trees, registry) + ";";
        return newick;
    }

    // -------------------------------------------------------------------------
    // DP
    // -------------------------------------------------------------------------

    private long solve(ClusterHash ch, DPTable dpTable, WeightTable weightTable) {
        Long memo = dpMemo.get(ch);
        if (memo != null) return memo;

        // Base case: singleton taxon
        if (ch.size == 1) {
            dpMemo.put(ch, 0L);
            return 0L;
        }

        // Base case: pair of taxa - no further split needed
        if (ch.size == 2) {
            dpMemo.put(ch, 0L);
            return 0L;
        }

        Set<BipartitionSplit> splits = dpTable.getSplits(ch);

        if (splits.isEmpty()) {
            // No splits available - edge case for clusters that don't appear as subtrees.
            // Return 0 so the path is still traversable.
            dpMemo.put(ch, 0L);
            return 0L;
        }

        long best = Long.MIN_VALUE;
        BipartitionSplit bestSp = null;

        for (BipartitionSplit split : splits) {
            long score = weightTable.getScore(split)
                       + solve(split.lo, dpTable, weightTable)
                       + solve(split.hi, dpTable, weightTable);
            if (score > best) {
                best  = score;
                bestSp = split;
            }
        }

        dpMemo.put(ch, best < 0 ? 0L : best);
        if (bestSp != null) bestSplits.put(ch, bestSp);
        return best < 0 ? 0L : best;
    }

    // -------------------------------------------------------------------------
    // Newick reconstruction
    // -------------------------------------------------------------------------

    private String buildNewick(ClusterHash ch, DPTable dpTable,
                                ClusterTable clusterTable,
                                List<Tree> trees, TaxonRegistry registry) {
        // Singleton
        if (ch.size == 1) {
            return taxonName(ch, clusterTable, trees, registry);
        }

        // Pair of taxa - build leaf nodes directly
        if (ch.size == 2) {
            BipartitionSplit split = bestSplits.get(ch);
            if (split != null) {
                String left  = buildNewick(split.lo, dpTable, clusterTable, trees, registry);
                String right = buildNewick(split.hi, dpTable, clusterTable, trees, registry);
                return "(" + left + "," + right + ")";
            }
            // Fallback: use the splits from dpTable
            Set<BipartitionSplit> splits = dpTable.getSplits(ch);
            if (!splits.isEmpty()) {
                BipartitionSplit sp = splits.iterator().next();
                String left  = buildNewick(sp.lo, dpTable, clusterTable, trees, registry);
                String right = buildNewick(sp.hi, dpTable, clusterTable, trees, registry);
                return "(" + left + "," + right + ")";
            }
            // Last resort polytomy
            return polytomy(ch, clusterTable, trees, registry);
        }

        BipartitionSplit split = bestSplits.get(ch);
        if (split == null) {
            // No split found: output polytomy with all taxa in this cluster
            return polytomy(ch, clusterTable, trees, registry);
        }

        String left  = buildNewick(split.lo, dpTable, clusterTable, trees, registry);
        String right = buildNewick(split.hi, dpTable, clusterTable, trees, registry);
        return "(" + left + "," + right + ")";
    }

    /** Get the taxon name for a singleton cluster. */
    private String taxonName(ClusterHash ch, ClusterTable ct,
                              List<Tree> trees, TaxonRegistry registry) {
        ClusterTable.Entry entry = ct.get(ch);
        if (entry == null) return "?";
        Cluster ex = entry.exemplar;
        Tree t = trees.get(ex.treeIndex);
        if (!ex.complement) {
            return registry.getName(t.postorderArray[ex.left]);
        } else {
            // Complement singleton: find the one taxon NOT in [left, right)
            for (int i = 0; i < t.leafCount; i++) {
                if (i < ex.left || i >= ex.right)
                    return registry.getName(t.postorderArray[i]);
            }
            return "?";
        }
    }

    /**
     * Fallback for clusters with no split: emit a star (polytomy) by listing
     * all taxa. Gets taxa from the exemplar cluster.
     */
    private String polytomy(ClusterHash ch, ClusterTable ct,
                             List<Tree> trees, TaxonRegistry registry) {
        ClusterTable.Entry entry = ct.get(ch);
        if (entry == null) return "?";
        Cluster ex = entry.exemplar;
        Tree t = trees.get(ex.treeIndex);
        StringBuilder sb = new StringBuilder("(");
        boolean first = true;
        if (!ex.complement) {
            for (int i = ex.left; i < ex.right; i++) {
                if (!first) sb.append(','); first = false;
                sb.append(registry.getName(t.postorderArray[i]));
            }
        } else {
            for (int i = 0; i < t.leafCount; i++) {
                if (i >= ex.left && i < ex.right) continue;
                if (!first) sb.append(','); first = false;
                sb.append(registry.getName(t.postorderArray[i]));
            }
        }
        sb.append(')');
        return sb.toString();
    }

    // -------------------------------------------------------------------------
    // Accessors for external use (e.g. scoring statistics)
    // -------------------------------------------------------------------------

    public long getDPScore(ClusterHash ch)              { return dpMemo.getOrDefault(ch, 0L); }
    public BipartitionSplit getBestSplit(ClusterHash ch) { return bestSplits.get(ch); }
}
