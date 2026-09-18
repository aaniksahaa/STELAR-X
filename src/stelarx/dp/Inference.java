package stelarx.dp;

import stelarx.Logging;
import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.partition.PartitionTable;
import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.NewickLabel;
import stelarx.tree.Tree;
import stelarx.util.Int128;
import stelarx.weight.WeightTable;

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

    private final Map<ClusterHash, Long>             dpMemo     = new HashMap<>();   // LONG score path
    private final Map<ClusterHash, Double>           dpMemoD    = new HashMap<>();   // DOUBLE score path
    private final Map<ClusterHash, Int128>           dpMemoI    = new HashMap<>();   // INT128 score path
    private final Map<ClusterHash, BipartitionSplit> bestSplits = new HashMap<>();
    private String lastTripletScore;
    /** Exact inverse for every size-one cluster hash in the analysis universe. */
    private Map<ClusterHash, Integer> singletonTaxa = Collections.emptyMap();
    private TaxonHasher reconstructionHasher;

    // -------------------------------------------------------------------------

    /**
     * Run the inference DP and return the optimal species tree as a Newick string.
     *
     * @param dpTable       DP search space (transitions)
     * @param weightTable   precomputed split scores
     * @param clusterTable  cluster exemplars used to recover split membership
     * @param exemplarTrees all trees referenced by ClusterTable exemplars,
     *                      including S3 consensus-emission snapshots
     * @param registry      taxon ID ↔ name
     * @param hasher        per-taxon hashes used to validate reconstruction
     * @return Newick string ending with ";"
     */
    public String run(DPTable dpTable, WeightTable weightTable,
                      ClusterTable clusterTable, List<Tree> exemplarTrees,
                      TaxonRegistry registry, TaxonHasher hasher) {
        long t0 = System.nanoTime();
        reconstructionHasher = hasher;
        indexSingletonTaxa(registry, hasher);
        validateExemplarCoverage(clusterTable, exemplarTrees);

        ClusterHash root = dpTable.getRootHash();

        // The score type mirrors the WeightTable's accumulation decision: exact
        // LONG for normal sizes; above the overflow threshold, exact INT128
        // (default) or approximate DOUBLE.  The [tag] makes the type explicit.
        if (weightTable.isInt128()) {
            Int128 totalScore = solveI(root, dpTable, weightTable);
            lastTripletScore = totalScore.toString();
            long ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("Inference DP: optimization-objective triplet score = %s  "
                + "[int128]  (%d ms)", totalScore, ms);
        } else if (weightTable.isDouble()) {
            double totalScore = solveD(root, dpTable, weightTable);
            lastTripletScore = String.format(java.util.Locale.ROOT, "%.0f", totalScore);
            long ms = (System.nanoTime() - t0) / 1_000_000;
            // %.0f keeps it a plain (huge) number for log parsers; the [double]
            // tag makes the active numeric type explicit in the logs.
            Logging.info("Inference DP: optimization-objective triplet score = %.0f  "
                + "[double]  (%d ms)", totalScore, ms);
        } else {
            long totalScore = solve(root, dpTable, weightTable);
            lastTripletScore = Long.toString(totalScore);
            long ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("Inference DP: optimization-objective triplet score = %d  "
                + "[long]  (%d ms)", totalScore, ms);
        }

        BitSet rootMembers = new BitSet(registry.size());
        rootMembers.set(0, registry.size());
        validateMembers(root, rootMembers);
        String newick = buildNewick(
            root, clusterTable, exemplarTrees, registry, rootMembers) + ";";
        return newick;
    }

    private void indexSingletonTaxa(TaxonRegistry registry, TaxonHasher hasher) {
        if (hasher.numTaxa() != registry.size()) {
            throw new IllegalArgumentException("Taxon hasher/registry size mismatch during "
                + "species-tree reconstruction");
        }
        Map<ClusterHash, Integer> indexed = new HashMap<>(
            Math.max(16, registry.size() * 2));
        int m = hasher.numSeeds();
        for (int taxonId = 0; taxonId < registry.size(); taxonId++) {
            long[] sums = new long[m];
            long[] xors = new long[m];
            for (int seed = 0; seed < m; seed++) {
                long value = hasher.get(seed, taxonId);
                sums[seed] = value;
                xors[seed] = value;
            }
            indexed.put(new ClusterHash(sums, xors, 1, m), taxonId);
        }
        singletonTaxa = indexed;
    }

    private static void validateExemplarCoverage(ClusterTable clusterTable,
                                                 List<Tree> exemplarTrees) {
        int maxTreeIndex = -1;
        for (ClusterTable.Entry entry : clusterTable.entries()) {
            maxTreeIndex = Math.max(maxTreeIndex, entry.exemplar.treeIndex);
        }
        if (maxTreeIndex >= exemplarTrees.size()) {
            throw new IllegalStateException("Cluster table references exemplar tree index "
                + maxTreeIndex + ", but reconstruction received only "
                + exemplarTrees.size() + " exemplar tree(s)");
        }
    }

    /** Raw rooted-triplet score from the most recent successful inference run. */
    public String getLastTripletScore() { return lastTripletScore; }

    /**
     * Score the fixed tree represented by {@code dpTable}; no tree reconstruction
     * or inferred species tree output is produced.
     *
     * @return formatted raw rooted-triplet score
     */
    public String scoreFixedTree(DPTable dpTable, WeightTable weightTable) {
        return scoreFixedTree(dpTable, weightTable, "Score-only");
    }

    /** Score a fixed tree with a caller-supplied log label. */
    public String scoreFixedTree(DPTable dpTable, WeightTable weightTable,
                                 String logLabel) {
        long t0 = System.nanoTime();
        ClusterHash root = dpTable.getRootHash();

        if (weightTable.isInt128()) {
            Int128 totalScore = solveI(root, dpTable, weightTable);
            long ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("%s: triplet score = %s  [int128]  (%d ms)",
                logLabel, totalScore, ms);
            return totalScore.toString();
        } else if (weightTable.isDouble()) {
            double totalScore = solveD(root, dpTable, weightTable);
            long ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("%s: triplet score = %.0f  [double]  (%d ms)",
                logLabel, totalScore, ms);
            return String.format("%.0f", totalScore);
        } else {
            long totalScore = solve(root, dpTable, weightTable);
            long ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("%s: triplet score = %d  [long]  (%d ms)",
                logLabel, totalScore, ms);
            return Long.toString(totalScore);
        }
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

        Set<BipartitionSplit> splits = dpTable.getSplits(ch);

        if (splits.isEmpty()) {
            // No splits available (should only happen for singletons; edge case for
            // complement clusters that don't appear as subtrees in any gene tree).
            // Return 0 so the path is still traversable (polytomy in output).
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

        dpMemo.put(ch, best);
        if (bestSp != null) bestSplits.put(ch, bestSp);
        return best;
    }

    /**
     * Floating-point mirror of {@link #solve} used when the WeightTable scores in
     * {@code double} (very large taxon sets).  Identical structure; only the
     * accumulator/comparison type differs.  Populates the shared {@code bestSplits}
     * map, so Newick reconstruction is unchanged.
     */
    private double solveD(ClusterHash ch, DPTable dpTable, WeightTable weightTable) {
        Double memo = dpMemoD.get(ch);
        if (memo != null) return memo;

        if (ch.size == 1) {
            dpMemoD.put(ch, 0.0);
            return 0.0;
        }

        Set<BipartitionSplit> splits = dpTable.getSplits(ch);
        if (splits.isEmpty()) {
            dpMemoD.put(ch, 0.0);
            return 0.0;
        }

        double best = Double.NEGATIVE_INFINITY;
        BipartitionSplit bestSp = null;

        for (BipartitionSplit split : splits) {
            double score = weightTable.getScoreD(split)
                         + solveD(split.lo, dpTable, weightTable)
                         + solveD(split.hi, dpTable, weightTable);
            if (score > best) {
                best  = score;
                bestSp = split;
            }
        }

        dpMemoD.put(ch, best);
        if (bestSp != null) bestSplits.put(ch, bestSp);
        return best;
    }

    /**
     * Exact 128-bit mirror of {@link #solve} (default path for very large taxon
     * sets).  Identical structure; accumulator/comparison use {@link Int128}.
     * Populates the shared {@code bestSplits} map, so Newick reconstruction is
     * unchanged.
     */
    private Int128 solveI(ClusterHash ch, DPTable dpTable, WeightTable weightTable) {
        Int128 memo = dpMemoI.get(ch);
        if (memo != null) return memo;

        if (ch.size == 1) {
            dpMemoI.put(ch, Int128.ZERO);
            return Int128.ZERO;
        }

        Set<BipartitionSplit> splits = dpTable.getSplits(ch);
        if (splits.isEmpty()) {
            dpMemoI.put(ch, Int128.ZERO);
            return Int128.ZERO;
        }

        Int128 best = null;   // null = lowest sentinel
        BipartitionSplit bestSp = null;

        for (BipartitionSplit split : splits) {
            Int128 score = weightTable.getScoreI(split)
                         .add(solveI(split.lo, dpTable, weightTable))
                         .add(solveI(split.hi, dpTable, weightTable));
            if (best == null || score.compareTo(best) > 0) {
                best  = score;
                bestSp = split;
            }
        }

        dpMemoI.put(ch, best);
        if (bestSp != null) bestSplits.put(ch, bestSp);
        return best;
    }

    // -------------------------------------------------------------------------
    // Newick reconstruction
    // -------------------------------------------------------------------------

    private String buildNewick(ClusterHash ch, ClusterTable clusterTable,
                               List<Tree> trees, TaxonRegistry registry,
                               BitSet members) {
        if (members.cardinality() != ch.size) {
            throw new IllegalStateException("Reconstructed cluster has "
                + members.cardinality() + " taxa, expected " + ch.size);
        }

        // Singleton
        if (ch.size == 1) {
            return taxonName(ch, registry, members);
        }

        BipartitionSplit split = bestSplits.get(ch);
        if (split == null) {
            // No chosen resolution: emit exactly the membership inherited from
            // the parent split.  This also covers valid DP residual clusters that
            // have no ClusterTable exemplar (common with incomplete gene trees).
            return polytomy(ch, registry, members);
        }

        BitSet[] childMembers = partitionMembers(
            ch, members, split, clusterTable, trees, registry.size());
        String left = buildNewick(
            split.lo, clusterTable, trees, registry, childMembers[0]);
        String right = buildNewick(
            split.hi, clusterTable, trees, registry, childMembers[1]);
        return "(" + left + "," + right + ")";
    }

    /** Get the taxon name for a singleton cluster. */
    private String taxonName(ClusterHash ch, TaxonRegistry registry,
                             BitSet members) {
        int inheritedId = members.nextSetBit(0);
        if (inheritedId >= 0 && members.nextSetBit(inheritedId + 1) < 0) {
            return NewickLabel.encode(registry.getName(inheritedId));
        }
        Integer taxonId = singletonTaxa.get(ch);
        if (taxonId == null) {
            throw new IllegalStateException(
                "Cannot resolve singleton cluster to a taxon: " + ch);
        }
        return NewickLabel.encode(registry.getName(taxonId));
    }

    /** Emit an unresolved cluster from its exact inherited membership. */
    private String polytomy(ClusterHash ch, TaxonRegistry registry,
                             BitSet members) {
        StringBuilder sb = new StringBuilder("(");
        boolean first = true;
        for (int taxonId = members.nextSetBit(0); taxonId >= 0;
                taxonId = members.nextSetBit(taxonId + 1)) {
            if (!first) sb.append(',');
            first = false;
            sb.append(NewickLabel.encode(registry.getName(taxonId)));
        }
        sb.append(')');
        return sb.toString();
    }

    /**
     * Recover one child from a concrete exemplar, then derive its sibling as the
     * exact difference from the already-known parent membership.  This avoids
     * requiring every valid DP residual cluster to have its own exemplar.
     */
    private BitSet[] partitionMembers(ClusterHash parentHash, BitSet parent,
                                      BipartitionSplit split,
                                      ClusterTable ct, List<Tree> trees,
                                      int numTaxa) {
        boolean tryLoFirst = split.lo.size <= split.hi.size;
        ClusterHash firstHash = tryLoFirst ? split.lo : split.hi;
        ClusterHash secondHash = tryLoFirst ? split.hi : split.lo;

        BitSet first = explicitMembers(firstHash, ct, trees, numTaxa);
        if (first == null || !isSubset(first, parent) || !matchesHash(first, firstHash)) {
            firstHash = secondHash;
            first = explicitMembers(firstHash, ct, trees, numTaxa);
        }
        if (first == null || !isSubset(first, parent) || !matchesHash(first, firstHash)) {
            throw new IllegalStateException("Cannot recover either child membership for "
                + "inferred split " + split + " within parent " + parent.cardinality());
        }

        BitSet second = (BitSet) parent.clone();
        second.andNot(first);
        ClusterHash derivedHash = firstHash.equals(split.lo) ? split.hi : split.lo;
        ClusterHash algebraicResidual = ClusterHash.residual(parentHash, firstHash);
        if (second.cardinality() != derivedHash.size
                || !algebraicResidual.equals(derivedHash)) {
            throw new IllegalStateException("Inferred split membership is inconsistent "
                + "with its parent: " + split);
        }

        if (firstHash.equals(split.lo)) return new BitSet[] { first, second };
        return new BitSet[] { second, first };
    }

    /** Materialize a cluster when it has a singleton identity or table exemplar. */
    private BitSet explicitMembers(ClusterHash ch, ClusterTable ct,
                                   List<Tree> trees, int numTaxa) {
        if (ch.size == 1) {
            Integer taxonId = singletonTaxa.get(ch);
            if (taxonId == null) return null;
            BitSet singleton = new BitSet(numTaxa);
            singleton.set(taxonId);
            return singleton;
        }

        ClusterTable.Entry entry = ct.get(ch);
        if (entry == null) return null;
        Cluster ex = entry.exemplar;
        if (ex.treeIndex < 0 || ex.treeIndex >= trees.size()) {
            throw new IllegalStateException("Cluster exemplar references tree index "
                + ex.treeIndex + ", but reconstruction received only " + trees.size()
                + " exemplar tree(s)");
        }
        Tree tree = trees.get(ex.treeIndex);
        BitSet explicit = new BitSet(numTaxa);
        if (ex.complement) explicit.set(0, numTaxa);

        if (ex.isMultiRange()) {
            for (int range = 0; range < ex.los.length; range++) {
                for (int pos = ex.los[range]; pos < ex.his[range]; pos++) {
                    int taxonId = tree.postorderArray[pos];
                    if (ex.complement) explicit.clear(taxonId);
                    else explicit.set(taxonId);
                }
            }
        } else {
            for (int pos = ex.left; pos < ex.right; pos++) {
                int taxonId = tree.postorderArray[pos];
                if (ex.complement) explicit.clear(taxonId);
                else explicit.set(taxonId);
            }
        }
        return explicit;
    }

    private static boolean isSubset(BitSet candidate, BitSet parent) {
        BitSet outside = (BitSet) candidate.clone();
        outside.andNot(parent);
        return outside.isEmpty();
    }

    private void validateMembers(ClusterHash expected, BitSet members) {
        if (members.cardinality() != expected.size || !matchesHash(members, expected)) {
            throw new IllegalStateException("Reconstructed cluster membership does not "
                + "match its hash: expected " + expected + ", recovered "
                + members.cardinality() + " taxa");
        }
    }

    private boolean matchesHash(BitSet members, ClusterHash expected) {
        int m = reconstructionHasher.numSeeds();
        if (expected.sums.length != m) return false;
        long[] sums = new long[m];
        long[] xors = new long[m];
        for (int taxonId = members.nextSetBit(0); taxonId >= 0;
                taxonId = members.nextSetBit(taxonId + 1)) {
            for (int seed = 0; seed < m; seed++) {
                long value = reconstructionHasher.get(seed, taxonId);
                sums[seed] += value;
                xors[seed] ^= value;
            }
        }
        return Arrays.equals(sums, expected.sums)
            && Arrays.equals(xors, expected.xors);
    }


    // -------------------------------------------------------------------------
    // Accessors for external use (e.g. scoring statistics)
    // -------------------------------------------------------------------------

    public long getDPScore(ClusterHash ch)          { return dpMemo.getOrDefault(ch, 0L); }
    public BipartitionSplit getBestSplit(ClusterHash ch) { return bestSplits.get(ch); }
}
