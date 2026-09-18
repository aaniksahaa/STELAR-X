package stelarx;

import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.dp.BipartitionSplit;
import stelarx.dp.DPTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;

import java.io.*;
import java.util.*;

/**
 * Verifies Phase-5 DP search space (tree-local transitions).
 *
 * Checks:
 *   1. Root hash has at least one split.
 *   2. For every split (A → B | C): size(B) + size(C) == size(A).
 *   3. For every split, both halves are in ClusterTable (or one of them is a
 *      root-child subtree cluster, which is always valid).
 *   4. Expected transition counts match rooted binary internal nodes.
 *   5. Small input: print all splits with taxon names.
 */
public class Phase5Verifier {

    public static void dump(List<Tree> trees, TaxonRegistry registry,
                            PrefixHashArrays pref,
                            ClusterTable clusterTable, DPTable dpTable,
                            String outFile) throws IOException {
        PrintStream out = (outFile != null)
            ? new PrintStream(new FileOutputStream(outFile)) : System.out;

        int n = registry.size();
        int k = trees.size();

        out.printf("=== Phase 5 DP Search Space Verification ===%n");
        out.printf("Taxa: %d  Trees: %d  Clusters in X: %d%n", n, k, clusterTable.size());
        out.printf("Clusters with splits: %d%n", dpTable.numClusters());
        out.printf("Unique splits total:  %d%n", dpTable.numUniqueSplits());
        out.printf("Transitions emitted:  %d (before dedup)%n%n", dpTable.numEmitted());

        int fails = 0;

        // ── Check 1: Root has splits ─────────────────────────────────────────
        ClusterHash root = dpTable.getRootHash();
        if (!dpTable.hasSplits(root)) {
            out.println("FAIL: root cluster has no splits");
            fails++;
        } else {
            out.printf("Root splits: %d%n", dpTable.getSplits(root).size());
        }

        // ── Check 2: Size consistency for every split ────────────────────────
        int sizeFailCount = 0;
        for (var entry : dpTable.entries()) {
            ClusterHash parent = entry.getKey();
            for (BipartitionSplit split : entry.getValue()) {
                int expected = parent.size;
                int actual   = split.lo.size + split.hi.size;
                if (actual != expected) {
                    out.printf("FAIL size: parent.size=%d but lo.size=%d + hi.size=%d = %d  in %s%n",
                        expected, split.lo.size, split.hi.size, actual, parent);
                    sizeFailCount++;
                    fails++;
                }
            }
        }
        if (sizeFailCount == 0) out.println("Check 2 (size consistency): PASSED");
        else                    out.printf("Check 2 (size consistency): %d FAILURES%n", sizeFailCount);

        // ── Check 3: Both halves in ClusterTable (warn on misses) ────────────
        // Halves should always be in X; if not, something is wrong with extraction.
        int membershipFails = 0;
        for (var entry : dpTable.entries()) {
            for (BipartitionSplit split : entry.getValue()) {
                boolean loInX = clusterTable.contains(split.lo);
                boolean hiInX = clusterTable.contains(split.hi);
                if (!loInX || !hiInX) {
                    out.printf("FAIL membership: lo=%s(%s) hi=%s(%s)%n",
                        split.lo, loInX ? "OK" : "MISSING",
                        split.hi, hiInX ? "OK" : "MISSING");
                    membershipFails++;
                    fails++;
                }
            }
        }
        if (membershipFails == 0) out.println("Check 3 (cluster membership): PASSED");
        else out.printf("Check 3 (cluster membership): %d FAILURES%n", membershipFails);

        // ── Check 4: Expected transition counts (tree-structural) ────────────
        // Type 1 comes from each resolved binary internal node.
        // Type 2 comes from each resolved non-root node (leaf OR internal) whose
        // parent is binary and whose parent's super-complement is nonempty.
        // Type 3 connects an incomplete tree's taxon boundary to S.
        int expType1 = 0;
        for (Tree t : trees) {
            expType1 += countType1(t.root, -1, false);
        }
        out.printf("%nExpected rooted child transitions (incl. root): %d%n", expType1);
        int expectedEmitted = expType1;
        out.printf("Total emitted: %d  (expected %d)%n",
            dpTable.numEmitted(), expectedEmitted);
        if (dpTable.numEmitted() != expectedEmitted) {
            out.println("FAIL: emitted count mismatch");
            fails++;
        } else {
            out.println("Check 4 (transition count): PASSED");
        }

        // ── Summary ──────────────────────────────────────────────────────────
        out.printf("%n--- Summary ---%n");
        if (fails == 0) out.println("ALL ASSERTIONS PASSED");
        else            out.printf("%d FAILURES%n", fails);

        // ── Split count distribution ─────────────────────────────────────────
        out.printf("%n--- Splits-per-cluster distribution ---%n");
        Map<Integer, Integer> dist = new TreeMap<>();
        for (var entry : dpTable.entries()) {
            dist.merge(entry.getValue().size(), 1, Integer::sum);
        }
        dist.forEach((cnt, num) -> out.printf("  splits=%2d : %d clusters%n", cnt, num));

        // ── Small input: print all splits with taxa names ─────────────────────
        if (n <= 8) {
            out.printf("%n--- All DP transitions (small input) ---%n");
            // Sort by parent size then parent hash
            var allEntries = new ArrayList<>(dpTable.entries());
            allEntries.sort(Comparator.comparingInt(e -> e.getKey().size));

            for (var entry : allEntries) {
                ClusterHash parent = entry.getKey();
                String parentName = clusterName(parent, clusterTable, trees, registry);
                for (BipartitionSplit split : entry.getValue()) {
                    String loName = clusterName(split.lo, clusterTable, trees, registry);
                    String hiName = clusterName(split.hi, clusterTable, trees, registry);
                    out.printf("  {%s} -> {%s} | {%s}%n", parentName, loName, hiName);
                }
            }
        }

        if (outFile != null) out.close();
    }

    // -------------------------------------------------------------------------

    /** Count resolved internal nodes that emit Type 1. */
    private static int countType1(stelarx.tree.TreeNode u, int anchorPos,
                                  boolean anchorFree) {
        if (u.isLeaf()) return 0;
        int count = 0;
        if (!u.isPolytomous()
                && (!anchorFree || !containsPosition(u, anchorPos))) {
            count = 1;
        }
        if (u.isPolytomous()) {
            for (var child : u.children) {
                count += countType1(child, anchorPos, anchorFree);
            }
        } else {
            count += countType1(u.left, anchorPos, anchorFree);
            count += countType1(u.right, anchorPos, anchorFree);
        }
        return count;
    }

    /** Count nodes (including leaves) that emit a nondegenerate Type 2. */
    private static int countType2(stelarx.tree.TreeNode u, int n, int anchorPos,
                                  boolean anchorFree) {
        int count = 0;
        if (!u.isRoot() && !u.isPolytomous() && !u.parent.isPolytomous()
                && n - u.parent.rangeSize() > 0
                && (!anchorFree || containsPosition(u, anchorPos))) {
            count = 1;
        }
        if (!u.isLeaf()) {
            if (u.isPolytomous()) {
                for (var child : u.children) {
                    count += countType2(child, n, anchorPos, anchorFree);
                }
            } else {
                count += countType2(u.left, n, anchorPos, anchorFree);
                count += countType2(u.right, n, anchorPos, anchorFree);
            }
        }
        return count;
    }

    private static boolean containsPosition(stelarx.tree.TreeNode u, int position) {
        return position >= u.rangeStart && position < u.rangeEnd;
    }

    /**
     * Return a comma-separated taxon name string for the given cluster hash.
     * Looks up the exemplar in ClusterTable; falls back to "(root)" for the
     * all-taxa cluster.
     */
    private static String clusterName(ClusterHash hash, ClusterTable ct,
                                       List<Tree> trees, TaxonRegistry registry) {
        ClusterTable.Entry entry = ct.get(hash);
        if (entry == null) return "(root)";
        Cluster ex = entry.exemplar;
        Tree t = trees.get(ex.treeIndex);
        StringBuilder sb = new StringBuilder();
        if (!ex.complement) {
            for (int i = ex.left; i < ex.right; i++) {
                if (sb.length() > 0) sb.append(',');
                sb.append(registry.getName(t.postorderArray[i]));
            }
        } else {
            for (int i = 0; i < t.leafCount; i++) {
                if (i >= ex.left && i < ex.right) continue;
                if (sb.length() > 0) sb.append(',');
                sb.append(registry.getName(t.postorderArray[i]));
            }
        }
        return sb.toString();
    }
}
