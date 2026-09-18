package stelarx;

import stelarx.cluster.ClusterHash;
import stelarx.hash.PrefixHashArrays;
import stelarx.partition.Partition;
import stelarx.partition.PartitionTable;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;

import java.io.*;
import java.util.*;

/**
 * Verifies Phase-4 rooted child-partition extraction.
 *
 * Key checks:
 *   1. Child sizes exactly cover the exemplar node range.
 *   2. Root partitions are retained.
 *   3. Total frequency equals the number of nodes capable of displaying a triplet.
 *   4. Child range hashes sum to the node-range hash.
 *   5. For small inputs: show taxa in each part.
 */
public class Phase4Verifier {

    public static void dump(List<Tree> trees, TaxonRegistry registry,
                            PrefixHashArrays pref, PartitionTable partTable,
                            String outFile) throws IOException {
        PrintStream out = (outFile != null)
            ? new PrintStream(new FileOutputStream(outFile)) : System.out;

        int n = registry.size();
        int k = trees.size();
        int m = pref.numSeeds();

        out.printf("=== Phase 4 Rooted Child-Partition Verification ===%n");
        out.printf("Taxa: %d  Trees: %d  Seeds: %d%n", n, k, m);
        out.printf("Unique rooted child partitions: %d%n%n", partTable.size());

        // Expected contributing internal nodes, root included. Cherries and
        // singleton-only polytomies are identically zero and are not stored.
        int expectedTotal = 0;
        for (Tree t : trees) expectedTotal += countContributingNodes(t.root);
        out.printf("Expected contributing rooted internal nodes: %d%n%n", expectedTotal);

        int fails = 0;

        // Check 1+2: child-size consistency
        for (PartitionTable.Entry e : partTable.entries()) {
            Partition p = e.exemplar;
            int totalSize;
            if (p.d == 3) {
                totalSize = p.size1 + p.size2;
            } else {
                totalSize = 0;
                for (int size : p.sizes) totalSize += size;
            }
            int nodeLo = p.d == 3 ? p.leftStart : p.partStarts[0];
            int nodeHi = p.d == 3 ? p.rightEnd : p.partEnds[p.d - 2];
            if (totalSize != nodeHi - nodeLo) {
                out.printf("FAIL: child sizes=%d != node range size=%d in %s%n",
                    totalSize, nodeHi - nodeLo, p);
                fails++;
            }
        }
        int observedTotal = partTable.entries().stream().mapToInt(e -> e.frequency).sum();
        if (observedTotal != expectedTotal) {
            out.printf("FAIL: partition frequency total=%d, expected=%d%n",
                observedTotal, expectedTotal);
            fails++;
        }

        // Check 4: additive hash consistency over actual children.
        for (PartitionTable.Entry e : partTable.entries()) {
            Partition p = e.exemplar;
            int ti = p.treeIndex;
            for (int s = 0; s < m; s++) {
                long partsSum = 0L;
                long partsXor = 0L;
                if (p.d == 3) {
                    partsSum += pref.rangeSum(ti, s, p.leftStart, p.leftEnd);
                    partsSum += pref.rangeSum(ti, s, p.rightStart, p.rightEnd);
                    partsXor ^= pref.rangeXor(ti, s, p.leftStart, p.leftEnd);
                    partsXor ^= pref.rangeXor(ti, s, p.rightStart, p.rightEnd);
                } else {
                    for (int i = 0; i < p.d - 1; i++) {
                        partsSum += pref.rangeSum(ti, s, p.partStarts[i], p.partEnds[i]);
                        partsXor ^= pref.rangeXor(ti, s, p.partStarts[i], p.partEnds[i]);
                    }
                }
                int nodeLo = p.d == 3 ? p.leftStart : p.partStarts[0];
                int nodeHi = p.d == 3 ? p.rightEnd : p.partEnds[p.d - 2];
                long nodeSum = pref.rangeSum(ti, s, nodeLo, nodeHi);
                long nodeXor = pref.rangeXor(ti, s, nodeLo, nodeHi);
                if (partsSum != nodeSum || partsXor != nodeXor) {
                    out.printf("FAIL hash s=%d: children=(sum=%x,xor=%x) "
                            + "!= node=(sum=%x,xor=%x) in %s%n",
                        s, partsSum, partsXor, nodeSum, nodeXor, p);
                    fails++;
                }
            }
        }

        out.printf("%n--- Summary ---%n");
        if (fails == 0) out.println("ALL ASSERTIONS PASSED");
        else            out.printf("%d FAILURES%n", fails);

        // Diagnostic only: outside size is derived, never stored or keyed.
        out.printf("%n--- derived outside-size distribution ---%n");
        Map<Integer, Integer> dist = new TreeMap<>();
        for (PartitionTable.Entry e : partTable.entries()) {
            Partition p = e.exemplar;
            int childSize = p.d == 3 ? p.size1 + p.size2
                : java.util.Arrays.stream(p.sizes).sum();
            dist.merge(trees.get(p.treeIndex).leafCount - childSize, 1, Integer::sum);
        }
        dist.forEach((sz, cnt) -> out.printf("  outside_size=%3d : %d partitions%n", sz, cnt));

        // Small-input: show all child partitions with taxon names
        if (n <= 8) {
            out.printf("%n--- All rooted child partitions (small input) ---%n");
            List<PartitionTable.Entry> sorted = new ArrayList<>(partTable.entries());
            sorted.sort(Comparator.comparingInt(e -> e.exemplar.size1));
            for (PartitionTable.Entry e : sorted) {
                Partition p = e.exemplar;
                Tree t = trees.get(p.treeIndex);
                if (p.d == 3) {
                    String s1 = rangeNames(t, p.leftStart, p.leftEnd, false, registry);
                    String s2 = rangeNames(t, p.rightStart, p.rightEnd, false, registry);
                    String s3 = rangeNames(t, p.leftStart, p.rightEnd, true, registry);
                    out.printf("  freq=%d  children={%s}|{%s}  outside={%s}%n",
                        e.frequency, s1, s2, s3);
                } else {
                    out.printf("  freq=%d  children=", e.frequency);
                    for (int i = 0; i < p.d - 1; i++) {
                        if (i > 0) out.print('|');
                        out.printf("{%s}", rangeNames(t, p.partStarts[i], p.partEnds[i], false, registry));
                    }
                    out.printf("  outside={%s}%n", rangeNames(t, p.partStarts[0],
                        p.partEnds[p.d - 2], true, registry));
                }
            }
        }

        if (outFile != null) out.close();
    }

    private static String rangeNames(Tree tree, int lo, int hi, boolean complement,
                                     TaxonRegistry registry) {
        StringBuilder sb = new StringBuilder();
        if (!complement) {
            for (int i = lo; i < hi; i++) {
                if (sb.length() > 0) sb.append(",");
                sb.append(registry.getName(tree.postorderArray[i]));
            }
        } else {
            for (int i = 0; i < tree.leafCount; i++) {
                if (i >= lo && i < hi) continue;
                if (sb.length() > 0) sb.append(",");
                sb.append(registry.getName(tree.postorderArray[i]));
            }
        }
        return sb.toString();
    }

    private static int countContributingNodes(stelarx.tree.TreeNode node) {
        if (node.isLeaf()) return 0;
        int count = 0;
        if (node.isPolytomous()) {
            boolean contributes = false;
            for (var child : node.children) {
                if (child.rangeEnd - child.rangeStart >= 2) contributes = true;
                count += countContributingNodes(child);
            }
            if (contributes) count++;
        } else {
            if (node.left.rangeEnd - node.left.rangeStart >= 2
                    || node.right.rangeEnd - node.right.rangeStart >= 2) count++;
            count += countContributingNodes(node.left);
            count += countContributingNodes(node.right);
        }
        return count;
    }
}
