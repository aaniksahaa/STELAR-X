package stelarx;

import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;

import java.io.*;
import java.util.*;

/**
 * Verifies Phase-3 cluster extraction.
 *
 * Key checks:
 *   1. Every cluster has size in [1, n-1].
 *   2. No stored rooted cluster is represented as a complement.
 *   3. Size bins are consistent with stored cluster sizes.
 */
public class Phase3Verifier {

    public static void dump(List<Tree> trees, TaxonRegistry registry,
                            PrefixHashArrays pref, ClusterTable clusterTable,
                            String outFile) throws IOException {
        PrintStream out = (outFile != null)
            ? new PrintStream(new FileOutputStream(outFile)) : System.out;

        int n   = registry.size();
        int k   = trees.size();
        int m   = pref.numSeeds();

        out.printf("=== Phase 3 Cluster Extraction Verification ===%n");
        out.printf("Taxa: %d  Trees: %d  Seeds: %d%n", n, k, m);
        out.printf("Unique clusters in X: %d%n", clusterTable.size());
        out.printf("All-taxa hash: %s%n%n", clusterTable.getAllTaxaHash());

        int fails = 0;

        // Check 1: all sizes in [1, n-1]
        for (ClusterTable.Entry e : clusterTable.entries()) {
            int sz = e.hash.size;
            if (sz < 1 || sz >= n) {
                out.printf("FAIL: cluster size %d out of [1,%d)%n", sz, n);
                fails++;
            }
        }

        // Check 2: rooted X contains descendant clades, never complement orientations.
        for (ClusterTable.Entry e : clusterTable.entries()) {
            if (e.exemplar.complement) {
                out.printf("FAIL: rooted X contains complement exemplar %s%n", e.hash);
                fails++;
            }
        }

        // Check 3: size-bin consistency
        for (int sz : clusterTable.sizes()) {
            for (ClusterHash h : clusterTable.getBySize(sz)) {
                ClusterTable.Entry e = clusterTable.get(h);
                if (e == null) {
                    out.printf("FAIL: size bin %d contains orphan hash%n", sz);
                    fails++;
                } else if (e.hash.size != sz) {
                    out.printf("FAIL: size bin %d has entry with size %d%n", sz, e.hash.size);
                    fails++;
                }
            }
        }

        // Summary
        out.printf("%n--- Summary ---%n");
        if (fails == 0) {
            out.println("ALL ASSERTIONS PASSED");
        } else {
            out.printf("%d FAILURES%n", fails);
        }

        // Size distribution
        out.printf("%n--- Size distribution ---%n");
        List<Integer> sizes = new ArrayList<>(clusterTable.sizes());
        Collections.sort(sizes);
        for (int sz : sizes) {
            out.printf("  size %3d: %d clusters%n", sz, clusterTable.getBySize(sz).size());
        }

        // Small-input: show all clusters with their taxon sets
        if (n <= 8) {
            out.printf("%n--- All clusters (small input) ---%n");
            List<ClusterTable.Entry> sorted = new ArrayList<>(clusterTable.entries());
            sorted.sort(Comparator.comparingInt(e -> e.hash.size));
            for (ClusterTable.Entry e : sorted) {
                String taxa = taxaInCluster(e.exemplar, trees.get(e.exemplar.treeIndex), registry);
                out.printf("  sz=%d freq=%d  taxa={%s}  %s%n",
                    e.hash.size, e.frequency, taxa, e.exemplar);
            }
        }

        if (outFile != null) out.close();
    }

    /** Enumerate the actual taxon names in a cluster from its exemplar. */
    static String taxaInCluster(Cluster c, Tree tree, TaxonRegistry registry) {
        StringBuilder sb = new StringBuilder();
        int[] arr = tree.postorderArray;
        if (!c.complement) {
            for (int i = c.left; i < c.right; i++) {
                if (sb.length() > 0) sb.append(",");
                sb.append(registry.getName(arr[i]));
            }
        } else {
            for (int i = 0; i < tree.leafCount; i++) {
                if (i >= c.left && i < c.right) continue;
                if (sb.length() > 0) sb.append(",");
                sb.append(registry.getName(arr[i]));
            }
        }
        return sb.toString();
    }
}
