package stelarx;

import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;

import java.io.*;
import java.util.List;

/**
 * Dumps and verifies Phase-2 output: taxon hashes and prefix arrays.
 *
 * Key checks:
 *   1. rangeSum(t,s,0,L) == totalSum(t,s)          -- prefix covers full tree
 *   2. rangeSum(t,s,l,r) + compSum(t,s,l,r)
 *      == totalSum(t,s)                             -- sum + complement = total
 *   3. rangeXor(t,s,l,r) ^ compXor(t,s,l,r)
 *      == totalXor(t,s)                             -- XOR complement identity
 *   4. Complete trees: totalSum == allTaxaSum        -- all-taxa hash consistent
 *   5. Linearity: rangeSum of union = sum of parts  -- additive decomposition
 */
public class Phase2Verifier {

    public static void dump(List<Tree> trees, TaxonRegistry registry,
                            TaxonHasher hasher, PrefixHashArrays pref,
                            String outFile) throws IOException {
        PrintStream out = (outFile != null)
            ? new PrintStream(new FileOutputStream(outFile)) : System.out;

        int m = hasher.numSeeds();
        int n = registry.size();

        out.printf("=== Phase 2 Hash Verification ===%n");
        out.printf("Taxa: %d  Trees: %d  Seeds: %d%n%n", n, trees.size(), m);

        // ── Taxon hash table ─────────────────────────────────────────────────
        out.println("--- Taxon hashes (first seed only) ---");
        for (int id = 0; id < n; id++) {
            out.printf("  %-6s id=%-2d  h0=%016x%n",
                registry.getName(id), id, hasher.get(0, id));
        }
        out.println();

        // ── Per-tree checks ───────────────────────────────────────────────────
        int totalFail = 0;
        for (Tree t : trees) {
            int L = t.leafCount;
            int ti = t.treeIndex;
            int fails = 0;

            // Check 1: full-range sum/xor == total
            for (int s = 0; s < m; s++) {
                long rSum = pref.rangeSum(ti, s, 0, L);
                long rXor = pref.rangeXor(ti, s, 0, L);
                if (rSum != pref.totalSum(ti, s)) {
                    out.printf("FAIL tree=%d s=%d: rangeSum(0,L)=%x != totalSum=%x%n",
                        ti, s, rSum, pref.totalSum(ti, s));
                    fails++;
                }
                if (rXor != pref.totalXor(ti, s)) {
                    out.printf("FAIL tree=%d s=%d: rangeXor(0,L)=%x != totalXor=%x%n",
                        ti, s, rXor, pref.totalXor(ti, s));
                    fails++;
                }
            }

            // Check 2+3: complement identity for every position split
            for (int split = 1; split < L; split++) {
                for (int s = 0; s < m; s++) {
                    long rs = pref.rangeSum(ti, s, 0, split);
                    long cs = pref.compSum(ti, s, 0, split);
                    if (rs + cs != pref.totalSum(ti, s)) {
                        out.printf("FAIL tree=%d s=%d split=%d: rs+cs != total%n", ti, s, split);
                        fails++;
                    }
                    long rx = pref.rangeXor(ti, s, 0, split);
                    long cx = pref.compXor(ti, s, 0, split);
                    if ((rx ^ cx) != pref.totalXor(ti, s)) {
                        out.printf("FAIL tree=%d s=%d split=%d: rx^cx != totalXor%n", ti, s, split);
                        fails++;
                    }
                }
            }

            // Check 4: complete trees have totalSum == allTaxaSum
            if (t.isComplete) {
                for (int s = 0; s < m; s++) {
                    if (pref.totalSum(ti, s) != pref.allTaxaSum(s)) {
                        out.printf("FAIL tree=%d s=%d: complete tree total != allTaxaSum%n", ti, s);
                        fails++;
                    }
                }
            }

            // Check 5: linearity -- sum of two sub-ranges equals full range
            if (L >= 3) {
                int mid = L / 2;
                for (int s = 0; s < m; s++) {
                    long left  = pref.rangeSum(ti, s, 0, mid);
                    long right = pref.rangeSum(ti, s, mid, L);
                    long full  = pref.rangeSum(ti, s, 0, L);
                    if (left + right != full) {
                        out.printf("FAIL tree=%d s=%d: linearity left+right != full%n", ti, s);
                        fails++;
                    }
                }
            }

            if (fails == 0) {
                if (Logging.isDebug())
                    out.printf("Tree %3d: ALL CHECKS PASSED%n", ti);
            } else {
                out.printf("Tree %3d: %d FAILURES%n", ti, fails);
                totalFail += fails;
            }
        }

        out.printf("%n--- Summary ---%n");
        if (totalFail == 0)
            out.println("ALL ASSERTIONS PASSED across all trees and seeds.");
        else
            out.printf("%d TOTAL FAILURES%n", totalFail);

        // ── Small-tree detailed dump ──────────────────────────────────────────
        if (trees.size() <= 6) {
            out.println();
            out.println("--- Detailed prefix arrays (small input) ---");
            for (Tree t : trees) {
                int ti = t.treeIndex;
                out.printf("Tree %d postorder: ", ti);
                for (int j = 0; j < t.leafCount; j++)
                    out.printf("%s ", registry.getName(t.postorderArray[j]));
                out.println();
                for (int s = 0; s < m; s++) {
                    out.printf("  seed%d prefSum: ", s);
                    for (int j = 0; j <= t.leafCount; j++)
                        out.printf("%016x ", pref.rangeSum(ti, s, 0, j));
                    out.println();
                    out.printf("  seed%d prefXor: ", s);
                    for (int j = 0; j <= t.leafCount; j++)
                        out.printf("%016x ", pref.rangeXor(ti, s, 0, j));
                    out.println();
                }
            }
        }

        if (outFile != null) out.close();
    }
}
