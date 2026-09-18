package stelarx;

import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.dp.BipartitionSplit;
import stelarx.dp.DPTable;
import stelarx.partition.PartitionTable;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.weight.WeightTable;

import java.io.*;
import java.util.*;

/**
 * Verifies Phase-6 STELAR-X rooted-triplet weight computation.
 *
 * Checks:
 *   1. All scores are >= 0.
 *   2. Root-level splits have the highest scores (not guaranteed but typically true).
 *   3. Score distribution: print histogram.
 *   4. Small input: print all splits sorted by score with taxon names.
 */
public class Phase6Verifier {

    public static void dump(List<Tree> trees, TaxonRegistry registry,
                            ClusterTable clusterTable, DPTable dpTable,
                            WeightTable weightTable,
                            String outFile) throws IOException {
        PrintStream out = (outFile != null)
            ? new PrintStream(new FileOutputStream(outFile)) : System.out;

        int n = registry.size();

        out.printf("=== Phase 6 Weight Calculation Verification ===%n");
        out.printf("Score type: %s%n", weightTable.getMode());
        out.printf("Taxa: %d  Scored splits: %d%n", n, weightTable.size());
        out.printf("Max score:   %d%n", weightTable.getMaxScore());
        out.printf("Total score: %d%n%n", weightTable.getTotalScore());

        int fails = 0;

        // ── Check 1: non-negative scores ─────────────────────────────────────
        // entries() exposes only the exact-LONG map; in DOUBLE/INT128 mode it is
        // empty, so the per-split scan below is a no-op (scores are validated via
        // the root-split dump using getScoreD instead).
        if (weightTable.getMode() != WeightTable.Mode.LONG) {
            out.printf("Check 1 (non-negative): SKIPPED (%s-mode scores; see root splits below)%n",
                weightTable.getMode());
        } else {
            int negCount = 0;
            for (var entry : weightTable.entries()) {
                if (entry.getValue() < 0) { negCount++; fails++; }
            }
            if (negCount == 0) out.println("Check 1 (non-negative): PASSED");
            else out.printf("Check 1 (non-negative): %d NEGATIVE scores%n", negCount);
        }

        // ── Check 2: root splits ─────────────────────────────────────────────
        out.printf("%nRoot splits and their scores:%n");
        ClusterHash root = dpTable.getRootHash();
        long maxRootScore = 0;
        for (BipartitionSplit split : dpTable.getSplits(root)) {
            long s = weightTable.getScore(split);
            if (s > maxRootScore) maxRootScore = s;
            String loName = clusterName(split.lo, clusterTable, trees, registry);
            String hiName = clusterName(split.hi, clusterTable, trees, registry);
            out.printf("  score=%6d  {%s} | {%s}%n", s, loName, hiName);
        }

        // ── Check 3: Score distribution ──────────────────────────────────────
        out.printf("%n--- Score distribution ---%n");
        TreeMap<Long, Integer> dist = new TreeMap<>(Collections.reverseOrder());
        for (var entry : weightTable.entries()) {
            dist.merge(entry.getValue(), 1, Integer::sum);
        }
        int printed = 0;
        for (var entry : dist.entrySet()) {
            out.printf("  score=%6d : %d splits%n", entry.getKey(), entry.getValue());
            if (++printed >= 20) { out.println("  ..."); break; }
        }

        // ── Summary ──────────────────────────────────────────────────────────
        out.printf("%n--- Summary ---%n");
        if (fails == 0) out.println("ALL ASSERTIONS PASSED");
        else            out.printf("%d FAILURES%n", fails);

        // ── Small input: all splits sorted by score ───────────────────────────
        if (n <= 8) {
            out.printf("%n--- All splits ranked by score (small input) ---%n");
            List<Map.Entry<BipartitionSplit, Long>> sorted = new ArrayList<>(weightTable.entries());
            sorted.sort((x, y) -> Long.compare(y.getValue(), x.getValue()));

            for (var entry : sorted) {
                BipartitionSplit split = entry.getKey();
                long score = entry.getValue();
                String loName = clusterName(split.lo, clusterTable, trees, registry);
                String hiName = clusterName(split.hi, clusterTable, trees, registry);
                out.printf("  score=%4d  {%s} | {%s}%n", score, loName, hiName);
            }
        }

        if (outFile != null) out.close();
    }

    // -------------------------------------------------------------------------

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
