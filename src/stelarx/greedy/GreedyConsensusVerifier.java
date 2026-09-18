package stelarx.greedy;

import stelarx.cluster.ClusterTable;
import stelarx.completion.SimilarityMatrix;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Verifier for the greedy-consensus phase (Part I).
 *
 * Walks every unique bipartition in frequency-descending order, applying it
 * to both the fast {@link LaminarForest}+{@link LaminarBuilder} path and the
 * brute-force {@link LaminarOracle}.  Asserts:
 *
 *   1. Per-INSERT outcome agreement (ACCEPT / SKIP / REJECT).
 *   2. At every threshold boundary AND at the very end, the fast snapshot
 *      and the oracle's current state share the same canonical leaf-set
 *      representation — i.e. the two trees encode the same laminar family.
 *
 * Also dumps the 7 snapshot trees in Newick form to {@code outFile} (or stdout)
 * so they can be diffed against ASTRAL-MP's {@code allGreedies[i]} output.
 */
public final class GreedyConsensusVerifier {

    /**
     * @param geneTrees     gene trees only (no UPGMA), per the ASTRAL-MP-faithful
     *                      counting path.  Used for both bipartition counting
     *                      and exemplar taxa enumeration.
     * @param hasher        per-taxon hashes; same instance used to build pref.
     * @param sim           species similarity matrix (Phase 1b output); pass
     *                      null when --autocomplete is off — Step A is skipped.
     */
    public static void dump(List<Tree> geneTrees, TaxonRegistry registry,
                            ClusterTable clusterTable, PrefixHashArrays pref,
                            TaxonHasher hasher, SimilarityMatrix sim,
                            String outFile) throws IOException {
        PrintStream out = (outFile != null)
            ? new PrintStream(new FileOutputStream(outFile)) : System.out;

        int n = registry.size();
        int k = geneTrees.size();
        out.printf("=== Phase GC (Greedy Consensus) Verification ===%n");
        out.printf("Taxa: %d  Gene trees: %d  (UPGMA excluded)%n", n, k);
        out.printf("Cluster-side entries in X (incl. UPGMA): %d%n%n", clusterTable.size());

        List<Bipartition> bps = BipartitionCounter.collectFromGeneTrees(
            geneTrees, pref, clusterTable.getAllTaxaHash(), n);
        out.printf("Unique bipartitions: %d%n", bps.size());

        // ── Lockstep run of fast path + oracle ──
        LaminarForest forest    = new LaminarForest(n);
        LaminarBuilder fast     = new LaminarBuilder(forest, geneTrees, n);
        LaminarOracle oracle    = new LaminarOracle(n, geneTrees);

        ConsensusTree[] fastSnaps   = new ConsensusTree[GreedyConsensus.THRESHOLDS.length];
        String[]        oracleSnaps = new String[GreedyConsensus.THRESHOLDS.length];

        int ti = GreedyConsensus.THRESHOLDS.length - 1;
        double threshold = GreedyConsensus.THRESHOLDS[ti];

        int outcomeMismatches = 0;
        int snapshotMismatches = 0;
        int processed = 0;

        for (Bipartition b : bps) {
            double ratio = (double) b.frequency / (double) k;
            while (ti >= 0 && threshold > ratio) {
                fastSnaps[ti]   = ConsensusTree.snapshot(forest, hasher);
                oracleSnaps[ti] = oracle.canonicalLeafSets();
                String fastSets = fastSnaps[ti].canonicalLeafSets();
                if (!fastSets.equals(oracleSnaps[ti])) {
                    snapshotMismatches++;
                    out.printf("FAIL snapshot ti=%d threshold=%.4f%n", ti, threshold);
                    printDiff(out, fastSets, oracleSnaps[ti]);
                }
                ti--;
                if (ti < 0) break;
                threshold = GreedyConsensus.THRESHOLDS[ti];
            }

            LaminarBuilder.Outcome fo = fast.insert(b);
            LaminarBuilder.Outcome oo = oracle.insert(b);
            if (fo != oo) {
                outcomeMismatches++;
                if (outcomeMismatches <= 10) {
                    out.printf("FAIL outcome  bp#%d  size=%d freq=%d  fast=%s oracle=%s%n",
                        processed, b.size, b.frequency, fo, oo);
                }
            }
            processed++;
        }

        // Drain remaining lower thresholds
        while (ti >= 0) {
            fastSnaps[ti]   = ConsensusTree.snapshot(forest, hasher);
            oracleSnaps[ti] = oracle.canonicalLeafSets();
            String fastSets = fastSnaps[ti].canonicalLeafSets();
            if (!fastSets.equals(oracleSnaps[ti])) {
                snapshotMismatches++;
                out.printf("FAIL snapshot ti=%d threshold=%.4f (drain)%n",
                    ti, GreedyConsensus.THRESHOLDS[ti]);
                printDiff(out, fastSets, oracleSnaps[ti]);
            }
            ti--;
        }

        // ── Per-threshold stats + Newick dump ──
        out.printf("%n--- Per-threshold snapshots ---%n");
        for (int i = 0; i < fastSnaps.length; i++) {
            ConsensusTree s = fastSnaps[i];
            out.printf("T[%d] threshold=%.4f  internal=%d  polytomies=%d%n",
                i, GreedyConsensus.THRESHOLDS[i],
                s.numInternalNodes(), s.numPolytomies());
        }

        // ── §13.5 Cross-source signature parity check ──
        // For each non-root internal node of the densest snapshot T[0], compute
        // (a) the prefix-scan signature  σ_prefix = (sigma1, sigma2)
        // (b) the direct-sum signature   σ_direct = sum/xor hasher.get(s, t)
        //                                            over taxa in the range
        // (c) lookup in ClusterTable        the same bipartition should be in X
        //                                    because the consensus tree was
        //                                    refined from gene-tree bipartitions
        // All three must agree.
        int parityChecked = 0, parityPrefixFails = 0, parityCtFails = 0;
        ConsensusTree dense = fastSnaps[0];
        int seedCount = hasher.numSeeds();
        final int[] ck = {parityChecked}, pf = {parityPrefixFails}, cf = {parityCtFails};
        dense.forEachInternalNode(node -> {
            if (node == dense.root()) return;       // skip virtual root
            int lo = node.rangeLo(), hi = node.rangeHi();
            int sz = hi - lo;
            if (sz < 2 || sz >= n) return;          // trivial
            ck[0]++;

            long[] rawSums = new long[seedCount];
            long[] rawXors = new long[seedCount];

            for (int s = 0; s < seedCount; s++) {
                long sigPref1 = dense.sigma1(s, lo, hi);
                long sigPref2 = dense.sigma2(s, lo, hi);
                long sigDir1 = 0L, sigDir2 = 0L;
                for (int p = lo; p < hi; p++) {
                    long h = hasher.get(s, dense.aCons()[p]);
                    sigDir1 += h;
                    sigDir2 ^= h;
                }
                if (sigPref1 != sigDir1 || sigPref2 != sigDir2) {
                    pf[0]++;
                    if (pf[0] <= 3) {
                        out.printf("FAIL prefix-scan parity  node lo=%d hi=%d seed=%d  "
                                   + "prefix=(%x,%x)  direct=(%x,%x)%n",
                                   lo, hi, s, sigPref1, sigPref2, sigDir1, sigDir2);
                    }
                }
                rawSums[s] = sigPref1;
                rawXors[s] = sigPref2;
            }

            stelarx.cluster.ClusterHash ch =
                new stelarx.cluster.ClusterHash(rawSums, rawXors, sz, seedCount);
            if (!clusterTable.contains(ch)) {
                cf[0]++;
                if (cf[0] <= 3) {
                    out.printf("FAIL cross-source: consensus bipartition not in X  "
                               + "lo=%d hi=%d size=%d  (one side)%n", lo, hi, sz);
                }
            }
        });
        out.printf("%nCross-source signature parity (T[0] non-root internals): "
                   + "checked=%d  prefix-scan-fails=%d  X-lookup-fails=%d%n",
                   ck[0], pf[0], cf[0]);

        // ── §8.2 Polytomy pool + size limit ──
        PolytomyPool pool = PolytomyPool.build(fastSnaps, n, /*explicitLimit=*/0);
        out.printf("%n--- Polytomy pool (§8.2) ---%n");
        out.printf("Budget N = %d + n*%d = %d%n",
            PolytomyPool.BUDGET_MIN, PolytomyPool.BUDGET_MULT, pool.budgetN);
        out.printf("Total polytomies across all 7 snapshots:  %d%n", pool.numTotal);
        out.printf("Size limit (max accepted degree):         %d%n", pool.sizeLimit);
        out.printf("Sum of squares accumulated:               %d%n", pool.sumSquaresAccumulated);
        out.printf("Accepted polytomies (degree ≤ limit):     %d%n", pool.numAccepted);
        out.printf("Skipped  polytomies (degree > limit):     %d%n", pool.numSkipped);
        if (pool.numTotal > 0) {
            out.printf("Degree histogram (degree: count):%n");
            for (int d = 3; d < pool.degreeHistogram.length; d++) {
                if (pool.degreeHistogram[d] > 0) {
                    out.printf("  d=%d: %d%n", d, pool.degreeHistogram[d]);
                }
            }
        }

        // ── §8.3 + §8.4 + §10: per-polytomy Step A + Step B in parallel ──
        out.printf("%n--- Polytomy resolution: Step A + Step B (parallel, LPT) ---%n");
        if (sim == null) {
            out.println("SKIPPED — SimilarityMatrix not available "
                + "(--autocomplete-incomplete-gene-trees was off).");
        } else {
            EmissionBuffer buffer = new EmissionBuffer();
            long t0 = System.nanoTime();
            int total = PolytomyResolver.runAllParallel(
                pool.tasks, geneTrees, sim, buffer, n, /*baseSeed=*/692L);
            long ms = (System.nanoTime() - t0) / 1_000_000;

            int countA = 0, countB = 0;
            for (EmittedBipartition b : buffer.all()) {
                if (b.source == 'A') countA++;
                else if (b.source == 'B') countB++;
            }
            out.printf("Polytomies processed:                %d%n", pool.tasks.size());
            out.printf("Total emissions (raw, sum A+B):      %d%n", total);
            out.printf("Buffer Step A signatures (deduped):  %d%n", countA);
            out.printf("Buffer Step B signatures (deduped):  %d%n", countB);
            out.printf("Buffer A∪B (deduped):                %d%n", buffer.size());
            out.printf("Parallel dispatch time:              %d ms%n", ms);

            int alreadyInX = 0;
            for (EmittedBipartition b : buffer.all()) {
                if (clusterTable.contains(b.signature)) alreadyInX++;
            }
            out.printf("Already in ClusterTable:             %d  /  net-new: %d%n",
                alreadyInX, buffer.size() - alreadyInX);

            // Dump A and B emissions separately so they can be diffed against ASTRAL-MP.
            out.printf("%n--- Step A bipartitions (canonical taxa sets) ---%n");
            List<String> linesA = new ArrayList<>();
            for (EmittedBipartition b : buffer.all()) {
                if (b.source == 'A') linesA.add(taxaSetLine("STEPA", b, registry));
            }
            java.util.Collections.sort(linesA);
            for (String l : linesA) out.println(l);

            out.printf("%n--- Step B bipartitions (canonical taxa sets) ---%n");
            List<String> linesB = new ArrayList<>();
            for (EmittedBipartition b : buffer.all()) {
                if (b.source == 'B') linesB.add(taxaSetLine("STEPB", b, registry));
            }
            java.util.Collections.sort(linesB);
            for (String l : linesB) out.println(l);
        }

        out.printf("%n--- Newick (for ASTRAL-MP head-to-head) ---%n");
        for (int i = 0; i < fastSnaps.length; i++) {
            out.printf("T[%d]_threshold_%.4f: %s%n",
                i, GreedyConsensus.THRESHOLDS[i],
                fastSnaps[i].toNewick(registry));
        }

        // ── Summary ──
        out.printf("%n--- Summary ---%n");
        out.printf("Bipartitions processed: %d%n", processed);
        out.printf("INSERT outcome mismatches:  %d%n", outcomeMismatches);
        out.printf("Snapshot leaf-set mismatches: %d%n", snapshotMismatches);
        if (outcomeMismatches == 0 && snapshotMismatches == 0) {
            out.println("ALL ASSERTIONS PASSED (fast == oracle)");
        } else {
            out.println("FAILURES PRESENT — see above");
        }

        if (outFile != null) out.close();
    }

    /** Build "[TAG] ti=... size=... {taxa}" line for one emission. */
    private static String taxaSetLine(String tag, EmittedBipartition b, TaxonRegistry registry) {
        StringBuilder sb = new StringBuilder("[").append(tag).append("] ti=");
        sb.append(b.thresholdIndex).append("  size=").append(b.size).append("  {");
        java.util.TreeSet<String> taxa = new java.util.TreeSet<>();
        int[] arr = b.canonicalSide.tree.aCons();
        for (int r = 0; r < b.canonicalSide.numRanges(); r++) {
            for (int p = b.canonicalSide.los[r]; p < b.canonicalSide.his[r]; p++) {
                taxa.add(registry.getName(arr[p]));
            }
        }
        sb.append(String.join(",", taxa)).append('}');
        return sb.toString();
    }

    /** Print up to a few diff lines so a mismatch isn't a wall of text. */
    private static void printDiff(PrintStream out, String fast, String oracle) {
        String[] fLines = fast.split("\n");
        String[] oLines = oracle.split("\n");
        java.util.Set<String> fSet = new java.util.HashSet<>(java.util.Arrays.asList(fLines));
        java.util.Set<String> oSet = new java.util.HashSet<>(java.util.Arrays.asList(oLines));

        int shown = 0;
        for (String s : fLines) {
            if (!oSet.contains(s) && shown < 10) { out.printf("  + (fast only)   %s%n", s); shown++; }
        }
        for (String s : oLines) {
            if (!fSet.contains(s) && shown < 20) { out.printf("  - (oracle only) %s%n", s); shown++; }
        }
    }

    private GreedyConsensusVerifier() {}
}
