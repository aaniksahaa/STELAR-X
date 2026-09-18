package stelarx.greedy;

import stelarx.Logging;
import stelarx.cluster.Cluster;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Gene-tree polytomy X-enrichment (ASTRAL-MP "mechanism B").
 *
 * Mirrors ASTRAL-MP {@code WQDataCollection.addBipartitionsFromSignleIndTreesToX}
 * (lines 172–227), which for single-individual data runs over each INPUT gene tree
 * (there {@code allGreedies[gt] = [input tree]}).  For every polytomous node:
 *   1. Build the d arms = the k child subtrees + the complement (Lg \ sub(u)).
 *   2. {@code ROUNDS} times: sample one random taxon per arm.
 *   3. Restrict the reference {@code guide} tree (ASTRAL-MP's UPGMA species tree {@code ST})
 *      to those d reps and take its induced bipartitions ({@code Utils.getBitsets}).
 *   4. Expand each induced bipartition back to the union of its arms
 *      ({@code addbackAfterSampling}) and add that taxon set to X.
 *
 * Each emitted arm-union is a (generally non-contiguous) taxon set, represented as a
 * MULTI-RANGE {@link Cluster} anchored on the gene tree's own postorder: selected child
 * arms are contiguous child ranges, and the complement arm is the gene tree's two end
 * ranges {@code [0,ownLo)} and {@code [ownHi,L)}.  Its signature is computed from the
 * gene tree's prefix-hash arrays over those ranges, so it is directly comparable to the
 * gene-tree-derived cluster signatures already in X.
 *
 * This is X enrichment only (it proposes candidate bipartitions); it computes no quartet
 * signal — that is the d-partition QI weight (polytomy-design.md §3.8), which is separate.
 * It is gated behind {@code --resolve-input-gene-tree-polytomies} (default OFF) because it
 * enlarges X.
 */
public final class GeneTreePolytomySampler {

    /** Sampling rounds per polytomy — matches ASTRAL-MP's {@code for(ii<3)} loop. */
    public static final int ROUNDS = 3;
    /** int rep-bitmap cap (first implementation; larger polytomies are skipped). */
    public static final int MAX_ARMS = 31;

    private GeneTreePolytomySampler() {}

    /**
     * @param trees         exemplar tree list (cluster lookups index into it); arm-union
     *                      clusters are anchored on {@code trees[g]} for g &lt; numGeneTrees.
     * @param numGeneTrees  number of actual input gene trees (exclude appended guide/exemplars).
     * @param guide         reference resolution tree (UPGMA guide), spanning all taxa.
     * @param pref          prefix-hash arrays built over {@code trees} (for signatures).
     */
    public static int run(List<Tree> trees, int numGeneTrees, Tree guide,
                          PrefixHashArrays pref, int numTaxa, int m,
                          long seed, ClusterTable clusterTable) {
        if (guide == null || guide.root == null) {
            Logging.info("Gene-tree polytomy enrichment: no usable guide tree — skipped");
            return 0;
        }
        Random rng = new Random(seed);
        int[] stats = {0, 0, 0};   // {added, polytomies, skippedBig}
        int[] armAtPos = new int[guide.leafCount];   // reused scratch
        for (int g = 0; g < numGeneTrees; g++) {
            Tree gt = trees.get(g);
            if (gt.root == null) continue;
            walk(gt.root, gt, g, guide, armAtPos, pref, numTaxa, m, rng, clusterTable, stats);
        }
        Logging.info("Gene-tree polytomy enrichment: %d polytomies sampled (%d skipped: >%d arms), "
            + "%d arm-union clusters added to X", stats[1], stats[2], MAX_ARMS, stats[0]);
        return stats[0];
    }

    private static void walk(TreeNode node, Tree gt, int g, Tree guide, int[] armAtPos,
                             PrefixHashArrays pref, int numTaxa, int m, Random rng,
                             ClusterTable clusterTable, int[] stats) {
        if (node.isLeaf()) return;
        if (node.isPolytomous()) {
            for (TreeNode c : node.children)
                walk(c, gt, g, guide, armAtPos, pref, numTaxa, m, rng, clusterTable, stats);
            samplePolytomy(node, gt, g, guide, armAtPos, pref, numTaxa, m, rng, clusterTable, stats);
        } else {
            walk(node.left,  gt, g, guide, armAtPos, pref, numTaxa, m, rng, clusterTable, stats);
            walk(node.right, gt, g, guide, armAtPos, pref, numTaxa, m, rng, clusterTable, stats);
        }
    }

    private static void samplePolytomy(TreeNode u, Tree gt, int g, Tree guide, int[] armAtPos,
                                       PrefixHashArrays pref, int numTaxa, int m, Random rng,
                                       ClusterTable clusterTable, int[] stats) {
        int k = u.children.length;
        int L = gt.leafCount;
        int ownLo = u.rangeStart, ownHi = u.rangeEnd;
        boolean hasComp = (ownHi - ownLo) < L;       // complement arm non-empty
        int d = k + (hasComp ? 1 : 0);
        stats[1]++;
        if (d < 4) return;                            // need ≥4 arms for a non-trivial induced split
        if (d > MAX_ARMS) { stats[2]++; return; }

        int[] armLo = new int[k], armHi = new int[k];
        for (int i = 0; i < k; i++) { armLo[i] = u.children[i].rangeStart; armHi[i] = u.children[i].rangeEnd; }
        int compSize = L - (ownHi - ownLo);

        for (int round = 0; round < ROUNDS; round++) {
            // sample one rep taxon per arm
            int[] rep = new int[d];
            for (int i = 0; i < k; i++)
                rep[i] = gt.postorderArray[armLo[i] + rng.nextInt(armHi[i] - armLo[i])];
            if (hasComp) {
                int idx = rng.nextInt(compSize);              // index into [0,ownLo) ∪ [ownHi,L)
                int pos = (idx < ownLo) ? idx : ownHi + (idx - ownLo);
                rep[k] = gt.postorderArray[pos];
            }

            // restrict the guide tree to the d reps → induced clades (arm bitmaps)
            java.util.Arrays.fill(armAtPos, -1);
            int present = 0;
            for (int i = 0; i < d; i++) {
                int p = guide.positionMap[rep[i]];
                if (p >= 0) { armAtPos[p] = i; present++; }
            }
            if (present < 4) continue;
            ArrayList<Integer> clades = new ArrayList<>();
            walkGuide(guide.root, true, armAtPos, d, clades);

            for (int bm : clades)
                if (emitArmUnion(bm, k, hasComp, armLo, armHi, ownLo, ownHi, gt, g,
                                 pref, numTaxa, m, clusterTable)) stats[0]++;
        }
    }

    /** Induced-clade enumeration on the binary guide tree (subset-of-arms bitmaps). */
    private static int walkGuide(TreeNode node, boolean isRoot, int[] armAtPos, int d,
                                 ArrayList<Integer> out) {
        if (node.isLeaf()) {
            int arm = armAtPos[node.rangeStart];
            return arm >= 0 ? (1 << arm) : 0;
        }
        int lb = walkGuide(node.left,  false, armAtPos, d, out);
        int rb = walkGuide(node.right, false, armAtPos, d, out);
        int bm = lb | rb;
        if (isRoot) return bm;
        int legit = (lb != 0 ? 1 : 0) + (rb != 0 ? 1 : 0);
        if (legit < 2) return bm;                     // not a branching point of the induced tree
        int sz = Integer.bitCount(bm);
        if (sz >= 2 && sz <= d - 1) out.add(bm);      // emit this induced clade
        return bm;
    }

    /** Expand an arm-bitmap to its gene-tree multi-range, compute the signature, add to X. */
    private static boolean emitArmUnion(int armBitmap, int k, boolean hasComp,
                                        int[] armLo, int[] armHi, int ownLo, int ownHi,
                                        Tree gt, int g, PrefixHashArrays pref,
                                        int numTaxa, int m, ClusterTable clusterTable) {
        List<int[]> ranges = new ArrayList<>();
        int size = 0;
        for (int i = 0; i < k; i++) {
            if ((armBitmap & (1 << i)) != 0) {
                ranges.add(new int[]{armLo[i], armHi[i]});
                size += armHi[i] - armLo[i];
            }
        }
        if (hasComp && (armBitmap & (1 << k)) != 0) {
            if (ownLo > 0)            { ranges.add(new int[]{0, ownLo});          size += ownLo; }
            if (ownHi < gt.leafCount) { ranges.add(new int[]{ownHi, gt.leafCount}); size += gt.leafCount - ownHi; }
        }
        if (size < 2 || size > numTaxa - 2) return false;
        if (ranges.size() > 1) ranges.sort((x, y) -> Integer.compare(x[0], y[0]));

        int[] los = new int[ranges.size()], his = new int[ranges.size()];
        for (int j = 0; j < ranges.size(); j++) { los[j] = ranges.get(j)[0]; his[j] = ranges.get(j)[1]; }

        long[] sums = new long[m], xors = new long[m];
        for (int s = 0; s < m; s++) {
            long su = 0, xo = 0;
            for (int j = 0; j < los.length; j++) {
                su += pref.rangeSum(g, s, los[j], his[j]);
                xo ^= pref.rangeXor(g, s, los[j], his[j]);
            }
            sums[s] = su; xors[s] = xo;
        }
        ClusterHash sig = new ClusterHash(sums, xors, size, m);
        Cluster c = new Cluster(g, los, his, /*complement=*/false, size);
        return clusterTable.addCluster(sig, c);
    }
}
