package stelarx.greedy;

/**
 * {@code long[]}-bitmap variant of {@link MiniGreedyBuilder}, used by Step B's
 * {@code resolveLinearly} for LARGE polytomies ({@code d > 31}) under
 * {@code --stepb-process-large-polytomies}.  Logic is byte-for-byte the same
 * as {@link MiniGreedyBuilder} — same {@code buildTreeFromClusters} accept rule
 * (LCA classification, ≥ 2 children moved, not redundant with the LCA) — only
 * the rep-subset representation changes from a single {@code int} to a
 * {@code long[]} word array of length {@code W = ceil(d/64)}.
 *
 * The small ({@code d ≤ 31}) path stays on {@link MiniGreedyBuilder} so existing
 * behaviour is unchanged.
 */
final class MiniGreedyBuilderLong {

    private final int d;
    private final int W;                 // words per bitmap
    private final long[] allBits;

    private final int[]    parent;
    private final long[][] bitmap;       // bits in subtree (W longs each)
    private int[][]        children;
    private final int[]    childCount;
    private final int      virtualRoot;
    private final int[]    repToLeafNode;
    private int            nextId;

    MiniGreedyBuilderLong(int d) {
        if (d < 1) throw new IllegalArgumentException("MiniGreedyBuilderLong: d ≥ 1");
        this.d = d;
        this.W = (d + 63) >>> 6;
        this.allBits = new long[W];
        for (int b = 0; b < d; b++) allBits[b >>> 6] |= (1L << (b & 63));

        int cap = 2 * d + 2;
        this.parent     = new int[cap];
        this.bitmap     = new long[cap][];
        this.children   = new int[cap][];
        this.childCount = new int[cap];

        this.virtualRoot = 0;
        parent[0]     = -1;
        bitmap[0]     = allBits.clone();
        children[0]   = new int[d];
        childCount[0] = d;
        for (int r = 0; r < d; r++) {
            int leaf = 1 + r;
            parent[leaf]   = virtualRoot;
            long[] bm = new long[W];
            bm[r >>> 6] = (1L << (r & 63));
            bitmap[leaf]   = bm;
            children[leaf] = new int[0];
            childCount[leaf] = 0;
            children[0][r] = leaf;
        }
        this.repToLeafNode = new int[d];
        for (int r = 0; r < d; r++) repToLeafNode[r] = 1 + r;
        this.nextId = 1 + d;
    }

    /** Attempt to insert cluster bitmap {@code bm}; true if accepted. */
    boolean tryInsert(long[] bm) {
        int sz = popcount(bm);
        if (sz < 2 || sz > d - 1) return false;

        int lca = findLCA(bm);
        int kc  = childCount[lca];
        int[] kids = children[lca];

        int[] movedIds = new int[kc];
        int movedN = 0, accounted = 0;

        for (int i = 0; i < kc; i++) {
            int c = kids[i];
            long[] cbm = bitmap[c];
            boolean anyAnd = false, subsetOfBm = true;
            for (int k = 0; k < W; k++) {
                long an = cbm[k] & bm[k];
                if (an != 0)        anyAnd = true;
                if (an != cbm[k])   subsetOfBm = false;
            }
            if (!anyAnd) continue;            // disjoint sibling
            if (subsetOfBm) {                 // child wholly inside bm → move
                movedIds[movedN++] = c;
                accounted += popcount(cbm);
                continue;
            }
            return false;                     // cross-cut
        }
        if (accounted != sz) return false;
        if (movedN < 2)      return false;
        if (movedN == kc)    return false;

        int nNode = nextId++;
        parent[nNode]   = lca;
        bitmap[nNode]   = bm.clone();
        int[] newKids   = new int[movedN];
        System.arraycopy(movedIds, 0, newKids, 0, movedN);
        children[nNode]   = newKids;
        childCount[nNode] = movedN;
        for (int j = 0; j < movedN; j++) parent[movedIds[j]] = nNode;

        int newLcaKc = kc - movedN + 1;
        int[] newLcaKids = new int[newLcaKc];
        int p = 0;
        outer:
        for (int i = 0; i < kc; i++) {
            int c = kids[i];
            for (int j = 0; j < movedN; j++) if (movedIds[j] == c) continue outer;
            newLcaKids[p++] = c;
        }
        newLcaKids[p++] = nNode;
        children[lca]   = newLcaKids;
        childCount[lca] = p;
        return true;
    }

    /** Smallest currently-accepted node whose subtree bitmap ⊇ bm. */
    private int findLCA(long[] bm) {
        int anyRep = lowestSetBit(bm);
        int node = repToLeafNode[anyRep];
        while (true) {
            if (subset(bm, bitmap[node])) return node;   // bm ⊆ node
            int pr = parent[node];
            if (pr < 0) return virtualRoot;
            node = pr;
        }
    }

    /** Iterate every accepted internal node (excluding root + d leaves). */
    void forEachAcceptedInternal(java.util.function.Consumer<long[]> visitor) {
        for (int id = 1 + d; id < nextId; id++) visitor.accept(bitmap[id]);
    }

    /** ASTRAL-MP resolveLinearly leftover step (long[] variant of
     *  {@link MiniGreedyBuilder#resolveLeftoverPolytomiesRandomly}). */
    void resolveLeftoverPolytomiesRandomly(java.util.Random rng,
                                           java.util.function.Consumer<long[]> emit) {
        for (int id = 0; id < nextId; id++) {
            if (childCount[id] < 3) continue;
            java.util.ArrayList<long[]> parts = new java.util.ArrayList<>(childCount[id] + 1);
            for (int j = 0; j < childCount[id]; j++) parts.add(bitmap[children[id][j]].clone());
            long[] rest = new long[W];
            boolean restNonEmpty = false;
            for (int k = 0; k < W; k++) {
                rest[k] = allBits[k] & ~bitmap[id][k];
                if (rest[k] != 0) restNonEmpty = true;
            }
            if (restNonEmpty) parts.add(rest);
            while (parts.size() > 2) {
                long[] c1 = parts.remove(rng.nextInt(parts.size()));
                long[] c2 = parts.remove(rng.nextInt(parts.size()));
                long[] merged = new long[W];
                for (int k = 0; k < W; k++) merged[k] = c1[k] | c2[k];
                emit.accept(merged);
                parts.add(merged);
            }
        }
    }

    // ── bitmap helpers ─────────────────────────────────────────────────────
    private int popcount(long[] bm) {
        int c = 0;
        for (int k = 0; k < W; k++) c += Long.bitCount(bm[k]);
        return c;
    }
    private int lowestSetBit(long[] bm) {
        for (int k = 0; k < W; k++) if (bm[k] != 0)
            return (k << 6) + Long.numberOfTrailingZeros(bm[k]);
        return 0;
    }
    /** true iff a ⊆ b. */
    private boolean subset(long[] a, long[] b) {
        for (int k = 0; k < W; k++) if ((a[k] & ~b[k]) != 0) return false;
        return true;
    }
}
