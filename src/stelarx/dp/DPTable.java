package stelarx.dp;

import stelarx.Config;
import stelarx.Logging;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.gpu.GPUDPBuilder;
import stelarx.hash.PrefixHashArrays;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.ProgressBar;
import stelarx.util.Threading;

import java.util.*;

/**
 * DP search space: maps each cluster hash to its set of candidate bipartition splits.
 *
 * Built via Mode 1 (tree-local transitions only), O(nk):
 *
 *   Type 1  -- for every resolved binary internal node u (incl. root):
 *              sub(u) → sub(left(u)) | sub(right(u))
 *
 *   Type 2  -- for every resolved non-root node u (leaf or internal) with a
 *              binary parent and a nonempty parent super-complement:
 *              [S \ sub(u)] → sub(sibling(u)) | [S \ sub(parent(u))]
 *
 * For complete trees (Lg == S), Types 3a/3b add nothing new and are skipped.
 *
 * The root of the DP is the all-taxa cluster (from ClusterTable.getAllTaxaHash()).
 * Its Type 1 transition(s) are stored at that hash key.
 */
public class DPTable {

    // transitions[parentHash] = set of distinct BipartitionSplits for that cluster
    private final Map<ClusterHash, Set<BipartitionSplit>> transitions = new HashMap<>();

    private final ClusterHash rootHash; // allTaxaHash -- starting point of the DP
    private final int m;
    private final int n; // total taxa count

    // stats
    private int totalEmitted = 0;   // total emitted (with duplicates across trees)
    private int uniqueSplits = 0;   // sum of set sizes after dedup

    // -------------------------------------------------------------------------

    // Anchored-outgroup mode: skip building the with-anchor-parent
    // transitions — they are unreachable once the root is anchored, so omitting them
    // (and the redundant per-cluster root-split search) is exact and saves memory.
    private final boolean anchorFreeX;
    private final int     anchor;

    public DPTable(List<Tree> trees, PrefixHashArrays pref, ClusterTable clusterTable) {
        long t0 = System.nanoTime();
        this.m        = pref.numSeeds();
        this.rootHash = clusterTable.getAllTaxaHash();
        this.n        = rootHash.size;
        Config cfg = Config.getInstance();
        this.anchorFreeX = cfg.isAnchorOutgroup();
        this.anchor = anchorFreeX ? cfg.getAnchorTaxon() : -1;

        int treesDone = 0;
        ProgressBar localBar = new ProgressBar("Local DP transitions", trees.size());
        for (Tree tree : trees) {
            extractFromTree(tree, pref);
            localBar.update(++treesDone);
        }
        localBar.done();

        // Count total unique splits
        for (Set<BipartitionSplit> s : transitions.values()) uniqueSplits += s.size();

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("DP table: %d clusters with splits, %d unique splits (%d emitted) in %d ms",
            transitions.size(), uniqueSplits, totalEmitted, ms);
    }

    // -------------------------------------------------------------------------
    // Tree traversal
    // -------------------------------------------------------------------------

    private void extractFromTree(Tree tree, PrefixHashArrays pref) {
        emit(tree.root, tree.treeIndex, -1, pref);
    }

    /** Post-order recursion: emit transitions for this node, then children. */
    private void emit(TreeNode u, int ti, int anchorPos, PrefixHashArrays pref) {
        // A leaf has no Type-1 subtree split, but it can still induce a Type-2
        // split of its complement.  If u={x}, sibling(u)=B, and the taxa outside
        // parent(u) are O, the valid unrooted rotation is
        //
        //     S\{x} = B ∪ O  →  B | O.
        //
        // Keep this separate from the internal-node traversal so leaves remain
        // ordinary DP base cases while their incident edge can contribute a
        // resolution for the opposite side of the bipartition.
        if (u.isLeaf()) return;

        // Polytomous node: recurse into all children, but add NO direct transitions of
        // its own.  A polytomy is an unresolved node — it must not force any binary
        // resolution into the DP search space; its quartet signal still enters via the
        // d-partition QI weight.  (polytomy-design.md §3.7.)
        if (u.isPolytomous()) {
            for (TreeNode child : u.children) emit(child, ti, anchorPos, pref);
            return;
        }

        emit(u.left,  ti, anchorPos, pref);
        emit(u.right, ti, anchorPos, pref);

        // In anchor-free mode, sub(u) contains the anchor iff the anchor's position in
        // this tree falls in [rangeStart,rangeEnd).  Exactly one of the two parents
        // below (sub(u) for Type 1, S\sub(u) for Type 2) is then anchor-free; the other
        // is unreachable from the anchored root, so we skip building it.
        ClusterHash hU     = hashRange(ti, u.rangeStart,       u.rangeEnd,       false, pref);
        ClusterHash hLeft  = hashRange(ti, u.left.rangeStart,  u.left.rangeEnd,  false, pref);
        ClusterHash hRight = hashRange(ti, u.right.rangeStart, u.right.rangeEnd, false, pref);
        addTransition(hU, hLeft, hRight);
    }

    /**
     * Emit the complement-side rotation induced by {@code u}. Unlike Type 1,
     * this is valid for leaves as well as binary internal nodes: it splits the
     * complement of {@code sub(u)}, not {@code sub(u)} itself.
     */
    private void emitType2(TreeNode u, int ti, boolean subHasAnchor,
                           PrefixHashArrays pref) {
        // For non-root u: if parent is root and the tree is complete,
        // S\sub(root)=empty, so hCompParent.size==0 and we skip. For incomplete
        // trees, S\sub(root)=S\Lg is nonempty and the transition is valid.
        //
        // GUARD: u must not be polytomous and its parent must be binary.
        // A child of a polytomous parent has no unique sibling, while a
        // polytomous u must not inject a binary resolution of its own.
        //
        // Anchor-free: S\sub(u) is anchor-free iff the anchor IS in sub(u).
        if (u.isRoot() || u.isPolytomous() || u.parent.isPolytomous()
                || (anchorFreeX && !subHasAnchor)) {
            return;
        }

        TreeNode sib    = u.getSibling();
        TreeNode parent = u.parent;

        ClusterHash hCompU      = hashRange(ti, u.rangeStart,      u.rangeEnd,      true,  pref);
        ClusterHash hSib        = hashRange(ti, sib.rangeStart,    sib.rangeEnd,    false, pref);
        ClusterHash hCompParent = hashRange(ti, parent.rangeStart, parent.rangeEnd, true,  pref);
        if (hCompParent.size > 0) {
            addTransition(hCompU, hSib, hCompParent);
        }
    }

    // -------------------------------------------------------------------------

    private void addTransition(ClusterHash parent, ClusterHash a, ClusterHash b) {
        totalEmitted++;
        BipartitionSplit split = new BipartitionSplit(a, b);
        transitions.computeIfAbsent(parent, k -> new LinkedHashSet<>()).add(split);
    }

    /**
     * Compute a finalized ClusterHash for the range [lo,hi) in tree ti.
     * complement=true gives the super-complement S\[lo,hi) (w.r.t. ALL n taxa).
     */
    private ClusterHash hashRange(int ti, int lo, int hi, boolean complement,
                                   PrefixHashArrays pref) {
        long[] rawSums = new long[m], rawXors = new long[m];
        for (int s = 0; s < m; s++) {
            rawSums[s] = complement ? pref.superCompSum(ti, s, lo, hi) : pref.rangeSum(ti, s, lo, hi);
            rawXors[s] = complement ? pref.superCompXor(ti, s, lo, hi) : pref.rangeXor(ti, s, lo, hi);
        }
        int sz = complement ? (n - (hi - lo)) : (hi - lo);
        return new ClusterHash(rawSums, rawXors, sz, m);
    }

    // -------------------------------------------------------------------------
    // Mode 2: Cross-tree DP transitions
    // -------------------------------------------------------------------------

    /**
     * Expand the search space with cross-tree splits (ASTRAL Mode 2 / "full" search).
     *
     * For every cluster A ∈ X and every cluster B ∈ X with |B| ≤ |A|/2:
     *   if hash(A) − hash(B) matches another cluster R ∈ X  →  add A → B | R.
     *
     * Also handles the all-taxa root cluster (not in X, but its transitions matter).
     *
     * @param clusterTable  the cluster set X
     * @param useGPU        true to use CUDA acceleration; false for parallel CPU
     */
    public void addCrossTreeTransitions(ClusterTable clusterTable, boolean useGPU) {
        long t0 = System.nanoTime();
        int beforeSplits = 0;
        for (Set<BipartitionSplit> s : transitions.values()) beforeSplits += s.size();

        if (useGPU) {
            addCrossTreeGPU(clusterTable);
        } else {
            addCrossTreeCPU(clusterTable);
        }

        // Recount
        uniqueSplits = 0;
        for (Set<BipartitionSplit> s : transitions.values()) uniqueSplits += s.size();

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Cross-tree transitions (Mode 2): +%d splits (%d total) in %d ms",
            uniqueSplits - beforeSplits, uniqueSplits, ms);
    }

    // ── CPU path ─────────────────────────────────────────────────────────────

    private void addCrossTreeCPU(ClusterTable clusterTable) {
        List<ClusterHash> allHashes = new ArrayList<>();
        for (var e : clusterTable.entries()) allHashes.add(e.hash);
        int N = allHashes.size();

        // Parallel over all clusters A: each thread writes to its own perCluster[idx] slot
        @SuppressWarnings("unchecked")
        Set<BipartitionSplit>[] perCluster = new Set[N];

        java.util.concurrent.atomic.AtomicInteger cpuDone = new java.util.concurrent.atomic.AtomicInteger(0);
        ProgressBar cpuBar = new ProgressBar("Cross-tree DP (CPU)", N);
        ThreadLocal<ClusterTable.ResidualLookup> residualLookups =
            ThreadLocal.withInitial(clusterTable::newResidualLookup);
        Threading.processRangeParallel(N, idx -> {
            ClusterHash hashA = allHashes.get(idx);
            int szA = hashA.size;
            // Every A is processed by exactly one worker. Existing transition
            // sets are distinct objects, so workers can extend them directly
            // without duplicating all splits in a temporary perCluster set.
            Set<BipartitionSplit> localSet = transitions.get(hashA);
            boolean detached = false;
            ClusterTable.ResidualLookup residualLookup = residualLookups.get();

            for (int sz = 1; sz <= szA / 2; sz++) {
                for (ClusterHash hashB : clusterTable.getBySize(sz)) {
                    ClusterHash residual = residualLookup.find(hashA, hashB);
                    if (residual != null) {
                        if (localSet == null) {
                            localSet = new LinkedHashSet<>();
                            detached = true;
                        }
                        localSet.add(new BipartitionSplit(hashB, residual));
                    }
                }
            }

            if (detached) perCluster[idx] = localSet;
            cpuBar.update(cpuDone.incrementAndGet());
        });
        cpuBar.done();

        // Serial merge into transitions (different A → different keys, no map contention)
        for (int idx = 0; idx < N; idx++) {
            if (perCluster[idx] != null) {
                ClusterHash hashA = allHashes.get(idx);
                transitions.computeIfAbsent(hashA, k -> new LinkedHashSet<>())
                           .addAll(perCluster[idx]);
            }
        }

        // Also handle the root (all-taxa) cluster — not in clusterTable but is the DP root
        searchRootTransitions(clusterTable);
    }

    // ── GPU path ─────────────────────────────────────────────────────────────

    private void addCrossTreeGPU(ClusterTable clusterTable) {
        List<ClusterTable.Entry> entries = new ArrayList<>(clusterTable.entries());
        int N = entries.size();
        if (N == 0) { searchRootTransitions(clusterTable); return; }

        int maxSize = clusterTable.sizes().stream().mapToInt(Integer::intValue).max().orElse(1);

        // ── Flatten cluster data ──────────────────────────────────────────────
        long[] clusterSums  = new long[N * m];
        long[] clusterXors  = new long[N * m];
        int[]  clusterSizes = new int[N];

        for (int c = 0; c < N; c++) {
            ClusterHash h = entries.get(c).hash;
            clusterSizes[c] = h.size;
            for (int s = 0; s < m; s++) {
                clusterSums[c * m + s] = h.sums[s];
                clusterXors[c * m + s] = h.xors[s];
            }
        }

        // ── Build sortedBySize and binStart ───────────────────────────────────
        // sortedBySize[i] = cluster index (into above arrays) ordered by size asc
        // binStart[sz]    = first index in sortedBySize with size >= sz
        Integer[] order = new Integer[N];
        for (int i = 0; i < N; i++) order[i] = i;
        Arrays.sort(order, Comparator.comparingInt(i -> clusterSizes[i]));
        int[] sortedBySize = new int[N];
        for (int i = 0; i < N; i++) sortedBySize[i] = order[i];

        int[] binStart = new int[maxSize + 2];
        int ptr = 0;
        for (int sz = 0; sz <= maxSize + 1; sz++) {
            while (ptr < N && clusterSizes[sortedBySize[ptr]] < sz) ptr++;
            binStart[sz] = ptr;
        }

        // ── Compute maxPerRound to bound GPU output buffer ────────────────────
        // Default 120 MB = 10M triples.  Configurable via --gpu-dp-state-space-construction-output-cap.
        // Sub-batching within each round normally ensures this is never exceeded.
        int    maxPerRound    = Config.getInstance().getGpuDpOutputCapTriples();
        double progressInterval  = Config.getInstance().getGpuDpProgressInterval();
        int    progressMaxSteps  = Config.getInstance().getGpuDpProgressMaxSteps();

        // ── Call GPU ──────────────────────────────────────────────────────────
        Logging.debug("  GPU cross-tree search: N=%d clusters, maxSize=%d", N, maxSize);
        int[] raw = GPUDPBuilder.findCrossTreeTransitionsGPU(
            clusterSums, clusterXors, clusterSizes,
            N, m,
            sortedBySize, binStart, maxSize,
            maxPerRound, progressInterval, progressMaxSteps);

        if (raw == null) {
            Logging.info("  GPU cross-tree search returned null, falling back to CPU");
            addCrossTreeCPU(clusterTable);
            return;
        }

        // ── Process GPU results ───────────────────────────────────────────────
        int count = raw[0];
        Logging.debug("  GPU cross-tree: %d raw pairs found", count);
        for (int i = 0; i < count; i++) {
            int idxA   = raw[1 + i * 3];
            int idxB   = raw[1 + i * 3 + 1];
            int idxRes = raw[1 + i * 3 + 2];
            ClusterHash hashA   = entries.get(idxA).hash;
            ClusterHash hashB   = entries.get(idxB).hash;
            ClusterHash hashRes = entries.get(idxRes).hash;
            transitions.computeIfAbsent(hashA, k -> new LinkedHashSet<>())
                       .add(new BipartitionSplit(hashB, hashRes));
        }

        // Root cluster transitions (handled on CPU, fast)
        searchRootTransitions(clusterTable);
    }

    /** Search transitions for the all-taxa root cluster (not in X itself). */
    private void searchRootTransitions(ClusterTable clusterTable) {
        // Anchor-free mode replaces the root's transitions with a single anchored split
        // (applyAnchoredRoot), so building all complementary root splits here is wasted.
        if (anchorFreeX) return;
        int szRoot = rootHash.size;
        Set<BipartitionSplit> rootSet =
            transitions.computeIfAbsent(rootHash, k -> new LinkedHashSet<>());
        ClusterTable.ResidualLookup residualLookup = clusterTable.newResidualLookup();
        for (int sz = 1; sz <= szRoot / 2; sz++) {
            for (ClusterHash hashB : clusterTable.getBySize(sz)) {
                ClusterHash residual = residualLookup.find(rootHash, hashB);
                if (residual != null) {
                    rootSet.add(new BipartitionSplit(hashB, residual));
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Queries
    // -------------------------------------------------------------------------

    /**
     * Anchored-outgroup root: replace the all-taxa root's ENTIRE transition set with
     * the single split {@code ({anchor} | S\{anchor})}.
     *
     * Legacy unrooted implementation retained for internal binary/JNI compatibility.
     * STELAR-X rejects this mode before construction, so this method is unreachable
     * from the public CLI.
     *
     * Must be called AFTER all local (Mode 1) and cross-tree (Mode 2) transitions are
     * built. A no-op (with a warning) if the anchor hash is null.
     */
    public void applyAnchoredRoot(ClusterHash anchorHash) {
        if (anchorHash == null) {
            Logging.info("Anchored root requested but anchor singleton hash is unavailable "
                + "(anchor taxon in no tree?) — root transitions left unchanged");
            return;
        }
        ClusterHash sAnchor = ClusterHash.residual(rootHash, anchorHash);
        Set<BipartitionSplit> old = transitions.get(rootHash);
        int removed = (old != null) ? old.size() : 0;

        Set<BipartitionSplit> anchored = new LinkedHashSet<>();
        anchored.add(new BipartitionSplit(anchorHash, sAnchor));
        transitions.put(rootHash, anchored);

        uniqueSplits = 0;
        for (Set<BipartitionSplit> s : transitions.values()) uniqueSplits += s.size();

        Logging.info("Anchored root (exact, unrooted-invariant): replaced %d root split(s) with 1 "
            + "({anchor}=%d | S\\{anchor}=%d taxa); %d total unique splits",
            removed, anchorHash.size, sAnchor.size, uniqueSplits);
    }

    /**
     * Clusters actually reachable from the root by the top-down inference DP
     * (see {@link stelarx.dp.Inference}#solve): the root is reachable, and both
     * halves of every split of a reachable cluster are reachable.
     *
     * The DP calls {@code weightTable.getScore(split)} ONLY for splits of clusters
     * it reaches, so a cluster never reached contributes nothing to the root score
     * — its splits need no weight and can be skipped without changing the result.
     * This mirrors {@code solve}'s recursion exactly (minus the arithmetic), so the
     * returned set is a superset of every cluster whose splits the DP actually
     * queries; filtering the weight step to these clusters is therefore
     * result-preserving.
     *
     * O(reachable splits) time; the set holds only {@link ClusterHash} references
     * already owned by the transitions map — no splits or clusters are copied.
     */
    public Set<ClusterHash> reachableClusters() {
        Set<ClusterHash> reachable = new HashSet<>();
        ArrayDeque<ClusterHash> queue = new ArrayDeque<>();
        reachable.add(rootHash);
        queue.add(rootHash);
        while (!queue.isEmpty()) {
            ClusterHash ch = queue.poll();
            for (BipartitionSplit sp : getSplits(ch)) {
                if (reachable.add(sp.lo)) queue.add(sp.lo);
                if (reachable.add(sp.hi)) queue.add(sp.hi);
            }
        }
        return reachable;
    }

    public ClusterHash getRootHash()                       { return rootHash; }
    public Set<BipartitionSplit> getSplits(ClusterHash h)  { return transitions.getOrDefault(h, Collections.emptySet()); }
    public boolean hasSplits(ClusterHash h)                { return transitions.containsKey(h); }
    public int numClusters()                               { return transitions.size(); }
    public int numUniqueSplits()                           { return uniqueSplits; }
    public int numEmitted()                                { return totalEmitted; }
    public boolean isAnchorFree()                          { return anchorFreeX; }
    public int getAnchorTaxon()                            { return anchor; }
    public Set<Map.Entry<ClusterHash, Set<BipartitionSplit>>> entries() { return transitions.entrySet(); }
}
