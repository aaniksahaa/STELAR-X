package stelarx.greedy;

import stelarx.cluster.ClusterHash;

import java.util.Collection;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Thread-safe accumulator for polytomy-resolution emissions.
 *
 * Keyed by {@link ClusterHash} signature. In parallel polytomy resolution each
 * task owns a private instance, so adaptive-round novelty cannot be contaminated
 * by unrelated tasks. Completed task buffers merge immediately into the caller's
 * shared instance. If two tasks emit the same signature, a stable descriptor
 * ordering—not arrival time—selects the retained provenance.
 *
 * Phase 5 integration of these emissions into the global ClusterTable (with
 * exemplars, either by gene-tree lookup or by synthesizing multi-range
 * exemplars) happens in a separate later pass; this buffer is the canonical
 * record of WHAT was emitted by Part II.
 */
public final class EmissionBuffer {
    private final ConcurrentHashMap<ClusterHash, EmittedBipartition> emitted =
        new ConcurrentHashMap<>();

    /**
     * Add an emission. Returns true iff this signature was absent. On a duplicate,
     * atomically retain the deterministic preferred descriptor.
     */
    public boolean add(EmittedBipartition b) {
        EmittedBipartition current = emitted.putIfAbsent(b.signature, b);
        if (current == null) return true;
        while (prefer(b, current) == b) {
            if (emitted.replace(b.signature, current, b)) break;
            current = emitted.get(b.signature);
        }
        return false;
    }

    private static EmittedBipartition prefer(EmittedBipartition a,
                                              EmittedBipartition b) {
        int c = Character.compare(a.source, b.source); // Step A before Step B
        if (c != 0) return c < 0 ? a : b;
        c = Integer.compare(a.thresholdIndex, b.thresholdIndex);
        if (c != 0) return c < 0 ? a : b;
        int[] alo = a.canonicalSide.los, ahi = a.canonicalSide.his;
        int[] blo = b.canonicalSide.los, bhi = b.canonicalSide.his;
        c = Integer.compare(alo.length, blo.length);
        if (c != 0) return c < 0 ? a : b;
        for (int i = 0; i < alo.length; i++) {
            c = Integer.compare(alo[i], blo[i]);
            if (c != 0) return c < 0 ? a : b;
            c = Integer.compare(ahi[i], bhi[i]);
            if (c != 0) return c < 0 ? a : b;
        }
        return b;
    }

    public boolean contains(ClusterHash sig) { return emitted.containsKey(sig); }
    public int size()                        { return emitted.size(); }
    public Collection<EmittedBipartition> all() { return emitted.values(); }
}
