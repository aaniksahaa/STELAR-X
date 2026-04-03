package hash;

/**
 * Computes m independent 64-bit sparsified hashes for each taxon.
 *
 * Uses SplitMix64 finalizer:  h(taxonId, seed) = mix64(taxonId + seed)
 * Seeds are derived from a base seed so the set is reproducible.
 *
 * The key property: outputs are uniformly distributed over [0, 2^64),
 * so sums/XORs of distinct taxon hashes have good collision resistance.
 */
public class TaxonHasher {

    private final int m;          // number of seeds
    private final long[] seeds;   // seeds[s] for s = 0..m-1

    /**
     * hashes[s][taxonId] = mix64(taxonId + seeds[s])
     * Shape: (m, n)
     */
    private final long[][] hashes;

    public TaxonHasher(int numTaxa, int numSeeds, long baseSeed) {
        this.m = numSeeds;
        this.seeds = new long[m];
        // Derive m seeds from baseSeed using the same mix function
        for (int s = 0; s < m; s++) {
            seeds[s] = mix64(baseSeed + s * 0x9E3779B97F4A7C15L);
        }

        hashes = new long[m][numTaxa];
        for (int s = 0; s < m; s++) {
            for (int tid = 0; tid < numTaxa; tid++) {
                hashes[s][tid] = mix64((long) tid + seeds[s]);
                // Guarantee non-zero (the all-zero hash is reserved for "absent")
                if (hashes[s][tid] == 0L) hashes[s][tid] = 1L;
            }
        }
    }

    /** Return hash for taxon id under seed index s. */
    public long get(int seedIdx, int taxonId) {
        return hashes[seedIdx][taxonId];
    }

    /** Number of seeds. */
    public int numSeeds() { return m; }

    /** Number of taxa. */
    public int numTaxa() { return hashes[0].length; }

    // -------------------------------------------------------------------------
    // SplitMix64 finalizer  (standard bijection on 64-bit integers)
    // -------------------------------------------------------------------------
    public static long mix64(long x) {
        x ^= (x >>> 30);
        x *= 0xBF58476D1CE4E5B9L;
        x ^= (x >>> 27);
        x *= 0x94D049BB133111EBL;
        x ^= (x >>> 31);
        return x;
    }
}
