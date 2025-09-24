package tree;

import java.util.List;

/**
 * Efficient hash calculator for range-based bipartitions.
 * Uses prefix sums to compute range hashes in O(1) time after O(n) preprocessing.
 * Supports multiple hash functions with different parameters for collision resistance.
 */
public class RangeHashCalculator {
    
    // Configuration for different hash functions
    public static class HashConfig {
        public final long multiplier;
        public final long modulus;
        public final long base;
        
        public HashConfig(long multiplier, long modulus, long base) {
            this.multiplier = multiplier;
            this.modulus = modulus;
            this.base = base;
        }
        
        // Predefined hash configurations for good distribution
        public static final HashConfig HASH1 = new HashConfig(31L, 1000000007L, 1L);
        public static final HashConfig HASH2 = new HashConfig(37L, 1000000009L, 1L);
        public static final HashConfig ROLLING_HASH = new HashConfig(256L, 1000000021L, 1L);
    }
    
    private final int numGeneTrees;
    private final int[][] geneTreeOrderings;  // [geneTreeIndex][taxonPosition] = taxonId
    private final long[][] prefixHash1;       // Prefix sums for first hash function
    private final long[][] prefixHash2;       // Prefix sums for second hash function
    private final long[][] basePowers1;      // Powers of base for first hash
    private final long[][] basePowers2;      // Powers of base for second hash
    
    private final HashConfig config1;
    private final HashConfig config2;
    
    public RangeHashCalculator(List<Tree> geneTrees, int[][] geneTreeOrderings) {
        this(geneTrees, geneTreeOrderings, HashConfig.HASH1, HashConfig.HASH2);
    }
    
    public RangeHashCalculator(List<Tree> geneTrees, int[][] geneTreeOrderings, 
                              HashConfig config1, HashConfig config2) {
        this.numGeneTrees = geneTrees.size();
        this.geneTreeOrderings = geneTreeOrderings;
        this.config1 = config1;
        this.config2 = config2;
        
        // Initialize arrays
        int maxTaxa = geneTreeOrderings[0].length;
        this.prefixHash1 = new long[numGeneTrees][maxTaxa + 1];
        this.prefixHash2 = new long[numGeneTrees][maxTaxa + 1];
        this.basePowers1 = new long[numGeneTrees][maxTaxa + 1];
        this.basePowers2 = new long[numGeneTrees][maxTaxa + 1];
        
        // Precompute prefix sums and powers
        precomputeHashes();
    }
    
    private void precomputeHashes() {
        for (int treeIdx = 0; treeIdx < numGeneTrees; treeIdx++) {
            int[] ordering = geneTreeOrderings[treeIdx];
            int n = ordering.length;
            
            // Initialize base cases
            prefixHash1[treeIdx][0] = 0;
            prefixHash2[treeIdx][0] = 0;
            basePowers1[treeIdx][0] = config1.base;
            basePowers2[treeIdx][0] = config2.base;
            
            // Compute prefix sums and powers
            for (int i = 1; i <= n; i++) {
                int taxonId = ordering[i - 1];  // 0-indexed ordering, 1-indexed prefix
                
                // First hash function
                long contribution1 = (taxonId * basePowers1[treeIdx][i - 1]) % config1.modulus;
                prefixHash1[treeIdx][i] = (prefixHash1[treeIdx][i - 1] + contribution1) % config1.modulus;
                basePowers1[treeIdx][i] = (basePowers1[treeIdx][i - 1] * config1.multiplier) % config1.modulus;
                
                // Second hash function
                long contribution2 = (taxonId * basePowers2[treeIdx][i - 1]) % config2.modulus;
                prefixHash2[treeIdx][i] = (prefixHash2[treeIdx][i - 1] + contribution2) % config2.modulus;
                basePowers2[treeIdx][i] = (basePowers2[treeIdx][i - 1] * config2.multiplier) % config2.modulus;
            }
        }
    }
    
    /**
     * Computes range hash using first hash function in O(1) time.
     * Hash of range [start, end) in the specified gene tree.
     */
    public long getRangeHash1(int geneTreeIndex, int start, int end) {
        if (start >= end || geneTreeIndex >= numGeneTrees) {
            return 0;
        }
        
        long rangeHash = (prefixHash1[geneTreeIndex][end] - prefixHash1[geneTreeIndex][start] + config1.modulus) % config1.modulus;
        
        // Normalize by dividing by base^start to get position-independent hash
        if (start > 0) {
            long baseInverse = modularInverse(basePowers1[geneTreeIndex][start], config1.modulus);
            rangeHash = (rangeHash * baseInverse) % config1.modulus;
        }
        
        return rangeHash;
    }
    
    /**
     * Computes range hash using second hash function in O(1) time.
     */
    public long getRangeHash2(int geneTreeIndex, int start, int end) {
        if (start >= end || geneTreeIndex >= numGeneTrees) {
            return 0;
        }
        
        long rangeHash = (prefixHash2[geneTreeIndex][end] - prefixHash2[geneTreeIndex][start] + config2.modulus) % config2.modulus;
        
        // Normalize by dividing by base^start
        if (start > 0) {
            long baseInverse = modularInverse(basePowers2[geneTreeIndex][start], config2.modulus);
            rangeHash = (rangeHash * baseInverse) % config2.modulus;
        }
        
        return rangeHash;
    }
    
    /**
     * Alternative simple sum-based hash for quick comparison (less collision-resistant).
     */
    public long getSimpleSumHash(int geneTreeIndex, int start, int end) {
        if (start >= end || geneTreeIndex >= numGeneTrees) {
            return 0;
        }
        
        long sum = 0;
        int[] ordering = geneTreeOrderings[geneTreeIndex];
        for (int i = start; i < end; i++) {
            sum += ordering[i];
        }
        return sum;
    }
    
    /**
     * Computes modular multiplicative inverse using extended Euclidean algorithm.
     */
    private long modularInverse(long a, long m) {
        if (a == 1) return 1;
        
        long m0 = m, x0 = 0, x1 = 1;
        while (a > 1) {
            long q = a / m;
            long t = m;
            m = a % m;
            a = t;
            t = x0;
            x0 = x1 - q * x0;
            x1 = t;
        }
        
        if (x1 < 0) x1 += m0;
        return x1;
    }
    
    /**
     * Get the taxa ordering for a specific gene tree.
     */
    public int[] getGeneTreeOrdering(int geneTreeIndex) {
        return geneTreeOrderings[geneTreeIndex];
    }
    
    /**
     * Get all gene tree orderings.
     */
    public int[][] getAllGeneTreeOrderings() {
        return geneTreeOrderings;
    }
    
    public int getNumGeneTrees() {
        return numGeneTrees;
    }
    
    public int getMaxTaxaCount() {
        return geneTreeOrderings.length > 0 ? geneTreeOrderings[0].length : 0;
    }
}
