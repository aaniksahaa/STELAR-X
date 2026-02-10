package tree;

/**
 * UnrootedCluster: Represents a cluster in an unrooted tree context.
 * 
 * In ASTRAL's unrooted tree treatment, each edge creates a bipartition with two sides.
 * A cluster can be represented as either:
 * 1. A contiguous range in the tree's leaf ordering (isComplement = false)
 * 2. The complement of a contiguous range (isComplement = true)
 * 
 * For example, in tree ((A,(B,C)),(D,E)) with ordering [A,B,C,D,E]:
 * - Cluster {B,C} = range [1,3), isComplement = false
 * - Cluster {A,D,E} = range [1,3), isComplement = true (complement of {B,C})
 * 
 * This representation allows O(1) storage while supporting all clusters needed
 * for ASTRAL-style quartet scoring on unrooted trees.
 */
public class UnrootedCluster {
    
    public final int treeIndex;        // Which gene tree this cluster belongs to
    public final int rangeStart;       // Start index of the range (inclusive)
    public final int rangeEnd;         // End index of the range (exclusive)
    public final boolean isComplement; // True if this represents the complement of the range
    
    // Cached hash values
    private long _sumHash = 0;
    private long _xorHash = 0;
    private boolean _hashComputed = false;
    
    // Cached size
    private int _size = -1;
    
    public UnrootedCluster(int treeIndex, int rangeStart, int rangeEnd, boolean isComplement) {
        this.treeIndex = treeIndex;
        this.rangeStart = rangeStart;
        this.rangeEnd = rangeEnd;
        this.isComplement = isComplement;
    }
    
    /**
     * Create a direct (non-complement) cluster from a range.
     */
    public static UnrootedCluster fromRange(int treeIndex, int rangeStart, int rangeEnd) {
        return new UnrootedCluster(treeIndex, rangeStart, rangeEnd, false);
    }
    
    /**
     * Create a complement cluster from a range.
     */
    public static UnrootedCluster complementOf(int treeIndex, int rangeStart, int rangeEnd) {
        return new UnrootedCluster(treeIndex, rangeStart, rangeEnd, true);
    }
    
    /**
     * Get the complement of this cluster.
     */
    public UnrootedCluster complement() {
        return new UnrootedCluster(treeIndex, rangeStart, rangeEnd, !isComplement);
    }
    
    /**
     * Calculate the size of this cluster.
     * 
     * @param treeTaxaCount Total number of taxa in this tree
     * @return Size of the cluster
     */
    public int size(int treeTaxaCount) {
        if (_size < 0) {
            int rangeSize = rangeEnd - rangeStart;
            _size = isComplement ? (treeTaxaCount - rangeSize) : rangeSize;
        }
        return _size;
    }
    
    /**
     * Calculate hash values using prefix sums and XORs.
     * For complements, we use: complement_hash = total_hash - range_hash (for sums)
     *                          complement_hash = total_hash ^ range_hash (for XORs)
     * 
     * @param prefixSums Prefix sums array for this tree
     * @param prefixXORs Prefix XORs array for this tree
     * @param totalSumHash Total sum hash for all taxa in this tree
     * @param totalXORHash Total XOR hash for all taxa in this tree
     */
    public void computeHash(long[] prefixSums, long[] prefixXORs, 
                           long totalSumHash, long totalXORHash) {
        if (_hashComputed) return;
        
        // Calculate range hash
        long rangeSumHash = rangeSum(prefixSums, rangeStart, rangeEnd);
        long rangeXORHash = rangeXOR(prefixXORs, rangeStart, rangeEnd);
        
        if (isComplement) {
            // Complement hash = total - range (for sums)
            // Complement hash = total ^ range (for XORs, since A ^ B ^ B = A)
            _sumHash = totalSumHash - rangeSumHash;
            _xorHash = totalXORHash ^ rangeXORHash;
        } else {
            _sumHash = rangeSumHash;
            _xorHash = rangeXORHash;
        }
        
        _hashComputed = true;
    }
    
    /**
     * Get the sum hash (must call computeHash first).
     */
    public long getSumHash() {
        if (!_hashComputed) {
            throw new IllegalStateException("Hash not computed. Call computeHash() first.");
        }
        return _sumHash;
    }
    
    /**
     * Get the XOR hash (must call computeHash first).
     */
    public long getXORHash() {
        if (!_hashComputed) {
            throw new IllegalStateException("Hash not computed. Call computeHash() first.");
        }
        return _xorHash;
    }
    
    /**
     * Get a combined hash pair for use in hash maps.
     */
    public ClusterHashPair getHashPair() {
        if (!_hashComputed) {
            throw new IllegalStateException("Hash not computed. Call computeHash() first.");
        }
        return new ClusterHashPair(_sumHash, _xorHash);
    }
    
    /**
     * Calculate range sum from prefix sums.
     */
    private static long rangeSum(long[] prefix, int start, int end) {
        if (prefix == null || end <= 0 || end > prefix.length) return 0;
        long endSum = prefix[end - 1];
        long startSum = (start > 0) ? prefix[start - 1] : 0;
        return endSum - startSum;
    }
    
    /**
     * Calculate range XOR from prefix XORs.
     */
    private static long rangeXOR(long[] prefix, int start, int end) {
        if (prefix == null || end <= 0 || end > prefix.length) return 0;
        long endXOR = prefix[end - 1];
        long startXOR = (start > 0) ? prefix[start - 1] : 0;
        return endXOR ^ startXOR;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof UnrootedCluster)) return false;
        UnrootedCluster other = (UnrootedCluster) obj;
        return this.treeIndex == other.treeIndex &&
               this.rangeStart == other.rangeStart &&
               this.rangeEnd == other.rangeEnd &&
               this.isComplement == other.isComplement;
    }
    
    @Override
    public int hashCode() {
        int result = 31 * treeIndex + rangeStart;
        result = 31 * result + rangeEnd;
        result = 31 * result + (isComplement ? 1 : 0);
        return result;
    }
    
    @Override
    public String toString() {
        String complement = isComplement ? "COMPLEMENT" : "DIRECT";
        return "UnrootedCluster[tree=" + treeIndex + 
               ", range=[" + rangeStart + "," + rangeEnd + "), " + complement + "]";
    }
    
    /**
     * Hash pair class for cluster identification across trees.
     */
    public static class ClusterHashPair {
        public final long sumHash;
        public final long xorHash;
        
        public ClusterHashPair(long sumHash, long xorHash) {
            this.sumHash = sumHash;
            this.xorHash = xorHash;
        }
        
        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof ClusterHashPair)) return false;
            ClusterHashPair other = (ClusterHashPair) obj;
            return this.sumHash == other.sumHash && this.xorHash == other.xorHash;
        }
        
        @Override
        public int hashCode() {
            return (int)(sumHash ^ (sumHash >>> 32)) ^ (int)(xorHash ^ (xorHash >>> 32));
        }
        
        @Override
        public String toString() {
            return "ClusterHashPair(" + sumHash + ", " + xorHash + ")";
        }
    }
}
