package tree;

/**
 * Double hash representation for clusters using sum and XOR hash components.
 * This provides collision-resistant hashing for memory-efficient cluster operations.
 * 
 * Similar to RangeBipartition.HashPair but specifically designed for cluster hashing.
 * A cluster represents a contiguous range of taxa in a gene tree ordering.
 */
public class ClusterHashPair {
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
    
    /**
     * Get a string representation suitable for debugging and logging.
     */
    public String toDebugString() {
        return String.format("ClusterHash[sum=0x%016x, xor=0x%016x]", sumHash, xorHash);
    }
}
