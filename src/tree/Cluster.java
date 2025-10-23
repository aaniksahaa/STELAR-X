package tree;

import utils.HashUtils;

/**
 * Represents a cluster as a contiguous range in a gene tree ordering.
 * This provides a memory-efficient alternative to BitSet representations.
 * 
 * A cluster is defined by:
 * - geneTreeIndex: which gene tree this cluster belongs to
 * - startPos: start position in the gene tree's left-to-right ordering (inclusive)
 * - endPos: end position in the gene tree's left-to-right ordering (exclusive)
 * 
 * Memory usage: O(1) per cluster instead of O(n) for BitSet
 */
public class Cluster {
    public final int geneTreeIndex;
    public final int startPos;
    public final int endPos;
    private ClusterHashPair cachedHash;
    
    public Cluster(int geneTreeIndex, int startPos, int endPos) {
        if (startPos < 0 || endPos < startPos) {
            throw new IllegalArgumentException("Invalid cluster range: [" + startPos + ", " + endPos + ")");
        }
        
        this.geneTreeIndex = geneTreeIndex;
        this.startPos = startPos;
        this.endPos = endPos;
    }
    
    /**
     * Compute double hash for this cluster using the provided hash utilities.
     * Results are cached for efficiency.
     */
    public ClusterHashPair computeHash(int[][] geneTreeOrderings,
                                      long[][] prefixSums, 
                                      long[][] prefixXORs) {
        if (cachedHash == null) {
            cachedHash = HashUtils.computeClusterHash(
                geneTreeIndex, startPos, endPos,
                geneTreeOrderings, prefixSums, prefixXORs);
        }
        return cachedHash;
    }
    
    /**
     * Get the size (number of taxa) in this cluster.
     */
    public int size() {
        return endPos - startPos;
    }
    
    /**
     * Check if this cluster is empty.
     */
    public boolean isEmpty() {
        return startPos >= endPos;
    }
    
    /**
     * Check if this cluster contains a single taxon.
     */
    public boolean isSingleton() {
        return size() == 1;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof Cluster)) return false;
        Cluster other = (Cluster) obj;
        
        return this.geneTreeIndex == other.geneTreeIndex &&
               this.startPos == other.startPos &&
               this.endPos == other.endPos;
    }
    
    @Override
    public int hashCode() {
        return 31 * (31 * geneTreeIndex + startPos) + endPos;
    }
    
    @Override
    public String toString() {
        return "Cluster[tree=" + geneTreeIndex + 
               ", range=[" + startPos + "," + endPos + "), size=" + size() + "]";
    }
    
    /**
     * Get a detailed string representation for debugging.
     */
    public String toDebugString() {
        return String.format("Cluster[tree=%d, range=[%d,%d), size=%d, hash=%s]", 
                           geneTreeIndex, startPos, endPos, size(),
                           cachedHash != null ? cachedHash.toDebugString() : "not computed");
    }
}
