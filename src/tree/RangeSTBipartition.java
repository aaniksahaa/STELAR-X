package tree;

import java.util.Arrays;
import java.util.Objects;

/**
 * Memory-efficient representation of subtree bipartitions using range-based approach.
 * Instead of storing BitSets which consume O(n) space per bipartition, we represent 
 * each bipartition as (geneTreeIndex, startPos, endPos) which consumes O(1) space.
 * 
 * This reduces total memory from O(n²k) to O(nk) where n=taxa, k=gene trees.
 */
public class RangeSTBipartition {
    public final int geneTreeIndex;  // Index of the gene tree this bipartition comes from
    public final int startPos;       // Start position in the left-to-right taxa array
    public final int endPos;         // End position (exclusive) in the left-to-right taxa array
    
    // Cached hash values for efficient comparison
    private long hash1 = 0;  // Primary hash
    private long hash2 = 0;  // Secondary hash for collision resolution
    private int _hashCode = 0;
    
    public RangeSTBipartition(int geneTreeIndex, int startPos, int endPos) {
        this.geneTreeIndex = geneTreeIndex;
        this.startPos = startPos;
        this.endPos = endPos;
    }
    
    /**
     * Computes the primary and secondary hash values using precomputed prefix sums.
     */
    public void computeHashes(RangeHashCalculator hashCalculator) {
        this.hash1 = hashCalculator.getRangeHash1(geneTreeIndex, startPos, endPos);
        this.hash2 = hashCalculator.getRangeHash2(geneTreeIndex, startPos, endPos);
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof RangeSTBipartition)) return false;
        
        RangeSTBipartition other = (RangeSTBipartition) obj;
        
        // If both are from the same gene tree, check exact range match
        if (this.geneTreeIndex == other.geneTreeIndex) {
            return this.startPos == other.startPos && this.endPos == other.endPos;
        }
        
        // If from different gene trees, compare using hash values
        return this.hash1 == other.hash1 && this.hash2 == other.hash2;
    }
    
    @Override
    public int hashCode() {
        if (_hashCode == 0) {
            _hashCode = Objects.hash(hash1, hash2);
        }
        return _hashCode;
    }
    
    @Override
    public String toString() {
        return String.format("RangeSTBip[tree=%d, range=%d-%d, hash1=%d, hash2=%d]", 
            geneTreeIndex, startPos, endPos, hash1, hash2);
    }
    
    /**
     * Gets the taxa IDs in this bipartition range from the specified gene tree ordering.
     */
    public int[] getTaxaIds(int[][] geneTreeOrderings) {
        int[] ordering = geneTreeOrderings[geneTreeIndex];
        return Arrays.copyOfRange(ordering, startPos, endPos);
    }
    
    /**
     * Gets the complement taxa IDs (not in this bipartition) from the specified gene tree ordering.
     */
    public int[] getComplementTaxaIds(int[][] geneTreeOrderings) {
        int[] ordering = geneTreeOrderings[geneTreeIndex];
        int[] result = new int[ordering.length - (endPos - startPos)];
        
        int idx = 0;
        // Add taxa before startPos
        for (int i = 0; i < startPos; i++) {
            result[idx++] = ordering[i];
        }
        // Add taxa after endPos
        for (int i = endPos; i < ordering.length; i++) {
            result[idx++] = ordering[i];
        }
        
        return result;
    }
    
    public long getHash1() { return hash1; }
    public long getHash2() { return hash2; }
    public int size() { return endPos - startPos; }
}
