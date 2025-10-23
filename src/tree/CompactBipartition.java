package tree;

/**
 * Compact bipartition representation for memory-efficient weight calculation.
 * Instead of storing BitSets, we represent bipartitions using range information
 * and provide methods for efficient intersection counting using inverse indices.
 */
public class CompactBipartition {
    public final int geneTreeIndex;    // Which gene tree this bipartition belongs to
    public final int leftStart;       // Start index of left cluster range (inclusive)
    public final int leftEnd;         // End index of left cluster range (exclusive)
    public final int rightStart;      // Start index of right cluster range (inclusive) 
    public final int rightEnd;        // End index of right cluster range (exclusive)
    
    // Cached hash for efficiency
    private int _hash = 0;
    private Object _permutationInvariantHash = null;
    
    public CompactBipartition(int geneTreeIndex, int leftStart, int leftEnd, int rightStart, int rightEnd) {
        this.geneTreeIndex = geneTreeIndex;
        this.leftStart = leftStart;
        this.leftEnd = leftEnd;
        this.rightStart = rightStart;
        this.rightEnd = rightEnd;
    }
    
    /**
     * Create CompactBipartition from RangeBipartition
     */
    public CompactBipartition(RangeBipartition rangeBip) {
        this(rangeBip.geneTreeIndex, rangeBip.leftStart, rangeBip.leftEnd, 
             rangeBip.rightStart, rangeBip.rightEnd);
    }
    
    /**
     * Calculate permutation-invariant hash for this bipartition using the same
     * hash function as RangeBipartition for consistency.
     */
    public Object calculatePermutationInvariantHash(RangeBipartition.HashFunction hashFunction, 
                                                   long[][] prefixSums, long[][] prefixXORs) {
        if (_permutationInvariantHash == null) {
            RangeBipartition rangeBip = new RangeBipartition(geneTreeIndex, leftStart, leftEnd, rightStart, rightEnd);
            _permutationInvariantHash = hashFunction.calculateHash(rangeBip, prefixSums, prefixXORs);
        }
        return _permutationInvariantHash;
    }
    
    /**
     * Calculate intersection size between two ranges using inverse index mapping.
     * This is the core optimization described in the mathematical document.
     * 
     * @param rangeA_start Start of range A (inclusive)
     * @param rangeA_end End of range A (exclusive)
     * @param rangeB_start Start of range B (inclusive)
     * @param rangeB_end End of range B (exclusive)
     * @param geneTreeA Gene tree index for range A
     * @param geneTreeB Gene tree index for range B
     * @param inverseIndex Inverse index matrix [geneTree][taxonId] = position
     * @param geneTreeOrderings Gene tree orderings [geneTree][position] = taxonId
     * @return Intersection size |A ∩ B|
     */
    public static int calculateRangeIntersection(int rangeA_start, int rangeA_end,
                                               int rangeB_start, int rangeB_end,
                                               int geneTreeA, int geneTreeB,
                                               int[][] inverseIndex,
                                               int[][] geneTreeOrderings) {
        
        // Determine which range is smaller to iterate over
        int sizeA = rangeA_end - rangeA_start;
        int sizeB = rangeB_end - rangeB_start;
        
        if (sizeA <= sizeB) {
            return countIntersection(rangeA_start, rangeA_end, geneTreeA,
                                   rangeB_start, rangeB_end, geneTreeB,
                                   inverseIndex, geneTreeOrderings);
        } else {
            return countIntersection(rangeB_start, rangeB_end, geneTreeB,
                                   rangeA_start, rangeA_end, geneTreeA,
                                   inverseIndex, geneTreeOrderings);
        }
    }
    
    /**
     * Count intersection by iterating over the smaller range.
     * For each element v in the smaller range, check if its position in the other
     * gene tree falls within the target range using O(1) inverse index lookup.
     */
    private static int countIntersection(int smallerStart, int smallerEnd, int smallerTreeIdx,
                                       int targetStart, int targetEnd, int targetTreeIdx,
                                       int[][] inverseIndex, int[][] geneTreeOrderings) {
        int count = 0;
        
        // Iterate over elements in the smaller range
        for (int pos = smallerStart; pos < smallerEnd; pos++) {
            if (pos >= geneTreeOrderings[smallerTreeIdx].length) break;
            
            int taxonId = geneTreeOrderings[smallerTreeIdx][pos];
            
            // Check if this taxon's position in the target tree falls within target range
            if (targetTreeIdx < inverseIndex.length && taxonId < inverseIndex[targetTreeIdx].length) {
                int positionInTarget = inverseIndex[targetTreeIdx][taxonId];
                
                if (positionInTarget >= targetStart && positionInTarget < targetEnd) {
                    count++;
                }
            }
        }
        
        return count;
    }
    
    /**
     * Calculate all four intersection counts needed for weight calculation between two bipartitions.
     * Returns array: [|A1 ∩ A2|, |A1 ∩ B2|, |B1 ∩ A2|, |B1 ∩ B2|]
     */
    public int[] calculateAllIntersections(CompactBipartition other,
                                         int[][] inverseIndex,
                                         int[][] geneTreeOrderings) {
        int[] intersections = new int[4];
        
        // |A1 ∩ A2| - left cluster of this with left cluster of other
        intersections[0] = calculateRangeIntersection(
            this.leftStart, this.leftEnd, other.leftStart, other.leftEnd,
            this.geneTreeIndex, other.geneTreeIndex, inverseIndex, geneTreeOrderings);
        
        // |A1 ∩ B2| - left cluster of this with right cluster of other
        intersections[1] = calculateRangeIntersection(
            this.leftStart, this.leftEnd, other.rightStart, other.rightEnd,
            this.geneTreeIndex, other.geneTreeIndex, inverseIndex, geneTreeOrderings);
        
        // |B1 ∩ A2| - right cluster of this with left cluster of other
        intersections[2] = calculateRangeIntersection(
            this.rightStart, this.rightEnd, other.leftStart, other.leftEnd,
            this.geneTreeIndex, other.geneTreeIndex, inverseIndex, geneTreeOrderings);
        
        // |B1 ∩ B2| - right cluster of this with right cluster of other
        intersections[3] = calculateRangeIntersection(
            this.rightStart, this.rightEnd, other.rightStart, other.rightEnd,
            this.geneTreeIndex, other.geneTreeIndex, inverseIndex, geneTreeOrderings);
        
        return intersections;
    }
    
    /**
     * Calculate the score between two bipartitions using the compact representation.
     * This implements the same scoring logic as the original BitSet-based method
     * but uses range-based intersection counting.
     */
    public double calculateScore(CompactBipartition other,
                               int[][] inverseIndex,
                               int[][] geneTreeOrderings) {
        int[] intersections = calculateAllIntersections(other, inverseIndex, geneTreeOrderings);
        
        // First configuration: (A1|B1) with (A2|B2)
        int p1 = intersections[0]; // |A1 ∩ A2|
        int p2 = intersections[3]; // |B1 ∩ B2|
        
        double score1 = 0.0;
        if (p1 + p2 >= 2) {
            score1 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        // Second configuration: (A1|B1) with (B2|A2) - cross configuration
        p1 = intersections[1]; // |A1 ∩ B2|
        p2 = intersections[2]; // |B1 ∩ A2|
        
        double score2 = 0.0;
        if (p1 + p2 >= 2) {
            score2 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        return score1 + score2;
    }
    
    public int leftSize() {
        return leftEnd - leftStart;
    }
    
    public int rightSize() {
        return rightEnd - rightStart;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof CompactBipartition)) return false;
        CompactBipartition other = (CompactBipartition) obj;
        
        return this.geneTreeIndex == other.geneTreeIndex &&
               this.leftStart == other.leftStart && this.leftEnd == other.leftEnd &&
               this.rightStart == other.rightStart && this.rightEnd == other.rightEnd;
    }
    
    @Override
    public int hashCode() {
        if (_hash == 0) {
            _hash = 31 * (31 * (31 * (31 * geneTreeIndex + leftStart) + leftEnd) + rightStart) + rightEnd;
        }
        return _hash;
    }
    
    @Override
    public String toString() {
        return "CompactBipartition[tree=" + geneTreeIndex + 
               ", left=[" + leftStart + "," + leftEnd + "), right=[" + rightStart + "," + rightEnd + ")]";
    }
}
