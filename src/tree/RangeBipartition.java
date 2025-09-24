package tree;

/**
 * Range-based bipartition representation for memory optimization.
 * Instead of storing BitSets (O(n) memory each), we represent bipartitions
 * using the left-to-right ordering of gene tree taxa with two ranges:
 * one for the left subtree and one for the right subtree.
 * 
 * Memory usage: O(1) per bipartition instead of O(n)
 */
public class RangeBipartition {
    public final int geneTreeIndex;    // Which gene tree this bipartition belongs to
    public final int leftStart;       // Start index of left subtree range (inclusive)
    public final int leftEnd;         // End index of left subtree range (exclusive)
    public final int rightStart;      // Start index of right subtree range (inclusive) 
    public final int rightEnd;        // End index of right subtree range (exclusive)
    
    // Cached hash for efficiency
    private int _hash = 0;
    
    public RangeBipartition(int geneTreeIndex, int leftStart, int leftEnd, int rightStart, int rightEnd) {
        this.geneTreeIndex = geneTreeIndex;
        this.leftStart = leftStart;
        this.leftEnd = leftEnd;
        this.rightStart = rightStart;
        this.rightEnd = rightEnd;
    }
    
    /**
     * Calculates hash using configurable hash functions.
     * Default: simple sum of taxon IDs in the range.
     */
    public long calculateHash(HashFunction hashFunction, long[][] prefixSums) {
        return hashFunction.calculateHash(this, prefixSums);
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof RangeBipartition)) return false;
        RangeBipartition other = (RangeBipartition) obj;
        
        // If same gene tree, compare ranges directly
        if (this.geneTreeIndex == other.geneTreeIndex) {
            return this.leftStart == other.leftStart && this.leftEnd == other.leftEnd &&
                   this.rightStart == other.rightStart && this.rightEnd == other.rightEnd;
        }
        
        // Different gene trees - equality will be determined by hash comparison
        // This is handled at a higher level where prefix sums are available
        return false;
    }
    
    @Override
    public int hashCode() {
        if (_hash == 0) {
            // Simple hash based on gene tree index and ranges
            // The actual hash used for equality will be calculated separately
            _hash = 31 * (31 * (31 * (31 * geneTreeIndex + leftStart) + leftEnd) + rightStart) + rightEnd;
        }
        return _hash;
    }
    
    @Override
    public String toString() {
        return "RangeBipartition[tree=" + geneTreeIndex + ", left=[" + leftStart + "," + leftEnd + "], right=[" + rightStart + "," + rightEnd + "]]";
    }
    
    /**
     * Returns the size of the left cluster.
     */
    public int leftSize() {
        return leftEnd - leftStart;
    }
    
    /**
     * Returns the size of the right cluster.
     */
    public int rightSize() {
        return rightEnd - rightStart;
    }
    
    /**
     * Interface for configurable hash functions.
     */
    public interface HashFunction {
        long calculateHash(RangeBipartition range, long[][] prefixSums);
        String getName();
    }
    
    /**
     * Simple sum hash function - sums all taxon IDs in both left and right ranges.
     */
    public static class SumHashFunction implements HashFunction {
        @Override
        public long calculateHash(RangeBipartition range, long[][] prefixSums) {
            if (prefixSums == null || range.geneTreeIndex >= prefixSums.length) {
                return 0;
            }
            
            long[] prefixSum = prefixSums[range.geneTreeIndex];
            if (prefixSum == null) {
                return 0;
            }
            
            // Calculate sum for left range
            long leftSum = 0;
            if (range.leftEnd > 0 && range.leftEnd <= prefixSum.length) {
                long leftStartSum = (range.leftStart > 0) ? prefixSum[range.leftStart - 1] : 0;
                long leftEndSum = prefixSum[range.leftEnd - 1];
                leftSum = leftEndSum - leftStartSum;
            }
            
            // Calculate sum for right range  
            long rightSum = 0;
            if (range.rightEnd > 0 && range.rightEnd <= prefixSum.length) {
                long rightStartSum = (range.rightStart > 0) ? prefixSum[range.rightStart - 1] : 0;
                long rightEndSum = prefixSum[range.rightEnd - 1];
                rightSum = rightEndSum - rightStartSum;
            }
            
            // Combine the sums - we use the left sum as the primary hash
            // since that represents cluster1 in the STBipartition
            // We could also combine both: return (leftSum << 32) | (rightSum & 0xFFFFFFFFL);
            return leftSum;
        }
        
        @Override
        public String getName() {
            return "SUM";
        }
    }
    
    /**
     * Polynomial hash function with configurable base and modulus.
     */
    public static class PolynomialHashFunction implements HashFunction {
        private final long base;
        private final long mod;
        
        public PolynomialHashFunction(long base, long mod) {
            this.base = base;
            this.mod = mod;
        }
        
        public PolynomialHashFunction() {
            this(31, 1000000007L); // Default values
        }
        
        @Override
        public long calculateHash(RangeBipartition range, long[][] prefixSums) {
            // This would require polynomial prefix sums - for now, fall back to sum
            SumHashFunction sumHash = new SumHashFunction();
            long sum = sumHash.calculateHash(range, prefixSums);
            return ((sum * base) % mod + range.leftSize() + range.rightSize()) % mod;
        }
        
        @Override
        public String getName() {
            return "POLYNOMIAL_" + base + "_" + mod;
        }
    }
    
    /**
     * Double hash function combining two different hash functions.
     */
    public static class DoubleHashFunction implements HashFunction {
        private final HashFunction hash1;
        private final HashFunction hash2;
        
        public DoubleHashFunction(HashFunction hash1, HashFunction hash2) {
            this.hash1 = hash1;
            this.hash2 = hash2;
        }
        
        @Override
        public long calculateHash(RangeBipartition range, long[][] prefixSums) {
            long h1 = hash1.calculateHash(range, prefixSums);
            long h2 = hash2.calculateHash(range, prefixSums);
            return (h1 << 32) | (h2 & 0xFFFFFFFFL);
        }
        
        @Override
        public String getName() {
            return "DOUBLE_" + hash1.getName() + "_" + hash2.getName();
        }
    }
}
