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
     */
    public Object calculateHash(HashFunction hashFunction, long[][] prefixSums, long[][] prefixXORs) {
        return hashFunction.calculateHash(this, prefixSums, prefixXORs);
    }
    
    /**
     * Backwards compatibility method for existing code.
     */
    public Object calculateHash(HashFunction hashFunction, long[][] prefixSums) {
        return hashFunction.calculateHash(this, prefixSums, null);
    }
    
    /**
     * Hash pair class to preserve both SUM and XOR hash values without information loss.
     */
    public static class HashPair {
        public final long sumHash;
        public final long xorHash;
        
        public HashPair(long sumHash, long xorHash) {
            this.sumHash = sumHash;
            this.xorHash = xorHash;
        }
        
        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof HashPair)) return false;
            HashPair other = (HashPair) obj;
            return this.sumHash == other.sumHash && this.xorHash == other.xorHash;
        }
        
        @Override
        public int hashCode() {
            return (int)(sumHash ^ (sumHash >>> 32)) ^ (int)(xorHash ^ (xorHash >>> 32));
        }
        
        @Override
        public String toString() {
            return "HashPair(" + sumHash + ", " + xorHash + ")";
        }
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof RangeBipartition)) return false;
        RangeBipartition other = (RangeBipartition) obj;
        
        if (this.geneTreeIndex == other.geneTreeIndex) {
            return this.leftStart == other.leftStart && this.leftEnd == other.leftEnd &&
                   this.rightStart == other.rightStart && this.rightEnd == other.rightEnd;
        }
        return false;
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
        return "RangeBipartition[tree=" + geneTreeIndex + 
               ", left=[" + leftStart + "," + leftEnd + "], right=[" + rightStart + "," + rightEnd + "]]";
    }
    
    public int leftSize() {
        return leftEnd - leftStart;
    }
    
    public int rightSize() {
        return rightEnd - rightStart;
    }
    
    // ================================
    // Hash Functions
    // ================================
    public interface HashFunction {
        Object calculateHash(RangeBipartition range, long[][] prefixSums, long[][] prefixXORs);
        String getName();
    }
    
    /**
     * Stronger Sum Hash Function:
     * Combines left and right sums of hashed taxon IDs + range sizes + a mixing step.
     * Normalizes left/right for symmetry (smaller size or sum on left).
     * Note: Works on prefix sums of hashed taxon IDs for better distribution.
     */
    public static class SumHashFunction implements HashFunction {
        private static final long MIX_CONST1 = 0x9e3779b97f4a7c15L; // golden ratio prime
        private static final long MIX_CONST2 = 0xc2b2ae3d27d4eb4fL; // Murmur-like
        
        @Override
        public Long calculateHash(RangeBipartition range, long[][] prefixSums, long[][] prefixXORs) {
            if (prefixSums == null || range.geneTreeIndex >= prefixSums.length) {
                return 0L;
            }
            long[] prefixSum = prefixSums[range.geneTreeIndex];
            if (prefixSum == null) {
                return 0L;
            }
            long leftSum = rangeSum(prefixSum, range.leftStart, range.leftEnd);
            long rightSum = rangeSum(prefixSum, range.rightStart, range.rightEnd);
            int leftSize = range.leftSize();
            int rightSize = range.rightSize();
            
            // Normalize: smaller size on left; if equal, smaller sum on left
            if (leftSize > rightSize || (leftSize == rightSize && leftSum > rightSum)) {
                long tempSum = leftSum;
                leftSum = rightSum;
                rightSum = tempSum;
                int tempSize = leftSize;
                leftSize = rightSize;
                rightSize = tempSize;
            }
            
            long combined = leftSum * MIX_CONST1 ^ rightSum * MIX_CONST2 ^ leftSize ^ (rightSize << 16);
            return Long.rotateLeft(combined, 27) ^ (combined >>> 33);
        }
        
        private long rangeSum(long[] prefix, int start, int end) {
            if (end == 0) return 0;
            long startSum = (start > 0) ? prefix[start - 1] : 0;
            return prefix[end - 1] - startSum;
        }
        
        @Override
        public String getName() {
            return "SUM_STRONG";
        }
    }
    
    /**
     * XOR-based hash function (permutation invariant).
     * Uses XOR of all hashed taxon IDs in each cluster, combined with cluster sizes for better distribution.
     * This provides a genuinely different hash from SUM-based, reducing collisions significantly.
     * Note: Works on prefix XORs of hashed taxon IDs for better distribution.
     */
    public static class XORHashFunction implements HashFunction {
        private static final long XOR_MIX_A = 0xbf58476d1ce4e5b9L;
        private static final long XOR_MIX_B = 0x94d049bb133111ebL;
        
        @Override
        public Long calculateHash(RangeBipartition range, long[][] prefixSums, long[][] prefixXORs) {
            if (prefixXORs == null || range.geneTreeIndex < 0 || range.geneTreeIndex >= prefixXORs.length) {
                // fallback deterministic hash from structural info
                long structural = ((long)range.leftStart << 48) ^ ((long)range.leftEnd << 32) ^ 
                                  ((long)range.rightStart << 16) ^ (long)range.rightEnd;
                return splitMix64Local(structural ^ XOR_MIX_A);
            }
            long[] prefixXOR = prefixXORs[range.geneTreeIndex];
            if (prefixXOR == null) {
                long structural = ((long)range.leftStart << 48) ^ ((long)range.leftEnd << 32) ^ 
                                  ((long)range.rightStart << 16) ^ (long)range.rightEnd;
                return splitMix64Local(structural ^ XOR_MIX_B);
            }
            
            long leftXOR = rangeXOR(prefixXOR, range.leftStart, range.leftEnd);
            long rightXOR = rangeXOR(prefixXOR, range.rightStart, range.rightEnd);
            int leftSize = range.leftSize();
            int rightSize = range.rightSize();
            
            // Normalize: smaller size on "left"; if equal, smaller XOR on "left" (unsigned comparison)
            if (leftSize > rightSize || (leftSize == rightSize && unsignedLessLocal(rightXOR, leftXOR))) {
                long tempXOR = leftXOR;
                leftXOR = rightXOR;
                rightXOR = tempXOR;
                int tempSize = leftSize;
                leftSize = rightSize;
                rightSize = tempSize;
            }
            
            long leftClusterHash = xorClusterMix(leftXOR, leftSize, XOR_MIX_A);
            long rightClusterHash = xorClusterMix(rightXOR, rightSize, XOR_MIX_B);
            
            return combineClusterHashesLocal(leftClusterHash, rightClusterHash, leftSize + rightSize);
        }
        
        private static long rangeXOR(long[] prefixXOR, int start, int end) {
            // Efficient O(1) range XOR using prefix XOR arrays
            if (prefixXOR == null || end <= 0) return 0L;
            
            int e = Math.min(end - 1, prefixXOR.length - 1);
            long endXOR = prefixXOR[e];
            long startXOR = (start > 0) ? prefixXOR[Math.min(start - 1, prefixXOR.length - 1)] : 0L;
            
            // For XOR: rangeXOR(start, end) = prefixXOR[end-1] ^ prefixXOR[start-1]
            return endXOR ^ startXOR;
        }
        
        private static long xorClusterMix(long xorVal, int size, long salt) {
            long v = xorVal ^ (((long) size) << 32) ^ salt;
            return splitMix64Local(v);
        }
        
        private static long combineClusterHashesLocal(long a, long b, int totalSize) {
            long mix = a ^ (b + 0x9e3779b97f4a7c15L + (a << 6) + (a >>> 2));
            mix ^= ((long) totalSize << 17) | ((long) totalSize << 3);
            return splitMix64Local(mix);
        }
        
        private static boolean unsignedLessLocal(long x, long y) {
            return (x + Long.MIN_VALUE) < (y + Long.MIN_VALUE);
        }
        
        private static long splitMix64Local(long z) {
            z += 0x9e3779b97f4a7c15L;
            z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
            z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
            return z ^ (z >>> 31);
        }
        
        @Override
        public String getName() {
            return "XOR_CANONICAL";
        }
    }
    
    // Polynomial hash function removed - not permutation invariant
    
    /**
     * Double Hash Function using SumHashFunction and XORHashFunction.
     * Provides genuine collision resistance by combining two different permutation-invariant hash functions.
     */
    public static class DoubleHashFunction implements HashFunction {
        private final SumHashFunction sumHash;
        private final XORHashFunction xorHash;
        
        public DoubleHashFunction() {
            this.sumHash = new SumHashFunction();
            this.xorHash = new XORHashFunction();
        }
        
        @Override
        public HashPair calculateHash(RangeBipartition range, long[][] prefixSums, long[][] prefixXORs) {
            Long h1 = sumHash.calculateHash(range, prefixSums, prefixXORs);
            Long h2 = xorHash.calculateHash(range, prefixSums, prefixXORs);
            
            // Return both hash values without information loss
            return new HashPair(h1, h2);
        }
        
        @Override
        public String getName() {
            return "DOUBLE_SUM_XOR";
        }
    }
}