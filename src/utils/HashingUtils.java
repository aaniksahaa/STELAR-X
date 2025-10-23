package utils;

import java.util.Set;

/**
 * Unified hashing utilities for consistent hash calculation across 
 * RangeBipartition and ClusterHash implementations.
 * 
 * This ensures that the same taxa produce the same hash regardless of 
 * whether they're represented as ranges or as explicit taxon sets.
 */
public class HashingUtils {
    
    // Constants for consistent hashing across all implementations
    private static final long MIX_CONST1 = 0x9e3779b97f4a7c15L; // golden ratio prime
    private static final long MIX_CONST2 = 0xc2b2ae3d27d4eb4fL; // Murmur-like
    private static final long XOR_MIX_A = 0x517cc1b727220a95L;
    private static final long XOR_MIX_B = 0x85ebca6b8e4a7c15L;
    
    /**
     * Hash a single taxon ID consistently across all implementations.
     */
    public static long hashSingleTaxon(int taxonId) {
        long x = taxonId;
        x ^= x >>> 16;
        x *= 0x85ebca6b;
        x ^= x >>> 13;
        x *= 0xc2b2ae35;
        x ^= x >>> 16;
        return x == 0 ? 1 : x;
    }
    
    /**
     * Calculate range sum from prefix sum array.
     */
    public static long rangeSum(long[] prefixSum, int start, int end) {
        if (prefixSum == null || end <= 0) return 0L;
        
        int e = Math.min(end - 1, prefixSum.length - 1);
        long endSum = prefixSum[e];
        long startSum = (start > 0) ? prefixSum[Math.min(start - 1, prefixSum.length - 1)] : 0L;
        
        return endSum - startSum;
    }
    
    /**
     * Calculate range XOR from prefix XOR array.
     */
    public static long rangeXOR(long[] prefixXOR, int start, int end) {
        if (prefixXOR == null || end <= 0) return 0L;
        
        int e = Math.min(end - 1, prefixXOR.length - 1);
        long endXOR = prefixXOR[e];
        long startXOR = (start > 0) ? prefixXOR[Math.min(start - 1, prefixXOR.length - 1)] : 0L;
        
        // For XOR: rangeXOR(start, end) = prefixXOR[end-1] ^ prefixXOR[start-1]
        return endXOR ^ startXOR;
    }
    
    /**
     * Calculate sum hash for a cluster using the SUM hash function approach.
     * This works for both ranges (via prefix arrays) and explicit taxon sets.
     */
    public static long calculateSumHash(long leftSum, long rightSum, int leftSize, int rightSize) {
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
    
    /**
     * Calculate XOR hash for a cluster using the XOR hash function approach.
     * This works for both ranges (via prefix arrays) and explicit taxon sets.
     */
    public static long calculateXORHash(long leftXOR, long rightXOR, int leftSize, int rightSize) {
        // Normalize: smaller size on left; if equal, smaller XOR on left
        if (leftSize > rightSize || (leftSize == rightSize && unsignedLess(rightXOR, leftXOR))) {
            long tempXOR = leftXOR;
            leftXOR = rightXOR;
            rightXOR = tempXOR;
            int tempSize = leftSize;
            leftSize = rightSize;
            rightSize = tempSize;
        }
        
        long leftClusterHash = xorClusterMix(leftXOR, leftSize, XOR_MIX_A);
        long rightClusterHash = xorClusterMix(rightXOR, rightSize, XOR_MIX_B);
        
        return combineClusterHashes(leftClusterHash, rightClusterHash, leftSize + rightSize);
    }
    
    /**
     * Calculate hash for a single cluster (no bipartition).
     * This is used when we have just one cluster and need its hash.
     */
    public static long calculateSingleClusterSumHash(long sum, int size) {
        // For a single cluster, we treat it as if it's the "left" cluster with empty "right"
        return calculateSumHash(sum, 0L, size, 0);
    }
    
    /**
     * Calculate hash for a single cluster using XOR approach.
     */
    public static long calculateSingleClusterXORHash(long xor, int size) {
        // For a single cluster, we treat it as if it's the "left" cluster with empty "right"
        return calculateXORHash(xor, 0L, size, 0);
    }
    
    /**
     * Calculate sum and XOR for a set of taxon IDs.
     * Returns [sum, xor] array.
     */
    public static long[] calculateTaxonSetSumAndXOR(Set<Integer> taxonIds) {
        long sum = 0;
        long xor = 0;
        
        for (int taxonId : taxonIds) {
            long hashedTaxon = hashSingleTaxon(taxonId);
            sum += hashedTaxon;
            xor ^= hashedTaxon;
        }
        
        return new long[]{sum, xor};
    }
    
    /**
     * Calculate sum and XOR for a range using prefix arrays.
     * Returns [sum, xor] array.
     */
    public static long[] calculateRangeSumAndXOR(long[] prefixSum, long[] prefixXOR, int start, int end) {
        long sum = rangeSum(prefixSum, start, end);
        long xor = rangeXOR(prefixXOR, start, end);
        return new long[]{sum, xor};
    }
    
    // Helper methods (private)
    
    private static long xorClusterMix(long xorVal, int size, long salt) {
        long v = xorVal ^ (((long) size) << 32) ^ salt;
        return splitMix64(v);
    }
    
    private static long combineClusterHashes(long a, long b, int totalSize) {
        long mix = a ^ (b + 0x9e3779b97f4a7c15L + (a << 6) + (a >>> 2));
        mix ^= ((long) totalSize << 17) | ((long) totalSize << 3);
        return splitMix64(mix);
    }
    
    private static boolean unsignedLess(long x, long y) {
        return (x + Long.MIN_VALUE) < (y + Long.MIN_VALUE);
    }
    
    private static long splitMix64(long z) {
        z += 0x9e3779b97f4a7c15L;
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }
}
