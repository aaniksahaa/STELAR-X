package utils;

import tree.ClusterHashPair;

/**
 * Hash utilities for memory-efficient cluster operations.
 * Provides cluster-specific hash computation using the same double-hashing scheme
 * as RangeBipartition but isolated for cluster use cases.
 * 
 * Uses sum and XOR hash functions for collision resistance, operating on
 * prefix arrays of hashed taxon IDs for better distribution.
 */
public class HashUtils {
    
    // Hash mixing constants (same as RangeBipartition for consistency)
    private static final long MIX_CONST1 = 0x9e3779b97f4a7c15L; // golden ratio prime
    private static final long XOR_MIX_A = 0xbf58476d1ce4e5b9L;
    private static final long XOR_MIX_B = 0x94d049bb133111ebL;
    
    /**
     * Compute double hash (sum + XOR) for a cluster defined by range in a gene tree.
     * 
     * @param geneTreeIndex Index of the gene tree
     * @param start Start position in gene tree ordering (inclusive)
     * @param end End position in gene tree ordering (exclusive)
     * @param orderings Gene tree orderings [tree][position] = taxonId
     * @param prefixSums Prefix sums of hashed taxon IDs [tree][position]
     * @param prefixXORs Prefix XORs of hashed taxon IDs [tree][position]
     * @return ClusterHashPair containing both sum and XOR hash components
     */
    public static ClusterHashPair computeClusterHash(
            int geneTreeIndex, int start, int end,
            int[][] orderings, long[][] prefixSums, long[][] prefixXORs) {
        
        // Compute sum hash (similar to RangeBipartition.SumHashFunction)
        long sumHash = computeSumHash(geneTreeIndex, start, end, prefixSums);
        
        // Compute XOR hash (similar to RangeBipartition.XORHashFunction)  
        long xorHash = computeXORHash(geneTreeIndex, start, end, prefixXORs);
        
        return new ClusterHashPair(sumHash, xorHash);
    }
    
    /**
     * Compute sum-based hash for a cluster range.
     * Uses prefix sums for O(1) range sum computation.
     */
    private static long computeSumHash(int treeIndex, int start, int end, long[][] prefixSums) {
        if (prefixSums == null || treeIndex >= prefixSums.length || treeIndex < 0) {
            return computeFallbackSumHash(treeIndex, start, end);
        }
        
        long[] prefixSum = prefixSums[treeIndex];
        if (prefixSum == null || prefixSum.length == 0) {
            return computeFallbackSumHash(treeIndex, start, end);
        }
        
        long rangeSum = rangeSum(prefixSum, start, end);
        int size = end - start;
        
        // Apply mixing similar to RangeBipartition.SumHashFunction
        // Normalize by including size information for better distribution
        long combined = rangeSum * MIX_CONST1 ^ ((long) size << 16);
        return Long.rotateLeft(combined, 27) ^ (combined >>> 33);
    }
    
    /**
     * Compute XOR-based hash for a cluster range.
     * Uses prefix XORs for O(1) range XOR computation.
     */
    private static long computeXORHash(int treeIndex, int start, int end, long[][] prefixXORs) {
        if (prefixXORs == null || treeIndex >= prefixXORs.length || treeIndex < 0) {
            return computeFallbackXORHash(treeIndex, start, end);
        }
        
        long[] prefixXOR = prefixXORs[treeIndex];
        if (prefixXOR == null || prefixXOR.length == 0) {
            return computeFallbackXORHash(treeIndex, start, end);
        }
        
        long rangeXOR = rangeXOR(prefixXOR, start, end);
        int size = end - start;
        
        // Apply mixing similar to RangeBipartition.XORHashFunction
        return xorClusterMix(rangeXOR, size, XOR_MIX_A);
    }
    
    /**
     * Compute range sum using prefix sum array: sum(start, end) = prefix[end-1] - prefix[start-1]
     */
    private static long rangeSum(long[] prefix, int start, int end) {
        if (end <= 0 || prefix.length == 0) return 0;
        
        int safeEnd = Math.min(end - 1, prefix.length - 1);
        long endSum = (safeEnd >= 0) ? prefix[safeEnd] : 0;
        long startSum = (start > 0 && start - 1 < prefix.length) ? prefix[start - 1] : 0;
        
        return endSum - startSum;
    }
    
    /**
     * Compute range XOR using prefix XOR array: xor(start, end) = prefix[end-1] ^ prefix[start-1]
     */
    private static long rangeXOR(long[] prefixXOR, int start, int end) {
        if (prefixXOR == null || end <= 0 || prefixXOR.length == 0) return 0L;
        
        int safeEnd = Math.min(end - 1, prefixXOR.length - 1);
        long endXOR = (safeEnd >= 0) ? prefixXOR[safeEnd] : 0L;
        long startXOR = (start > 0 && start - 1 < prefixXOR.length) ? prefixXOR[start - 1] : 0L;
        
        // For XOR: rangeXOR(start, end) = prefixXOR[end-1] ^ prefixXOR[start-1]
        return endXOR ^ startXOR;
    }
    
    /**
     * Apply XOR-based mixing for cluster hash.
     */
    private static long xorClusterMix(long xorVal, int size, long salt) {
        long v = xorVal ^ (((long) size) << 32) ^ salt;
        return splitMix64(v);
    }
    
    /**
     * SplitMix64 hash function for final mixing.
     */
    private static long splitMix64(long z) {
        z += 0x9e3779b97f4a7c15L;
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }
    
    /**
     * Fallback sum hash when prefix arrays are not available.
     * Uses structural information for deterministic hashing.
     */
    private static long computeFallbackSumHash(int treeIndex, int start, int end) {
        long structural = ((long) treeIndex << 32) ^ ((long) start << 16) ^ (long) end;
        return splitMix64(structural ^ MIX_CONST1);
    }
    
    /**
     * Fallback XOR hash when prefix arrays are not available.
     * Uses structural information for deterministic hashing.
     */
    private static long computeFallbackXORHash(int treeIndex, int start, int end) {
        long structural = ((long) treeIndex << 32) ^ ((long) start << 16) ^ (long) end;
        return splitMix64(structural ^ XOR_MIX_B);
    }
    
    /**
     * Hash a single taxon ID using the same method as MemoryEfficientBipartitionManager.
     * This ensures consistency with existing hash computations.
     */
    public static long hashSingleTaxon(int taxonId) {
        // Use the same hash mixing function as MemoryEfficientBipartitionManager
        long x = taxonId;
        x ^= x >>> 16;
        x *= 0x85ebca6b;
        x ^= x >>> 13;
        x *= 0xc2b2ae35;
        x ^= x >>> 16;
        
        // Ensure we never return 0 to avoid issues with prefix operations
        return x == 0 ? 1 : x;
    }
    
    /**
     * Compute a deterministic hash for a set of taxon IDs.
     * Used as fallback when range-based hashing is not available.
     */
    public static ClusterHashPair computeFallbackClusterHash(java.util.Set<Integer> taxonSet) {
        long sumHash = 0;
        long xorHash = 0;
        
        // Process taxon IDs in sorted order for deterministic results
        java.util.List<Integer> sortedTaxa = new java.util.ArrayList<>(taxonSet);
        java.util.Collections.sort(sortedTaxa);
        
        for (int taxonId : sortedTaxa) {
            long hashedTaxon = hashSingleTaxon(taxonId);
            sumHash += hashedTaxon;
            xorHash ^= hashedTaxon;
        }
        
        // Apply final mixing
        sumHash = Long.rotateLeft(sumHash * MIX_CONST1, 27);
        xorHash = splitMix64(xorHash ^ (taxonSet.size() << 16));
        
        return new ClusterHashPair(sumHash, xorHash);
    }
}
