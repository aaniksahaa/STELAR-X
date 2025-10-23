package tree;

import java.util.Set;
import java.util.HashSet;
import utils.HashingUtils;

/**
 * Permutation-invariant hash representation for clusters in DP.
 * Instead of using BitSet as keys in DP maps, we use hash-based cluster representation
 * that is both memory-efficient and provides fast equality checking.
 */
public class ClusterHash {
    private final Object hash;
    private final int size;
    
    // Cache for reconstructed taxon set (lazy initialization)
    private Set<Integer> taxonSet = null;
    
    public ClusterHash(Object hash, int size) {
        this.hash = hash;
        this.size = size;
    }
    
    /**
     * Create ClusterHash from a CompactBipartition's union (left + right clusters).
     * This correctly calculates the hash for the actual union of taxa.
     */
    public static ClusterHash fromCompactBipartition(CompactBipartition compactBip,
                                                   RangeBipartition.HashFunction hashFunction,
                                                   long[][] prefixSums, long[][] prefixXORs,
                                                   int[][] geneTreeOrderings) {
        // Get the actual taxa from both left and right clusters
        int[] ordering = geneTreeOrderings[compactBip.geneTreeIndex];
        Set<Integer> unionTaxa = new HashSet<>();
        
        // Add taxa from left cluster
        for (int i = compactBip.leftStart; i < compactBip.leftEnd && i < ordering.length; i++) {
            unionTaxa.add(ordering[i]);
        }
        
        // Add taxa from right cluster
        for (int i = compactBip.rightStart; i < compactBip.rightEnd && i < ordering.length; i++) {
            unionTaxa.add(ordering[i]);
        }
        
        // Calculate hash using unified hashing utility
        return fromTaxonSet(unionTaxa, hashFunction);
    }
    
    
    /**
     * Create ClusterHash from a set of taxon IDs.
     * Uses the unified hashing utility to ensure consistency with RangeBipartition hashing.
     */
    public static ClusterHash fromTaxonSet(Set<Integer> taxonIds,
                                         RangeBipartition.HashFunction hashFunction) {
        // Calculate sum and XOR for the taxon set
        long[] sumAndXOR = HashingUtils.calculateTaxonSetSumAndXOR(taxonIds);
        long sum = sumAndXOR[0];
        long xor = sumAndXOR[1];
        int size = taxonIds.size();
        
        // Calculate hash based on the hash function type
        long hashValue;
        if (hashFunction instanceof RangeBipartition.SumHashFunction) {
            // For single cluster, treat as left cluster with empty right
            hashValue = HashingUtils.calculateSingleClusterSumHash(sum, size);
        } else if (hashFunction instanceof RangeBipartition.XORHashFunction) {
            // For single cluster, treat as left cluster with empty right
            hashValue = HashingUtils.calculateSingleClusterXORHash(xor, size);
        } else {
            // Default to sum-based hashing
            hashValue = HashingUtils.calculateSingleClusterSumHash(sum, size);
        }
        
        return new ClusterHash(hashValue, size);
    }
    
    /**
     * Create ClusterHash for the left cluster of a CompactBipartition.
     * Uses unified hashing for consistency with RangeBipartition calculations.
     */
    public static ClusterHash fromLeftCluster(CompactBipartition compactBip,
                                            RangeBipartition.HashFunction hashFunction,
                                            long[][] prefixSums, long[][] prefixXORs,
                                            int[][] geneTreeOrderings) {
        // Option 1: Use range-based calculation (more efficient)
        if (prefixSums != null && prefixXORs != null && 
            compactBip.geneTreeIndex < prefixSums.length && 
            compactBip.geneTreeIndex < prefixXORs.length) {
            
            long[] sumAndXOR = HashingUtils.calculateRangeSumAndXOR(
                prefixSums[compactBip.geneTreeIndex], 
                prefixXORs[compactBip.geneTreeIndex],
                compactBip.leftStart, 
                compactBip.leftEnd
            );
            
            long sum = sumAndXOR[0];
            long xor = sumAndXOR[1];
            int size = compactBip.leftEnd - compactBip.leftStart;
            
            long hashValue;
            if (hashFunction instanceof RangeBipartition.SumHashFunction) {
                hashValue = HashingUtils.calculateSingleClusterSumHash(sum, size);
            } else if (hashFunction instanceof RangeBipartition.XORHashFunction) {
                hashValue = HashingUtils.calculateSingleClusterXORHash(xor, size);
            } else {
                hashValue = HashingUtils.calculateSingleClusterSumHash(sum, size);
            }
            
            return new ClusterHash(hashValue, size);
        }
        
        // Option 2: Fall back to taxon set approach
        int[] ordering = geneTreeOrderings[compactBip.geneTreeIndex];
        Set<Integer> leftTaxa = new HashSet<>();
        
        for (int i = compactBip.leftStart; i < compactBip.leftEnd && i < ordering.length; i++) {
            leftTaxa.add(ordering[i]);
        }
        
        return fromTaxonSet(leftTaxa, hashFunction);
    }
    
    /**
     * Create ClusterHash for the right cluster of a CompactBipartition.
     * Uses unified hashing for consistency with RangeBipartition calculations.
     */
    public static ClusterHash fromRightCluster(CompactBipartition compactBip,
                                             RangeBipartition.HashFunction hashFunction,
                                             long[][] prefixSums, long[][] prefixXORs,
                                             int[][] geneTreeOrderings) {
        // Option 1: Use range-based calculation (more efficient)
        if (prefixSums != null && prefixXORs != null && 
            compactBip.geneTreeIndex < prefixSums.length && 
            compactBip.geneTreeIndex < prefixXORs.length) {
            
            long[] sumAndXOR = HashingUtils.calculateRangeSumAndXOR(
                prefixSums[compactBip.geneTreeIndex], 
                prefixXORs[compactBip.geneTreeIndex],
                compactBip.rightStart, 
                compactBip.rightEnd
            );
            
            long sum = sumAndXOR[0];
            long xor = sumAndXOR[1];
            int size = compactBip.rightEnd - compactBip.rightStart;
            
            long hashValue;
            if (hashFunction instanceof RangeBipartition.SumHashFunction) {
                hashValue = HashingUtils.calculateSingleClusterSumHash(sum, size);
            } else if (hashFunction instanceof RangeBipartition.XORHashFunction) {
                hashValue = HashingUtils.calculateSingleClusterXORHash(xor, size);
            } else {
                hashValue = HashingUtils.calculateSingleClusterSumHash(sum, size);
            }
            
            return new ClusterHash(hashValue, size);
        }
        
        // Option 2: Fall back to taxon set approach
        int[] ordering = geneTreeOrderings[compactBip.geneTreeIndex];
        Set<Integer> rightTaxa = new HashSet<>();
        
        for (int i = compactBip.rightStart; i < compactBip.rightEnd && i < ordering.length; i++) {
            rightTaxa.add(ordering[i]);
        }
        
        return fromTaxonSet(rightTaxa, hashFunction);
    }
    
    /**
     * Get the size (number of taxa) in this cluster.
     */
    public int getSize() {
        return size;
    }
    
    /**
     * Get the hash value for this cluster.
     */
    public Object getHash() {
        return hash;
    }
    
    /**
     * Check if this cluster represents a single taxon.
     */
    public boolean isSingleTaxon() {
        return size == 1;
    }
    
    /**
     * Check if this cluster represents two taxa.
     */
    public boolean isTwoTaxa() {
        return size == 2;
    }
    
    /**
     * Reconstruct the taxon set for this cluster (expensive operation, use sparingly).
     * This is mainly used for tree reconstruction when we need the actual taxon IDs.
     */
    public Set<Integer> reconstructTaxonSet(int[][] geneTreeOrderings, 
                                          CompactBipartition sourceBipartition) {
        if (taxonSet == null) {
            taxonSet = new HashSet<>();
            
            if (sourceBipartition != null) {
                // Add taxa from both left and right clusters
                int[] ordering = geneTreeOrderings[sourceBipartition.geneTreeIndex];
                
                for (int i = sourceBipartition.leftStart; i < sourceBipartition.leftEnd && i < ordering.length; i++) {
                    taxonSet.add(ordering[i]);
                }
                
                for (int i = sourceBipartition.rightStart; i < sourceBipartition.rightEnd && i < ordering.length; i++) {
                    taxonSet.add(ordering[i]);
                }
            }
        }
        
        return taxonSet;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof ClusterHash)) return false;
        ClusterHash other = (ClusterHash) obj;
        
        // Two clusters are equal if they have the same hash and size
        return this.size == other.size && this.hash.equals(other.hash);
    }
    
    @Override
    public int hashCode() {
        return hash.hashCode() ^ (size << 16);
    }
    
    @Override
    public String toString() {
        return "ClusterHash[hash=" + hash + ", size=" + size + "]";
    }
}
