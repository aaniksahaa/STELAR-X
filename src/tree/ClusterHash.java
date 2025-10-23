package tree;

import java.util.Set;
import java.util.HashSet;

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
     * This correctly calculates the hash for the actual union of taxa, not just the range.
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
        
        // Calculate hash using the same method as fromTaxonSet
        return fromTaxonSet(unionTaxa, hashFunction);
    }
    
    /**
     * Create ClusterHash from a range representing a cluster.
     */
    public static ClusterHash fromRange(int geneTreeIndex, int start, int end,
                                      RangeBipartition.HashFunction hashFunction,
                                      long[][] prefixSums, long[][] prefixXORs) {
        // Create range bipartition with the cluster as left and empty right
        RangeBipartition range = new RangeBipartition(geneTreeIndex, start, end, end, end);
        Object hash = hashFunction.calculateHash(range, prefixSums, prefixXORs);
        int size = end - start;
        
        return new ClusterHash(hash, size);
    }
    
    /**
     * Create ClusterHash from a set of taxon IDs.
     * This is used when we need to create cluster hashes from BitSets during transition.
     */
    public static ClusterHash fromTaxonSet(Set<Integer> taxonIds,
                                         RangeBipartition.HashFunction hashFunction) {
        // For taxon sets, we need to create a permutation-invariant hash
        // We'll use a simple approach: sort the taxon IDs and hash them
        int[] sortedTaxa = taxonIds.stream().mapToInt(Integer::intValue).sorted().toArray();
        
        // Calculate hash using the same approach as RangeBipartition hash functions
        long sum = 0;
        long xor = 0;
        for (int taxonId : sortedTaxa) {
            long hashedTaxon = hashSingleTaxon(taxonId);
            sum += hashedTaxon;
            xor ^= hashedTaxon;
        }
        
        // Use the same mixing approach as in RangeBipartition.SumHashFunction
        long combined = sum * 0x9e3779b97f4a7c15L ^ xor * 0xc2b2ae3d27d4eb4fL ^ sortedTaxa.length;
        Object hash = Long.rotateLeft(combined, 27) ^ (combined >>> 33);
        
        return new ClusterHash(hash, taxonIds.size());
    }
    
    private static long hashSingleTaxon(int taxonId) {
        long x = taxonId;
        x ^= x >>> 16;
        x *= 0x85ebca6b;
        x ^= x >>> 13;
        x *= 0xc2b2ae35;
        x ^= x >>> 16;
        return x == 0 ? 1 : x;
    }
    
    /**
     * Create ClusterHash for the left cluster of a CompactBipartition.
     */
    public static ClusterHash fromLeftCluster(CompactBipartition compactBip,
                                            RangeBipartition.HashFunction hashFunction,
                                            long[][] prefixSums, long[][] prefixXORs,
                                            int[][] geneTreeOrderings) {
        // Get the actual taxa from the left cluster
        int[] ordering = geneTreeOrderings[compactBip.geneTreeIndex];
        Set<Integer> leftTaxa = new HashSet<>();
        
        for (int i = compactBip.leftStart; i < compactBip.leftEnd && i < ordering.length; i++) {
            leftTaxa.add(ordering[i]);
        }
        
        return fromTaxonSet(leftTaxa, hashFunction);
    }
    
    /**
     * Create ClusterHash for the right cluster of a CompactBipartition.
     */
    public static ClusterHash fromRightCluster(CompactBipartition compactBip,
                                             RangeBipartition.HashFunction hashFunction,
                                             long[][] prefixSums, long[][] prefixXORs,
                                             int[][] geneTreeOrderings) {
        // Get the actual taxa from the right cluster
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
