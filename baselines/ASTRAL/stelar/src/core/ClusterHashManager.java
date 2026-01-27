package core;

import java.util.*;
import tree.ClusterHashPair;
import tree.RangeBipartition;
import utils.BitSet;
import utils.HashUtils;

/**
 * Manages cluster hash computation and lookup for memory-optimized DP operations.
 * 
 * This class provides:
 * 1. Double hash computation for clusters (sum + XOR hashes)
 * 2. Mapping between BitSet clusters and hash representations
 * 3. Optional cluster equality verification (can be disabled for performance)
 * 4. Extraction of cluster hashes from bipartitions
 * 
 * The goal is to replace BitSet-based cluster keys in DP with compact hash keys,
 * reducing memory usage from O(n) per cluster to O(1) per cluster.
 */
public class ClusterHashManager {
    
    // Lightweight cluster range information
    public static class ClusterRange {
        public final int geneTreeIndex;
        public final int start;
        public final int end;
        
        public ClusterRange(int geneTreeIndex, int start, int end) {
            this.geneTreeIndex = geneTreeIndex;
            this.start = start;
            this.end = end;
        }
        
        public int size() {
            return end - start;
        }
        
        @Override
        public String toString() {
            return "ClusterRange[tree=" + geneTreeIndex + ", range=[" + start + "," + end + ")]";
        }
    }
    
    // Lightweight mapping: hash -> cluster range (just 3 integers)
    private final Map<ClusterHashPair, ClusterRange> hashToRange;
    
    private final int[][] geneTreeOrderings;
    private final long[][] prefixSums;
    private final long[][] prefixXORs;
    private final int realTaxaCount;
    
    // Statistics for monitoring
    private long hashComputations = 0;
    private long hashCollisions = 0;
    private long expensiveVerifications = 0;
    private long cacheHits = 0;
    
    public ClusterHashManager(int[][] geneTreeOrderings, long[][] prefixSums, 
                             long[][] prefixXORs, int realTaxaCount) {
        System.out.println("==== INITIALIZING CLUSTER HASH MANAGER ====");
        System.out.println("Real taxa count: " + realTaxaCount);
        System.out.println("Gene trees: " + (geneTreeOrderings != null ? geneTreeOrderings.length : 0));
        System.out.println("Storage mode: Lightweight range mapping (3 integers per cluster)");
        
        this.hashToRange = new HashMap<>();
        this.geneTreeOrderings = geneTreeOrderings;
        this.prefixSums = prefixSums;
        this.prefixXORs = prefixXORs;
        this.realTaxaCount = realTaxaCount;
        
        System.out.println("Cluster hash manager initialized (lightweight range mode)");
        System.out.println("==== CLUSTER HASH MANAGER READY ====");
    }
    
    /**
     * Compute hash for a cluster defined by range in a gene tree.
     */
    public ClusterHashPair getClusterHash(int geneTreeIndex, int start, int end) {
        hashComputations++;
        
        ClusterHashPair hash = HashUtils.computeClusterHash(geneTreeIndex, start, end,
                                                           geneTreeOrderings, prefixSums, prefixXORs);
        
        // Store lightweight range mapping for O(1) size queries and on-demand taxon derivation
        hashToRange.put(hash, new ClusterRange(geneTreeIndex, start, end));
        
        return hash;
    }
    
    /**
     * Compute hash for existing BitSet cluster (for compatibility).
     * This is the main interface for converting existing DP code.
     */
    public ClusterHashPair getClusterHash(BitSet cluster) {
        hashComputations++;
        
        // Convert BitSet to taxon set
        Set<Integer> taxonSet = bitSetToTaxonSet(cluster);
        
        // Try to find a range representation for this cluster
        ClusterHashPair rangeHash = findRangeBasedHash(taxonSet);
        if (rangeHash != null) {
            return rangeHash;
        }
        
        // Fallback to deterministic hash based on taxon set content
        ClusterHashPair fallbackHash = HashUtils.computeFallbackClusterHash(taxonSet);
        
        // Store a dummy range for fallback hashes (we can't derive range info from them)
        // This is mainly for compatibility - ideally all clusters should be range-based
        hashToRange.put(fallbackHash, new ClusterRange(-1, 0, taxonSet.size()));
        
        return fallbackHash;
    }
    
    /**
     * Get left cluster hash from a bipartition.
     */
    public ClusterHashPair getLeftClusterHash(RangeBipartition bip) {
        ClusterHashPair hash = getClusterHash(bip.geneTreeIndex, bip.leftStart, bip.leftEnd);
        
        // System.out.println("Left cluster hash for bipartition " + bip + ": " + hash.toDebugString());
        
        return hash;
    }
    
    /**
     * Get right cluster hash from a bipartition.
     */
    public ClusterHashPair getRightClusterHash(RangeBipartition bip) {
        ClusterHashPair hash = getClusterHash(bip.geneTreeIndex, bip.rightStart, bip.rightEnd);
        
        // System.out.println("Right cluster hash for bipartition " + bip + ": " + hash.toDebugString());
        
        return hash;
    }
    
    /**
     * Get union cluster hash from a bipartition (left ∪ right).
     */
    public ClusterHashPair getUnionClusterHash(RangeBipartition bip) {
        // By definition, subtree bipartitions must have contiguous ranges
        // The union is simply [leftStart, rightEnd]
        
        if (bip.rightStart != bip.leftEnd) {
            throw new IllegalArgumentException(
                "Non-contiguous bipartition detected in tree " + bip.geneTreeIndex + 
                ": left=[" + bip.leftStart + "," + bip.leftEnd + "], right=[" + bip.rightStart + "," + bip.rightEnd + "]. " +
                "Subtree bipartitions must have contiguous ranges where rightStart == leftEnd."
            );
        }
        
        // Contiguous ranges - union is simply [leftStart, rightEnd]
        hashComputations++;
        ClusterHashPair unionHash = HashUtils.computeClusterHash(bip.geneTreeIndex, 
                                                                bip.leftStart, bip.rightEnd,
                                                                geneTreeOrderings, prefixSums, prefixXORs);
        
        // Store the range mapping for the union cluster
        hashToRange.put(unionHash, new ClusterRange(bip.geneTreeIndex, bip.leftStart, bip.rightEnd));
        
        return unionHash;
    }
    
    /**
     * Verify cluster equality by trusting the hash completely.
     * In lightweight mode, we don't do expensive verification.
     */
    public boolean clustersEqual(ClusterHashPair hash1, ClusterHashPair hash2) {
        // Simply trust the hash - different hashes = different clusters
        return hash1.equals(hash2);
    }
    
    /**
     * Get the taxon set size for a cluster hash in O(1) time.
     * This is needed for DP operations that require cluster size.
     */
    public int getClusterSize(ClusterHashPair clusterHash) {
        ClusterRange range = hashToRange.get(clusterHash);
        if (range != null) {
            return range.size(); // O(1) calculation: end - start
        }
        
        System.err.println("Warning: Cluster size not available for hash " + clusterHash.toDebugString());
        return -1; // Unknown size
    }
    
    /**
     * Get the actual taxon set for a cluster hash by deriving it from gene tree orderings.
     * This is used for tree reconstruction and debugging.
     */
    public Set<Integer> getTaxonSet(ClusterHashPair clusterHash) {
        ClusterRange range = hashToRange.get(clusterHash);
        if (range != null) {
            return getRangeTaxonSet(range.geneTreeIndex, range.start, range.end);
        }
        
        System.err.println("Warning: Cannot derive taxon set for hash " + clusterHash.toDebugString());
        return new HashSet<>();
    }
    
    // Removed findRangeBasedHashOptimized - not needed for contiguous subtree bipartitions
    
    /**
     * Try to find a range-based hash for the given taxon set.
     * This attempts to find a gene tree range that matches the taxon set.
     */
    private ClusterHashPair findRangeBasedHash(Set<Integer> taxonSet) {
        // This is an expensive operation - search through gene tree orderings
        // to find a range that matches the taxon set
        System.out.println("WARNING: Using expensive O(n³) findRangeBasedHash for taxon set: " + taxonSet);
        
        if (geneTreeOrderings == null) {
            return null;
        }
        
        for (int treeIdx = 0; treeIdx < geneTreeOrderings.length; treeIdx++) {
            int[] ordering = geneTreeOrderings[treeIdx];
            if (ordering == null) continue;
            
            // Try to find a contiguous range that matches the taxon set
            for (int start = 0; start < ordering.length; start++) {
                for (int end = start + 1; end <= ordering.length; end++) {
                    Set<Integer> rangeSet = getRangeTaxonSet(treeIdx, start, end);
                    if (rangeSet.equals(taxonSet)) {
                        // Found a matching range - compute its hash
                        return HashUtils.computeClusterHash(treeIdx, start, end,
                                                           geneTreeOrderings, prefixSums, prefixXORs);
                    }
                    
                    // Early termination if range is already larger than target
                    if (rangeSet.size() > taxonSet.size()) {
                        break;
                    }
                }
            }
        }
        
        return null; // No range representation found
    }
    
    /**
     * Get taxon set for a range in a gene tree.
     */
    private Set<Integer> getRangeTaxonSet(int geneTreeIndex, int start, int end) {
        Set<Integer> taxonSet = new HashSet<>();
        
        if (geneTreeOrderings != null && geneTreeIndex >= 0 && geneTreeIndex < geneTreeOrderings.length) {
            int[] ordering = geneTreeOrderings[geneTreeIndex];
            if (ordering != null) {
                for (int i = start; i < end && i < ordering.length; i++) {
                    taxonSet.add(ordering[i]);
                }
            }
        }
        
        return taxonSet;
    }
    
    /**
     * Convert BitSet to set of taxon IDs.
     */
    private Set<Integer> bitSetToTaxonSet(BitSet bitSet) {
        Set<Integer> taxonSet = new HashSet<>();
        for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i + 1)) {
            taxonSet.add(i);
        }
        return taxonSet;
    }
    
    /**
     * Create a BitSet from a cluster hash by deriving the taxon set on-demand.
     */
    public BitSet createBitSetFromHash(ClusterHashPair clusterHash) {
        Set<Integer> taxonSet = getTaxonSet(clusterHash);
        if (taxonSet == null || taxonSet.isEmpty()) {
            System.err.println("Warning: Cannot create BitSet for hash " + clusterHash.toDebugString() + 
                             " - taxon set not available");
            return new BitSet(realTaxaCount);
        }
        
        BitSet bitSet = new BitSet(realTaxaCount);
        for (int taxonId : taxonSet) {
            if (taxonId >= 0 && taxonId < realTaxaCount) {
                bitSet.set(taxonId);
            }
        }
        
        return bitSet;
    }
    
    /**
     * Precompute hashes for all possible clusters from the given bipartitions.
     * This can improve performance by avoiding repeated hash computations.
     */
    public void precomputeClusterHashes(Map<Object, List<RangeBipartition>> hashToBipartitions) {
        System.out.println("Precomputing cluster hashes from range bipartitions...");
        
        int processedGroups = 0;
        int totalClusters = 0;
        
        for (Map.Entry<Object, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
            List<RangeBipartition> ranges = entry.getValue();
            
            for (RangeBipartition range : ranges) {
                // Precompute hashes for left, right, and union clusters
                getLeftClusterHash(range);
                getRightClusterHash(range);
                getUnionClusterHash(range);
                totalClusters += 3;
            }
            
            processedGroups++;
            
            // Log progress for large datasets
            if (processedGroups % 100 == 0 || processedGroups == hashToBipartitions.size()) {
                // System.out.println("Precomputed hashes for " + processedGroups + "/" + 
                //                 hashToBipartitions.size() + " bipartition groups");
            }
        }
        
        System.out.println("Precomputation completed: " + totalClusters + " cluster hashes computed");
        System.out.println("Cache contains " + hashToRange.size() + " lightweight cluster range mappings");
    }
    
    /**
     * Get statistics about cluster hash operations.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Cluster Hash Manager Statistics:\n");
        sb.append("  Hash computations: ").append(hashComputations).append("\n");
        sb.append("  Cache hits: ").append(cacheHits).append("\n");
        sb.append("  Hash collisions: 0 (trusting hash completely)\n");
        sb.append("  Expensive verifications: 0 (trusting hash completely)\n");
        sb.append("  Cached cluster ranges: ").append(hashToRange.size()).append("\n");
        
        if (hashComputations > 0) {
            sb.append("  Cache hit rate: ").append(String.format("%.2f%%", 
                100.0 * cacheHits / (hashComputations + cacheHits))).append("\n");
        }
        sb.append("  Collision rate: 0.000000% (trusting hash completely)\n");
        
        // Calculate memory usage
        long memoryBytes = hashToRange.size() * (64 + 12); // Rough estimate: hash object + 3 ints
        sb.append("  Estimated memory usage: ").append(memoryBytes / 1024).append(" KB for range mappings\n");
        
        return sb.toString();
    }
    
    /**
     * Reset statistics counters.
     */
    public void resetStatistics() {
        hashComputations = 0;
        hashCollisions = 0;
        expensiveVerifications = 0;
        cacheHits = 0;
    }
    
    /**
     * Clear all cached range mappings (for memory management).
     */
    public void clearCache() {
        System.out.println("Clearing cluster range cache...");
        int cachedCount = hashToRange.size();
        
        hashToRange.clear();
        
        System.out.println("Cleared " + cachedCount + " cached cluster range mappings");
    }
    
    /**
     * Get hash for all taxa cluster (root of a complete gene tree).
     * 
     * Find a gene tree that contains all taxa in the dataset and use its root cluster.
     */
    public ClusterHashPair getAllTaxaClusterHash() {
        System.out.println("==== COMPUTING ALL TAXA CLUSTER HASH ====");
        System.out.println("Real taxa count (dataset): " + realTaxaCount);
        
        // Find a gene tree that contains ALL taxa (complete tree)
        int completeTreeIndex = -1;
        
        for (int treeIdx = 0; treeIdx < geneTreeOrderings.length; treeIdx++) {
            int taxaInThisTree = geneTreeOrderings[treeIdx].length;
            
            // Log first few trees for debugging
            if (treeIdx < 5) {
                System.out.println("Tree " + treeIdx + " has " + taxaInThisTree + " taxa");
            }
            
            // Check if this tree contains all taxa
            if (taxaInThisTree == realTaxaCount) {
                completeTreeIndex = treeIdx;
                System.out.println("Found complete tree at index " + treeIdx + " with all " + realTaxaCount + " taxa");
                break; // Use the first complete tree we find
            }
        }
        
        if (completeTreeIndex == -1) {
            // Fallback: find the tree with the most taxa
            int maxTaxaInTree = 0;
            for (int treeIdx = 0; treeIdx < geneTreeOrderings.length; treeIdx++) {
                int taxaInThisTree = geneTreeOrderings[treeIdx].length;
                if (taxaInThisTree > maxTaxaInTree) {
                    maxTaxaInTree = taxaInThisTree;
                    completeTreeIndex = treeIdx;
                }
            }
            
            System.out.println("WARNING: No tree contains all " + realTaxaCount + " taxa!");
            System.out.println("Using tree " + completeTreeIndex + " with " + maxTaxaInTree + " taxa as fallback");
            System.out.println("This may cause DP issues - the algorithm expects at least one complete tree");
        }
        
        // Use the entire range of the complete tree (represents the root cluster)
        int numTaxaInCompleteTree = geneTreeOrderings[completeTreeIndex].length;
        ClusterHashPair result = getClusterHash(completeTreeIndex, 0, numTaxaInCompleteTree);
        
        System.out.println("All taxa cluster hash from tree " + completeTreeIndex + ": " + result.toDebugString());
        System.out.println("Cluster size: " + numTaxaInCompleteTree + " taxa");
        System.out.println("==== ALL TAXA CLUSTER HASH COMPUTED ====");
        
        return result;
    }
}
