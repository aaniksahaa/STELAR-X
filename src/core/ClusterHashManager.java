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
 * 3. Cluster equality verification with fallback to expensive checks
 * 4. Extraction of cluster hashes from bipartitions
 * 
 * The goal is to replace BitSet-based cluster keys in DP with compact hash keys,
 * reducing memory usage from O(n) per cluster to O(1) per cluster.
 */
public class ClusterHashManager {
    
    private final Map<ClusterHashPair, Set<Integer>> hashToTaxonSets;
    private final Map<Set<Integer>, ClusterHashPair> taxonSetToHash;
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
        
        this.hashToTaxonSets = new HashMap<>();
        this.taxonSetToHash = new HashMap<>();
        this.geneTreeOrderings = geneTreeOrderings;
        this.prefixSums = prefixSums;
        this.prefixXORs = prefixXORs;
        this.realTaxaCount = realTaxaCount;
        
        System.out.println("Cluster hash manager initialized");
        System.out.println("==== CLUSTER HASH MANAGER READY ====");
    }
    
    /**
     * Compute hash for a cluster defined by range in a gene tree.
     */
    public ClusterHashPair getClusterHash(int geneTreeIndex, int start, int end) {
        hashComputations++;
        
        ClusterHashPair hash = HashUtils.computeClusterHash(geneTreeIndex, start, end,
                                                           geneTreeOrderings, prefixSums, prefixXORs);
        
        // Cache the taxon set for verification if needed
        Set<Integer> taxonSet = getRangeTaxonSet(geneTreeIndex, start, end);
        cacheClusterMapping(hash, taxonSet);
        
        return hash;
    }
    
    /**
     * Compute hash for existing BitSet cluster (for compatibility).
     * This is the main interface for converting existing DP code.
     */
    public ClusterHashPair getClusterHash(BitSet cluster) {
        // Convert BitSet to taxon set
        Set<Integer> taxonSet = bitSetToTaxonSet(cluster);
        
        // Check if we already have a hash for this taxon set
        ClusterHashPair cachedHash = taxonSetToHash.get(taxonSet);
        if (cachedHash != null) {
            cacheHits++;
            return cachedHash;
        }
        
        hashComputations++;
        
        // Try to find a range representation for this cluster
        ClusterHashPair rangeHash = findRangeBasedHash(taxonSet);
        if (rangeHash != null) {
            cacheClusterMapping(rangeHash, taxonSet);
            return rangeHash;
        }
        
        // Fallback to deterministic hash based on taxon set content
        ClusterHashPair fallbackHash = HashUtils.computeFallbackClusterHash(taxonSet);
        cacheClusterMapping(fallbackHash, taxonSet);
        
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
        
        // Cache this result
        Set<Integer> unionSet = getRangeTaxonSet(bip.geneTreeIndex, bip.leftStart, bip.rightEnd);
        cacheClusterMapping(unionHash, unionSet);
        
        return unionHash;
    }
    
    /**
     * Verify cluster equality with fallback to expensive verification.
     * This is the key method for trusting double hashes vs. doing expensive checks.
     */
    public boolean clustersEqual(ClusterHashPair hash1, ClusterHashPair hash2) {
        if (hash1.equals(hash2)) {
            return true; // Trust the double hash - collision probability is extremely low
        }
        
        // Hash collision detected - need expensive verification
        hashCollisions++;
        expensiveVerifications++;
        
        System.out.println("Hash collision detected between " + hash1.toDebugString() + 
                         " and " + hash2.toDebugString() + " - performing expensive verification");
        
        Set<Integer> set1 = hashToTaxonSets.get(hash1);
        Set<Integer> set2 = hashToTaxonSets.get(hash2);
        
        if (set1 != null && set2 != null) {
            boolean equal = set1.equals(set2);
            System.out.println("Expensive verification result: " + equal);
            return equal;
        }
        
        System.out.println("Cannot verify - taxon sets not available, assuming different");
        return false; // Cannot verify - assume different
    }
    
    /**
     * Get the taxon set size for a cluster hash.
     * This is needed for DP operations that require cluster size.
     */
    public int getClusterSize(ClusterHashPair clusterHash) {
        Set<Integer> taxonSet = hashToTaxonSets.get(clusterHash);
        if (taxonSet != null) {
            return taxonSet.size();
        }
        
        // If not cached, we can't determine size without expensive lookup
        System.err.println("Warning: Cluster size not available for hash " + clusterHash.toDebugString());
        return -1; // Indicate unknown size
    }
    
    /**
     * Get the actual taxon set for a cluster hash (for debugging/verification).
     */
    public Set<Integer> getTaxonSet(ClusterHashPair clusterHash) {
        return hashToTaxonSets.get(clusterHash);
    }
    
    /**
     * Cache the mapping between cluster hash and taxon set.
     */
    private void cacheClusterMapping(ClusterHashPair hash, Set<Integer> taxonSet) {
        // Check for hash collisions
        Set<Integer> existingSet = hashToTaxonSets.get(hash);
        if (existingSet != null && !existingSet.equals(taxonSet)) {
            hashCollisions++;
            System.err.println("Hash collision detected! Hash " + hash.toDebugString() + 
                             " maps to both " + existingSet + " and " + taxonSet);
        }
        
        hashToTaxonSets.put(hash, new HashSet<>(taxonSet)); // Defensive copy
        taxonSetToHash.put(new HashSet<>(taxonSet), hash);   // Defensive copy
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
     * Create a BitSet from a taxon set (for compatibility with existing code).
     */
    public BitSet createBitSetFromHash(ClusterHashPair clusterHash) {
        Set<Integer> taxonSet = hashToTaxonSets.get(clusterHash);
        if (taxonSet == null) {
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
                System.out.println("Precomputed hashes for " + processedGroups + "/" + 
                                 hashToBipartitions.size() + " bipartition groups");
            }
        }
        
        System.out.println("Precomputation completed: " + totalClusters + " cluster hashes computed");
        System.out.println("Cache contains " + hashToTaxonSets.size() + " unique cluster hashes");
    }
    
    /**
     * Get statistics about cluster hash operations.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Cluster Hash Manager Statistics:\n");
        sb.append("  Hash computations: ").append(hashComputations).append("\n");
        sb.append("  Cache hits: ").append(cacheHits).append("\n");
        sb.append("  Hash collisions: ").append(hashCollisions).append("\n");
        sb.append("  Expensive verifications: ").append(expensiveVerifications).append("\n");
        sb.append("  Cached cluster hashes: ").append(hashToTaxonSets.size()).append("\n");
        
        if (hashComputations > 0) {
            sb.append("  Cache hit rate: ").append(String.format("%.2f%%", 
                100.0 * cacheHits / (hashComputations + cacheHits))).append("\n");
            sb.append("  Collision rate: ").append(String.format("%.6f%%", 
                100.0 * hashCollisions / hashComputations)).append("\n");
        }
        
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
     * Clear all cached mappings (for memory management).
     */
    public void clearCache() {
        System.out.println("Clearing cluster hash cache...");
        int cachedCount = hashToTaxonSets.size();
        
        hashToTaxonSets.clear();
        taxonSetToHash.clear();
        
        System.out.println("Cleared " + cachedCount + " cached cluster mappings");
    }
    
    /**
     * Get hash for all taxa cluster (entire range in first gene tree).
     */
    public ClusterHashPair getAllTaxaClusterHash() {
        // Use the entire range of the first gene tree (0 to numTaxa)
        return getClusterHash(0, 0, realTaxaCount);
    }
}
