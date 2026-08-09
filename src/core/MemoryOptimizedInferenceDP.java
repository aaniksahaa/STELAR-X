package core;

import java.util.*;
import preprocessing.GeneTrees;
import tree.*;

/**
 * Memory-optimized inference DP that uses cluster hashes instead of BitSet keys.
 * 
 * This implementation replaces:
 * - Map<BitSet, ...> with Map<ClusterHashPair, ...> for O(1) vs O(n) memory per key
 * - BitSet cluster operations with hash-based operations where possible
 * - Expensive BitSet equality with fast hash equality (with fallback verification)
 * 
 * The DP algorithm remains the same, but memory usage is significantly reduced
 * for large numbers of taxa and clusters.
 */
public class MemoryOptimizedInferenceDP {
    
    private final GeneTrees geneTrees;
    private final List<RangeBipartition> candidateRangeBips;
    private final Map<ClusterHashPair, List<RangeBipartition>> clusterHashToRangeBips;
    private Map<RangeBipartition, Double> rangeBipWeights;
    private final Map<ClusterHashPair, Double> dpMemo;
    private final Map<ClusterHashPair, Object> dpChoice;  // Can be RangeBipartition or MixedBipartition
    private final ClusterHashManager clusterHashManager;
    private final MemoryEfficientBipartitionManager bipartitionManager;
    
    // Mixed bipartition support (cross-tree recombination)
    private List<MixedBipartition> mixedBipartitions;
    private Map<ClusterHashPair, List<MixedBipartition>> clusterHashToMixedBips;
    private Map<MixedBipartition, Double> mixedBipWeights;
    private boolean useMixedBipartitions = false;
    
    // Statistics for monitoring
    private long dpCalls = 0;
    private long memoHits = 0;
    private long clusterValidations = 0;
    private long totalProcessingTime = 0;
    private long mixedBipChoices = 0;  // Count how many times mixed bipartition was chosen
    
    public MemoryOptimizedInferenceDP(GeneTrees geneTrees, List<RangeBipartition> candidateRangeBips) {
        System.out.println("==== INITIALIZING MEMORY-OPTIMIZED INFERENCE DP ====");
        System.out.println("Number of candidate bipartitions: " + candidateRangeBips.size());
        System.out.println("Real taxa count: " + geneTrees.realTaxaCount);
        
        this.geneTrees = geneTrees;
        this.candidateRangeBips = candidateRangeBips;
        this.clusterHashToRangeBips = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
        
        // Reuse the range data already built during GeneTrees preprocessing.
        System.out.println("Reusing preprocessed MemoryEfficientBipartitionManager...");
        this.bipartitionManager = geneTrees.getBipartitionManager();
        
        System.out.println("Initializing ClusterHashManager...");
        this.clusterHashManager = new ClusterHashManager(
            bipartitionManager.getGeneTreeTaxaOrdering(),
            bipartitionManager.getPrefixSums(),
            bipartitionManager.getPrefixXORs(),
            geneTrees.realTaxaCount);
        
        // Precompute cluster hashes for better performance
        System.out.println("Precomputing cluster hashes...");
        clusterHashManager.precomputeClusterHashes(bipartitionManager.getHashToBipartitions());
        
        System.out.println("Preprocessing candidates with cluster hashes...");
        preprocessCandidatesWithHashes();
        
        System.out.println("Calculating weights using memory-optimized approach...");
        calculateWeights();
        
        System.out.println("==== MEMORY-OPTIMIZED INFERENCE DP READY ====");
        System.out.println("Cluster hash groups: " + clusterHashToRangeBips.size());
        System.out.println("Cached cluster hashes: " + clusterHashManager.getStatistics().split("\\n")[4]);
    }
    
    /**
     * Enable mixed bipartitions (cross-tree recombination) for extended candidate set.
     * This should be called before solve() if mixed bipartitions are to be used.
     * 
     * @param mixedBips List of mixed bipartitions from CandidateExtender
     */
    public void enableMixedBipartitions(List<MixedBipartition> mixedBips) {
        if (mixedBips == null || mixedBips.isEmpty()) {
            System.out.println("No mixed bipartitions provided, using only gene tree bipartitions");
            return;
        }
        
        System.out.println("\n==== ENABLING MIXED BIPARTITIONS ====");
        System.out.println("Mixed bipartitions to add: " + mixedBips.size());
        
        this.mixedBipartitions = mixedBips;
        this.useMixedBipartitions = true;
        this.clusterHashToMixedBips = new HashMap<>();
        
        // Preprocess mixed bipartitions to create cluster hash mappings
        preprocessMixedBipartitions();
        
        // Calculate weights for mixed bipartitions
        calculateMixedWeights();
        
        System.out.println("Mixed bipartitions enabled");
        System.out.println("Mixed cluster hash groups: " + clusterHashToMixedBips.size());
        
        // Debug: Print first few cluster hashes from mixed bipartitions
        System.out.println("\nDEBUG: Sample cluster hashes in clusterHashToMixedBips:");
        int count = 0;
        for (ClusterHashPair hash : clusterHashToMixedBips.keySet()) {
            System.out.println("  Mixed cluster hash: " + hash.toDebugString());
            if (++count >= 5) break;
        }
        
        // Debug: Compare with RangeBipartition cluster hashes
        System.out.println("\nDEBUG: Sample cluster hashes in clusterHashToRangeBips:");
        count = 0;
        for (ClusterHashPair hash : clusterHashToRangeBips.keySet()) {
            System.out.println("  Range cluster hash: " + hash.toDebugString());
            if (++count >= 5) break;
        }
        
        // Check overlap between the two maps
        int overlapCount = 0;
        for (ClusterHashPair mixedHash : clusterHashToMixedBips.keySet()) {
            if (clusterHashToRangeBips.containsKey(mixedHash)) {
                overlapCount++;
            }
        }
        System.out.println("\nDEBUG: Overlap between mixed and range cluster hashes: " + overlapCount + 
                         " / " + clusterHashToMixedBips.size() + " mixed clusters");
        
        System.out.println("==== MIXED BIPARTITIONS READY ====\n");
    }
    
    /**
     * Preprocess mixed bipartitions to create hash-based cluster mappings.
     */
    private void preprocessMixedBipartitions() {
        System.out.println("Creating hash-based mappings for mixed bipartitions...");
        
        int processed = 0;
        int uniqueClusters = 0;
        
        for (MixedBipartition mixed : mixedBipartitions) {
            // Compute union cluster hash for the mixed bipartition
            // Union = left side ∪ right side
            ClusterHashPair unionHash = getMixedUnionHash(mixed);
            
            // Add to hash-based mapping
            List<MixedBipartition> candidates = clusterHashToMixedBips.get(unionHash);
            if (candidates == null) {
                candidates = new ArrayList<>();
                clusterHashToMixedBips.put(unionHash, candidates);
                uniqueClusters++;
            }
            candidates.add(mixed);
            
            processed++;
        }
        
        System.out.println("Mixed bipartition preprocessing completed");
        System.out.println("Created mappings for " + uniqueClusters + " unique cluster hashes");
    }
    
    /**
     * Get the union cluster hash for a mixed bipartition.
     * Since sides may be from different trees, we compute hash from raw prefix values.
     * 
     * IMPORTANT: We cannot simply add mixed hash values because mixing functions
     * are not additive. Instead, we must:
     * 1. Compute raw sums/XORs for each side from prefix arrays
     * 2. Add raw sums, XOR raw XORs to get union's raw values
     * 3. Apply mixing functions to get final ClusterHashPair
     */
    private ClusterHashPair getMixedUnionHash(MixedBipartition mixed) {
        // Get raw prefix arrays from bipartition manager
        long[][] prefixSums = bipartitionManager.getPrefixSums();
        long[][] prefixXORs = bipartitionManager.getPrefixXORs();
        
        // Compute raw (unmixed) sums for left and right sides
        long leftRawSum = rangeSum(prefixSums[mixed.leftTreeIndex], mixed.leftStart, mixed.leftEnd);
        long rightRawSum = rangeSum(prefixSums[mixed.rightTreeIndex], mixed.rightStart, mixed.rightEnd);
        long unionRawSum = leftRawSum + rightRawSum;
        
        // Compute raw (unmixed) XORs for left and right sides
        long leftRawXor = rangeXOR(prefixXORs[mixed.leftTreeIndex], mixed.leftStart, mixed.leftEnd);
        long rightRawXor = rangeXOR(prefixXORs[mixed.rightTreeIndex], mixed.rightStart, mixed.rightEnd);
        long unionRawXor = leftRawXor ^ rightRawXor;
        
        // Total size for mixing
        int unionSize = mixed.totalSize();
        
        // Apply mixing functions (same as HashUtils.computeClusterHash)
        long mixedSum = mixSumHash(unionRawSum, unionSize);
        long mixedXor = mixXorHash(unionRawXor, unionSize);
        
        return new ClusterHashPair(mixedSum, mixedXor);
    }
    
    // Helper: Compute range sum from prefix array
    private long rangeSum(long[] prefix, int start, int end) {
        if (end <= 0 || prefix.length == 0) return 0;
        int safeEnd = Math.min(end - 1, prefix.length - 1);
        long endSum = (safeEnd >= 0) ? prefix[safeEnd] : 0;
        long startSum = (start > 0 && start - 1 < prefix.length) ? prefix[start - 1] : 0;
        return endSum - startSum;
    }
    
    // Helper: Compute range XOR from prefix array
    private long rangeXOR(long[] prefixXOR, int start, int end) {
        if (prefixXOR == null || end <= 0 || prefixXOR.length == 0) return 0L;
        int safeEnd = Math.min(end - 1, prefixXOR.length - 1);
        long endXOR = (safeEnd >= 0) ? prefixXOR[safeEnd] : 0L;
        long startXOR = (start > 0 && start - 1 < prefixXOR.length) ? prefixXOR[start - 1] : 0L;
        return endXOR ^ startXOR;
    }
    
    // Helper: Apply sum hash mixing (same as HashUtils)
    private static final long MIX_CONST1 = 0x9e3779b97f4a7c15L;
    private long mixSumHash(long rangeSum, int size) {
        long combined = rangeSum * MIX_CONST1 ^ ((long) size << 16);
        return Long.rotateLeft(combined, 27) ^ (combined >>> 33);
    }
    
    // Helper: Apply XOR hash mixing (same as HashUtils)
    private static final long XOR_MIX_A = 0xbf58476d1ce4e5b9L;
    private long mixXorHash(long rangeXOR, int size) {
        long v = rangeXOR ^ (((long) size) << 32) ^ XOR_MIX_A;
        return splitMix64(v);
    }
    
    // SplitMix64 mixing function
    private long splitMix64(long z) {
        z += 0x9e3779b97f4a7c15L;
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }
    
    /**
     * Calculate weights for mixed bipartitions using the weight calculator.
     */
    private void calculateMixedWeights() {
        System.out.println("Calculating weights for mixed bipartitions...");
        
        MemoryOptimizedWeightCalculator calculator = new MemoryOptimizedWeightCalculator(geneTrees);
        mixedBipWeights = calculator.calculateMixedWeights(mixedBipartitions);
        
        // Print weight statistics
        if (!mixedBipWeights.isEmpty()) {
            double minWeight = mixedBipWeights.values().stream().mapToDouble(Double::doubleValue).min().orElse(0.0);
            double maxWeight = mixedBipWeights.values().stream().mapToDouble(Double::doubleValue).max().orElse(0.0);
            double avgWeight = mixedBipWeights.values().stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
            
            System.out.println("Mixed weight statistics: min=" + minWeight + ", max=" + maxWeight + ", avg=" + avgWeight);
        }
    }
    
    /**
     * Preprocess candidates to create hash-based cluster mappings.
     * This replaces the BitSet-based preprocessing in the original implementation.
     */
    private void preprocessCandidatesWithHashes() {
        System.out.println("Creating hash-based cluster mappings for candidates...");
        
        int processedCandidates = 0;
        int uniqueClusters = 0;
        
        for (RangeBipartition rangeBip : candidateRangeBips) {
            // Compute union cluster hash (left ∪ right) for range bipartition
            // For a RangeBipartition, the union is the entire range from leftStart to rightEnd
            ClusterHashPair unionHash = clusterHashManager.getUnionClusterHash(rangeBip);
            
            // Add to hash-based mapping
            List<RangeBipartition> candidates = clusterHashToRangeBips.get(unionHash);
            if (candidates == null) {
                candidates = new ArrayList<>();
                clusterHashToRangeBips.put(unionHash, candidates);
                uniqueClusters++;
            }
            candidates.add(rangeBip);
            
            processedCandidates++;
            
            // // Log progress for large datasets
            // if (processedCandidates % 1000 == 0 || processedCandidates == candidateRangeBips.size()) {
            //     System.out.println("Processed " + processedCandidates + "/" + candidateRangeBips.size() + 
            //                      " candidates, unique clusters: " + uniqueClusters);
            // }
        }
        
        System.out.println("Hash-based preprocessing completed");
        System.out.println("Created mappings for " + uniqueClusters + " unique cluster hashes");
        System.out.println("Average candidates per cluster: " + 
                         (uniqueClusters > 0 ? candidateRangeBips.size() / (double) uniqueClusters : 0));
    }
    
    /**
     * Calculate weights using memory-optimized weight calculator.
     */
    private void calculateWeights() {
        MemoryOptimizedWeightCalculator calculator = new MemoryOptimizedWeightCalculator(geneTrees);
        rangeBipWeights = calculator.calculateWeights(candidateRangeBips);
        
        System.out.println("Weight calculation completed for " + rangeBipWeights.size() + " bipartitions");
        
        // Print weight statistics
        double minWeight = rangeBipWeights.values().stream().mapToDouble(Double::doubleValue).min().orElse(0.0);
        double maxWeight = rangeBipWeights.values().stream().mapToDouble(Double::doubleValue).max().orElse(0.0);
        double avgWeight = rangeBipWeights.values().stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        
        System.out.println("Weight statistics: min=" + minWeight + ", max=" + maxWeight + ", avg=" + avgWeight);
    }
    
    /**
     * Solve the DP problem using hash-based cluster keys.
     */
    public double solve() {
        System.out.println("==== STARTING MEMORY-OPTIMIZED DP SOLUTION ====");
        
        long startTime = System.currentTimeMillis();
        
        // Create hash for all taxa cluster (entire range in first gene tree)
        ClusterHashPair allTaxaHash = clusterHashManager.getAllTaxaClusterHash();
        System.out.println("All taxa cluster hash: " + allTaxaHash.toDebugString());
        System.out.println("All taxa cluster size: " + clusterHashManager.getClusterSize(allTaxaHash));
        
        double result = dp(allTaxaHash);
        
        long endTime = System.currentTimeMillis();
        totalProcessingTime = endTime - startTime;
        
        System.out.println("==== DP SOLUTION COMPLETED ====");
        System.out.println("Optimal score: " + result);
        System.out.println("Processing time: " + totalProcessingTime + " ms");
        System.out.println("DP calls: " + dpCalls);
        System.out.println("Memo hits: " + memoHits + " (" + 
                         (dpCalls > 0 ? String.format("%.2f%%", 100.0 * memoHits / dpCalls) : "0%") + ")");
        System.out.println("Cluster validations: " + clusterValidations);
        System.out.println("Unique clusters in DP state space: " + clusterHashToRangeBips.size());
        System.out.println("Unique clusters visited (memoized): " + dpMemo.size());
        
        // Print mixed bipartition statistics if enabled
        if (useMixedBipartitions) {
            System.out.println("Mixed bipartitions enabled: YES");
            System.out.println("Mixed cluster hash groups: " + 
                             (clusterHashToMixedBips != null ? clusterHashToMixedBips.size() : 0));
            System.out.println("Mixed bipartition choices in optimal tree: " + mixedBipChoices);
        }
        
        // Print cluster hash manager statistics
        System.out.println("\n" + clusterHashManager.getStatistics());
        
        return result;
    }
    
    /**
     * Main DP function using cluster hashes instead of BitSets.
     * This is the core optimization - replacing O(n) BitSet keys with O(1) hash keys.
     * 
     * Now also considers MixedBipartition candidates if enabled.
     */
    private double dp(ClusterHashPair clusterHash) {
        dpCalls++;
        
        // Check memoization
        if (dpMemo.containsKey(clusterHash)) {
            memoHits++;
            return dpMemo.get(clusterHash);
        }
        
        // Get cluster size for base case check
        int taxaCount = clusterHashManager.getClusterSize(clusterHash);
        if (taxaCount < 0) {
            // Size not available - this shouldn't happen with proper preprocessing
            System.err.println("Warning: Cluster size not available for hash " + clusterHash.toDebugString());
            dpMemo.put(clusterHash, 0.0);
            return 0.0;
        }
        
        // Base case: clusters with ≤ 2 taxa
        if (taxaCount <= 2) {
            dpMemo.put(clusterHash, 0.0);
            return 0.0;
        }
        
        double maxScore = Double.NEGATIVE_INFINITY;
        Object bestChoice = null;  // Can be RangeBipartition or MixedBipartition
        
        // ========================================
        // Check RangeBipartition candidates (existing logic)
        // ========================================
        List<RangeBipartition> rangeCandidates = clusterHashToRangeBips.get(clusterHash);
        if (rangeCandidates != null) {
            for (RangeBipartition rangeBip : rangeCandidates) {
                if (isValidPartition(rangeBip, clusterHash)) {
                    clusterValidations++;
                    
                    // Get left and right cluster hashes
                    ClusterHashPair leftHash = clusterHashManager.getLeftClusterHash(rangeBip);
                    ClusterHashPair rightHash = clusterHashManager.getRightClusterHash(rangeBip);
                    
                    // Recursive DP calls
                    double leftScore = dp(leftHash);
                    double rightScore = dp(rightHash);
                    double rangeBipScore = rangeBipWeights.getOrDefault(rangeBip, 0.0);
                    
                    double totalScore = leftScore + rightScore + rangeBipScore;
                    
                    if (totalScore > maxScore) {
                        maxScore = totalScore;
                        bestChoice = rangeBip;
                    }
                }
            }
        }
        
        // ========================================
        // Check MixedBipartition candidates (if enabled)
        // ========================================
        if (useMixedBipartitions && clusterHashToMixedBips != null) {
            // System.out.println("***** Checking mixed bipartition candidates for cluster " + clusterHash.toDebugString());
            List<MixedBipartition> mixedCandidates = clusterHashToMixedBips.get(clusterHash);
            if (mixedCandidates != null) {
                // System.out.println("***** Found " + mixedCandidates.size() + " mixed bipartition candidates");
                for (MixedBipartition mixed : mixedCandidates) {
                    // System.out.println("***** Checking mixed bipartition " + mixed.toDebugString());
                    clusterValidations++;
                    
                    // Get left and right cluster hashes (from potentially different trees)
                    ClusterHashPair leftHash = clusterHashManager.getClusterHash(
                        mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd);
                    ClusterHashPair rightHash = clusterHashManager.getClusterHash(
                        mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd);
                    
                    // Recursive DP calls
                    double leftScore = dp(leftHash);
                    double rightScore = dp(rightHash);
                    double mixedScore = mixedBipWeights.getOrDefault(mixed, 0.0);
                    
                    double totalScore = leftScore + rightScore + mixedScore;
                    
                    if (totalScore > maxScore) {
                        maxScore = totalScore;
                        bestChoice = mixed;
                        mixedBipChoices++;
                    }
                }
            }
        }
        
        // Handle case where no valid bipartition found
        if (bestChoice == null) {
            if (rangeCandidates == null && (!useMixedBipartitions || clusterHashToMixedBips.get(clusterHash) == null)) {
                System.out.println("No candidates found for cluster " + clusterHash.toDebugString() + " (size " + taxaCount + ")");
            }
            maxScore = 0.0;
            System.out.println("No valid bipartition found for cluster " + clusterHash.toDebugString());
        }
        
        // Memoize result
        dpMemo.put(clusterHash, maxScore);
        dpChoice.put(clusterHash, bestChoice);
        
        return maxScore;
    }
    
    /**
     * Validate that a range bipartition is valid for the given cluster.
     * Uses hash-based comparison for efficiency.
     */
    private boolean isValidPartition(RangeBipartition rangeBip, ClusterHashPair clusterHash) {
        // Verify that the union of left and right ranges equals the cluster represented by clusterHash
        ClusterHashPair unionHash = clusterHashManager.getUnionClusterHash(rangeBip);
        
        // Use cluster hash manager's equality check (with collision handling)
        return clusterHashManager.clustersEqual(unionHash, clusterHash);
    }
    
    /**
     * Reconstruct the optimal tree using hash-based DP results.
     * This converts back to the traditional Tree structure for output.
     */
    public Tree reconstructTree() {
        System.out.println("==== RECONSTRUCTING OPTIMAL TREE ====");
        
        // Create hash for all taxa cluster
        ClusterHashPair allTaxaHash = clusterHashManager.getAllTaxaClusterHash();
        
        Tree tree = new Tree();
        tree.taxaMap = geneTrees.taxaMap;
        tree.root = buildTreeNode(allTaxaHash, tree);
        tree.isRooted = true;
        
        // Initialize the leaves array and count
        tree.leaves = new TreeNode[geneTrees.realTaxaCount];
        tree.leavesCount = 0;
        
        // Populate leaves array and count leaves
        for (TreeNode node : tree.nodes) {
            if (node.isLeaf()) {
                tree.leaves[node.taxon.id] = node;
                tree.leavesCount++;
            }
        }
        
        System.out.println("Tree reconstruction completed");
        System.out.println("Tree nodes: " + tree.nodes.size());
        System.out.println("Leaf nodes: " + tree.leavesCount);
        
        return tree;
    }
    
    /**
     * Build tree node recursively using hash-based DP choices.
     * Handles both RangeBipartition and MixedBipartition choices.
     */
    private TreeNode buildTreeNode(ClusterHashPair clusterHash, Tree tree) {
        int taxaCount = clusterHashManager.getClusterSize(clusterHash);
        
        if (taxaCount == 1) {
            // Single taxon - create leaf node
            Set<Integer> taxonSet = clusterHashManager.getTaxonSet(clusterHash);
            if (taxonSet != null && taxonSet.size() == 1) {
                int taxonId = taxonSet.iterator().next();
                return tree.addLeaf(geneTrees.taxa[taxonId]);
            } else {
                throw new RuntimeException("Cannot determine taxon ID for singleton cluster: " + 
                                         clusterHash.toDebugString());
            }
        }
        
        Object choice = dpChoice.get(clusterHash);
        if (choice == null) {
            if (taxaCount == 2) {
                // Two taxa - create binary internal node with two leaves
                Set<Integer> taxonSet = clusterHashManager.getTaxonSet(clusterHash);
                if (taxonSet != null && taxonSet.size() == 2) {
                    ArrayList<TreeNode> children = new ArrayList<>();
                    for (int taxonId : taxonSet) {
                        children.add(tree.addLeaf(geneTrees.taxa[taxonId]));
                    }
                    return tree.addInternalNode(children);
                } else {
                    throw new RuntimeException("Cannot determine taxa for binary cluster: " + 
                                             clusterHash.toDebugString());
                }
            }
            throw new RuntimeException("No valid bipartition found for cluster: " + 
                                     clusterHash.toDebugString() + ", taxaCount: " + taxaCount);
        }
        
        // Get child cluster hashes based on choice type
        ClusterHashPair leftHash;
        ClusterHashPair rightHash;
        
        if (choice instanceof RangeBipartition) {
            RangeBipartition rangeBip = (RangeBipartition) choice;
            leftHash = clusterHashManager.getLeftClusterHash(rangeBip);
            rightHash = clusterHashManager.getRightClusterHash(rangeBip);
        } else if (choice instanceof MixedBipartition) {
            System.out.println("Mixed bipartition choice found");
            MixedBipartition mixed = (MixedBipartition) choice;
            leftHash = clusterHashManager.getClusterHash(
                mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd);
            rightHash = clusterHashManager.getClusterHash(
                mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd);
        } else {
            throw new RuntimeException("Unknown choice type: " + choice.getClass().getName());
        }
        
        // Recursive construction
        TreeNode leftChild = buildTreeNode(leftHash, tree);
        TreeNode rightChild = buildTreeNode(rightHash, tree);
        
        ArrayList<TreeNode> children = new ArrayList<>();
        children.add(leftChild);
        children.add(rightChild);
        
        return tree.addInternalNode(children);
    }
    
    /**
     * Print DP table for debugging (limited output for hash-based keys).
     */
    public void printDPTable() {
        System.out.println("Memory-Optimized DP Results:");
        System.out.println("Total memoized clusters: " + dpMemo.size());
        
        // Print a sample of results (not all, as there could be many)
        int printed = 0;
        for (Map.Entry<ClusterHashPair, Double> entry : dpMemo.entrySet()) {
            ClusterHashPair clusterHash = entry.getKey();
            double score = entry.getValue();
            int size = clusterHashManager.getClusterSize(clusterHash);
            
            System.out.println("Cluster " + clusterHash.toDebugString() + 
                             " (size " + size + "): score = " + score);
            
            printed++;
            if (printed >= 10) {
                System.out.println("... (showing first 10 of " + dpMemo.size() + " total)");
                break;
            }
        }
    }
    
    /**
     * Get processing statistics for performance analysis.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Memory-Optimized Inference DP Statistics:\n");
        sb.append("  Total processing time: ").append(totalProcessingTime).append(" ms\n");
        sb.append("  DP calls: ").append(dpCalls).append("\n");
        sb.append("  Memo hits: ").append(memoHits).append(" (").append(
            dpCalls > 0 ? String.format("%.2f%%", 100.0 * memoHits / dpCalls) : "0%").append(")\n");
        sb.append("  Cluster validations: ").append(clusterValidations).append("\n");
        sb.append("  Memoized clusters: ").append(dpMemo.size()).append("\n");
        sb.append("  Cluster hash groups (range): ").append(clusterHashToRangeBips.size()).append("\n");
        
        if (useMixedBipartitions) {
            sb.append("  Mixed bipartitions enabled: YES\n");
            sb.append("  Cluster hash groups (mixed): ").append(
                clusterHashToMixedBips != null ? clusterHashToMixedBips.size() : 0).append("\n");
            sb.append("  Mixed bipartition choices: ").append(mixedBipChoices).append("\n");
        }
        
        if (dpCalls > 0) {
            sb.append("  Average time per DP call: ").append(totalProcessingTime / (double) dpCalls).append(" ms\n");
        }
        
        return sb.toString();
    }
}
