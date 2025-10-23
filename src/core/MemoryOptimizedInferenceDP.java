package core;

import java.util.*;
import preprocessing.GeneTrees;
import tree.*;
import utils.BitSet;

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
    private final List<STBipartition> candidateSTBips;
    private final Map<ClusterHashPair, List<STBipartition>> clusterHashToSTBips;
    private Map<STBipartition, Double> stbipWeights;
    private final Map<ClusterHashPair, Double> dpMemo;
    private final Map<ClusterHashPair, STBipartition> dpChoice;
    private final ClusterHashManager clusterHashManager;
    private final MemoryEfficientBipartitionManager bipartitionManager;
    
    // Statistics for monitoring
    private long dpCalls = 0;
    private long memoHits = 0;
    private long clusterValidations = 0;
    private long totalProcessingTime = 0;
    
    public MemoryOptimizedInferenceDP(GeneTrees geneTrees, List<STBipartition> candidateSTBips) {
        System.out.println("==== INITIALIZING MEMORY-OPTIMIZED INFERENCE DP ====");
        System.out.println("Number of candidate bipartitions: " + candidateSTBips.size());
        System.out.println("Real taxa count: " + geneTrees.realTaxaCount);
        
        this.geneTrees = geneTrees;
        this.candidateSTBips = candidateSTBips;
        this.clusterHashToSTBips = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
        
        // Initialize bipartition manager and cluster hash manager
        System.out.println("Initializing MemoryEfficientBipartitionManager...");
        this.bipartitionManager = new MemoryEfficientBipartitionManager(
            geneTrees.geneTrees, geneTrees.realTaxaCount);
        bipartitionManager.processGeneTreesParallel();
        
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
        System.out.println("Cluster hash groups: " + clusterHashToSTBips.size());
        System.out.println("Cached cluster hashes: " + clusterHashManager.getStatistics().split("\\n")[4]);
    }
    
    /**
     * Preprocess candidates to create hash-based cluster mappings.
     * This replaces the BitSet-based preprocessing in the original implementation.
     */
    private void preprocessCandidatesWithHashes() {
        System.out.println("Creating hash-based cluster mappings for candidates...");
        
        int processedCandidates = 0;
        int uniqueClusters = 0;
        
        for (STBipartition stbip : candidateSTBips) {
            // Compute union cluster hash (left ∪ right)
            BitSet union = (BitSet) stbip.cluster1.clone();
            union.or(stbip.cluster2);
            
            ClusterHashPair unionHash = clusterHashManager.getClusterHash(union);
            
            // Add to hash-based mapping
            List<STBipartition> candidates = clusterHashToSTBips.get(unionHash);
            if (candidates == null) {
                candidates = new ArrayList<>();
                clusterHashToSTBips.put(unionHash, candidates);
                uniqueClusters++;
            }
            candidates.add(stbip);
            
            processedCandidates++;
            
            // Log progress for large datasets
            if (processedCandidates % 1000 == 0 || processedCandidates == candidateSTBips.size()) {
                System.out.println("Processed " + processedCandidates + "/" + candidateSTBips.size() + 
                                 " candidates, unique clusters: " + uniqueClusters);
            }
        }
        
        System.out.println("Hash-based preprocessing completed");
        System.out.println("Created mappings for " + uniqueClusters + " unique cluster hashes");
        System.out.println("Average candidates per cluster: " + 
                         (uniqueClusters > 0 ? candidateSTBips.size() / (double) uniqueClusters : 0));
    }
    
    /**
     * Calculate weights using memory-optimized weight calculator.
     */
    private void calculateWeights() {
        MemoryOptimizedWeightCalculator calculator = new MemoryOptimizedWeightCalculator(geneTrees);
        stbipWeights = calculator.calculateWeights(candidateSTBips);
        
        System.out.println("Weight calculation completed for " + stbipWeights.size() + " bipartitions");
        
        // Print weight statistics
        double minWeight = stbipWeights.values().stream().mapToDouble(Double::doubleValue).min().orElse(0.0);
        double maxWeight = stbipWeights.values().stream().mapToDouble(Double::doubleValue).max().orElse(0.0);
        double avgWeight = stbipWeights.values().stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        
        System.out.println("Weight statistics: min=" + minWeight + ", max=" + maxWeight + ", avg=" + avgWeight);
    }
    
    /**
     * Solve the DP problem using hash-based cluster keys.
     */
    public double solve() {
        System.out.println("==== STARTING MEMORY-OPTIMIZED DP SOLUTION ====");
        
        long startTime = System.currentTimeMillis();
        
        // Create hash for all taxa cluster
        BitSet allTaxa = new BitSet(geneTrees.realTaxaCount);
        for (int i = 0; i < geneTrees.realTaxaCount; i++) {
            allTaxa.set(i);
        }
        
        ClusterHashPair allTaxaHash = clusterHashManager.getClusterHash(allTaxa);
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
        
        // Print cluster hash manager statistics
        System.out.println("\n" + clusterHashManager.getStatistics());
        
        return result;
    }
    
    /**
     * Main DP function using cluster hashes instead of BitSets.
     * This is the core optimization - replacing O(n) BitSet keys with O(1) hash keys.
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
        STBipartition bestChoice = null;
        
        // Get candidate bipartitions for this cluster
        List<STBipartition> candidates = clusterHashToSTBips.get(clusterHash);
        if (candidates != null) {
            System.out.println("Processing " + candidates.size() + " candidates for cluster " + 
                             clusterHash.toDebugString() + " (size " + taxaCount + ")");
            
            for (STBipartition stbip : candidates) {
                if (isValidPartition(stbip, clusterHash)) {
                    clusterValidations++;
                    
                    // Get left and right cluster hashes
                    ClusterHashPair leftHash = clusterHashManager.getClusterHash(stbip.cluster1);
                    ClusterHashPair rightHash = clusterHashManager.getClusterHash(stbip.cluster2);
                    
                    // Recursive DP calls
                    double leftScore = dp(leftHash);
                    double rightScore = dp(rightHash);
                    double stbipScore = stbipWeights.getOrDefault(stbip, 0.0);
                    
                    double totalScore = leftScore + rightScore + stbipScore;
                    
                    if (totalScore > maxScore) {
                        maxScore = totalScore;
                        bestChoice = stbip;
                    }
                }
            }
        } else {
            System.out.println("No candidates found for cluster " + clusterHash.toDebugString() + " (size " + taxaCount + ")");
        }
        
        if (bestChoice == null) {
            maxScore = 0.0;
            System.out.println("No valid bipartition found for cluster " + clusterHash.toDebugString());
        } else {
            System.out.println("Best choice for cluster " + clusterHash.toDebugString() + 
                             ": score=" + maxScore + ", bipartition=" + bestChoice);
        }
        
        // Memoize result
        dpMemo.put(clusterHash, maxScore);
        dpChoice.put(clusterHash, bestChoice);
        
        return maxScore;
    }
    
    /**
     * Validate that a bipartition is valid for the given cluster.
     * Uses hash-based comparison with fallback to expensive verification.
     */
    private boolean isValidPartition(STBipartition stbip, ClusterHashPair clusterHash) {
        // Verify that stbip.cluster1 ∪ stbip.cluster2 == cluster represented by clusterHash
        BitSet union = (BitSet) stbip.cluster1.clone();
        union.or(stbip.cluster2);
        
        ClusterHashPair unionHash = clusterHashManager.getClusterHash(union);
        
        // Use cluster hash manager's equality check (with collision handling)
        return clusterHashManager.clustersEqual(unionHash, clusterHash);
    }
    
    /**
     * Reconstruct the optimal tree using hash-based DP results.
     * This converts back to the traditional Tree structure for output.
     */
    public Tree reconstructTree() {
        System.out.println("==== RECONSTRUCTING OPTIMAL TREE ====");
        
        // Create hash for all taxa
        BitSet allTaxa = new BitSet(geneTrees.realTaxaCount);
        for (int i = 0; i < geneTrees.realTaxaCount; i++) {
            allTaxa.set(i);
        }
        ClusterHashPair allTaxaHash = clusterHashManager.getClusterHash(allTaxa);
        
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
        
        STBipartition choice = dpChoice.get(clusterHash);
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
        
        // Recursive construction
        ClusterHashPair leftHash = clusterHashManager.getClusterHash(choice.cluster1);
        ClusterHashPair rightHash = clusterHashManager.getClusterHash(choice.cluster2);
        
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
        sb.append("  Cluster hash groups: ").append(clusterHashToSTBips.size()).append("\n");
        
        if (dpCalls > 0) {
            sb.append("  Average time per DP call: ").append(totalProcessingTime / (double) dpCalls).append(" ms\n");
        }
        
        return sb.toString();
    }
}
