package core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.Arrays;

import preprocessing.MemoryOptimizedGeneTrees;
import tree.RangeSTBipartition;
import tree.Tree;
import tree.TreeNode;

/**
 * Memory-optimized inference algorithm using range-based bipartitions.
 * Adapts the dynamic programming approach to work with hash-based cluster representation
 * while maintaining the same algorithmic complexity.
 */
public class MemoryOptimizedInferenceDP {
    
    private MemoryOptimizedGeneTrees geneTrees;
    private List<RangeSTBipartition> candidateRangeBips;
    private Map<String, List<RangeSTBipartition>> clusterHashToRangeBips;  // Hash-based clustering
    private Map<RangeSTBipartition, Double> rangeBipWeights;
    private Map<String, Double> dpMemo;                                    // Hash-based DP memoization
    private Map<String, RangeSTBipartition> dpChoice;                      // Hash-based DP choices
    
    // For cluster representation using sorted taxon arrays
    private Map<String, int[]> hashToTaxaArray;                           // Hash to taxa array mapping
    
    public MemoryOptimizedInferenceDP(MemoryOptimizedGeneTrees geneTrees, 
                                     List<RangeSTBipartition> candidateRangeBips) {
        this.geneTrees = geneTrees;
        this.candidateRangeBips = candidateRangeBips;
        this.clusterHashToRangeBips = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
        this.hashToTaxaArray = new HashMap<>();
        
        preprocessCandidates();
        calculateWeights();
    }
    
    /**
     * Groups candidate bipartitions by their union cluster for efficient DP lookup.
     */
    private void preprocessCandidates() {
        System.out.println("Preprocessing " + candidateRangeBips.size() + " range-based candidates...");
        
        int[][] orderings = geneTrees.getGeneTreeOrderings();
        
        for (RangeSTBipartition rangeBip : candidateRangeBips) {
            // Get union of cluster1 and cluster2 (all taxa in the bipartition)
            int[] cluster1Taxa = rangeBip.getTaxaIds(orderings);
            int[] cluster2Taxa = rangeBip.getComplementTaxaIds(orderings);
            
            // Combine to get the full cluster
            Set<Integer> unionSet = new HashSet<>();
            for (int taxon : cluster1Taxa) unionSet.add(taxon);
            for (int taxon : cluster2Taxa) unionSet.add(taxon);
            
            // Convert to sorted array for consistent hashing
            int[] unionArray = unionSet.stream().mapToInt(Integer::intValue).sorted().toArray();
            String clusterHash = Arrays.toString(unionArray);
            
            // Store the mapping
            hashToTaxaArray.put(clusterHash, unionArray);
            clusterHashToRangeBips.computeIfAbsent(clusterHash, k -> new ArrayList<>()).add(rangeBip);
        }
        
        System.out.println("Created " + clusterHashToRangeBips.size() + " unique cluster groups");
    }
    
    /**
     * Calculate weights for all candidate range bipartitions.
     */
    private void calculateWeights() {
        System.out.println("Calculating weights for range-based bipartitions...");
        MemoryOptimizedWeightCalculatorCPUOnly calculator = new MemoryOptimizedWeightCalculatorCPUOnly(geneTrees);
        rangeBipWeights = calculator.calculateWeights(candidateRangeBips);
        System.out.println("Weight calculation completed");
    }
    
    /**
     * Solve the inference problem using dynamic programming.
     */
    public double solve() {
        // Create cluster for all taxa
        int[] allTaxa = new int[geneTrees.realTaxaCount];
        for (int i = 0; i < geneTrees.realTaxaCount; i++) {
            allTaxa[i] = i;
        }
        String allTaxaHash = Arrays.toString(allTaxa);
        hashToTaxaArray.put(allTaxaHash, allTaxa);
        
        return dp(allTaxaHash);
    }
    
    /**
     * Dynamic programming function using hash-based cluster representation.
     */
    private double dp(String clusterHash) {
        if (dpMemo.containsKey(clusterHash)) {
            return dpMemo.get(clusterHash);
        }
        
        int[] clusterTaxa = hashToTaxaArray.get(clusterHash);
        if (clusterTaxa == null) {
            System.err.println("Warning: Unknown cluster hash: " + clusterHash);
            return 0.0;
        }
        
        int taxaCount = clusterTaxa.length;
        if (taxaCount <= 2) {
            dpMemo.put(clusterHash, 0.0);
            return 0.0;
        }
        
        double maxScore = Double.NEGATIVE_INFINITY;
        RangeSTBipartition bestChoice = null;
        
        List<RangeSTBipartition> candidates = clusterHashToRangeBips.get(clusterHash);
        if (candidates != null) {
            int[][] orderings = geneTrees.getGeneTreeOrderings();
            
            for (RangeSTBipartition rangeBip : candidates) {
                if (isValidPartition(rangeBip, clusterTaxa, orderings)) {
                    // Get subclusters
                    int[] leftCluster = rangeBip.getTaxaIds(orderings);
                    int[] rightCluster = rangeBip.getComplementTaxaIds(orderings);
                    
                    // Create hashes for subclusters
                    Arrays.sort(leftCluster);
                    Arrays.sort(rightCluster);
                    String leftHash = Arrays.toString(leftCluster);
                    String rightHash = Arrays.toString(rightCluster);
                    
                    // Store cluster mappings
                    hashToTaxaArray.put(leftHash, leftCluster);
                    hashToTaxaArray.put(rightHash, rightCluster);
                    
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
        
        if (bestChoice == null) {
            maxScore = 0.0;
        }
        
        dpMemo.put(clusterHash, maxScore);
        dpChoice.put(clusterHash, bestChoice);
        return maxScore;
    }
    
    /**
     * Check if a range bipartition is valid for the given cluster.
     */
    private boolean isValidPartition(RangeSTBipartition rangeBip, int[] clusterTaxa, int[][] orderings) {
        int[] cluster1Taxa = rangeBip.getTaxaIds(orderings);
        int[] cluster2Taxa = rangeBip.getComplementTaxaIds(orderings);
        
        // Check if union of cluster1 and cluster2 equals the original cluster
        Set<Integer> unionSet = new HashSet<>();
        for (int taxon : cluster1Taxa) unionSet.add(taxon);
        for (int taxon : cluster2Taxa) unionSet.add(taxon);
        
        if (unionSet.size() != clusterTaxa.length) {
            return false;
        }
        
        for (int taxon : clusterTaxa) {
            if (!unionSet.contains(taxon)) {
                return false;
            }
        }
        
        return true;
    }
    
    /**
     * Reconstruct the final tree from DP choices.
     */
    public Tree reconstructTree() {
        int[] allTaxa = new int[geneTrees.realTaxaCount];
        for (int i = 0; i < geneTrees.realTaxaCount; i++) {
            allTaxa[i] = i;
        }
        String allTaxaHash = Arrays.toString(allTaxa);
        
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
        
        return tree;
    }
    
    /**
     * Recursively build tree nodes from DP choices.
     */
    private TreeNode buildTreeNode(String clusterHash, Tree tree) {
        int[] clusterTaxa = hashToTaxaArray.get(clusterHash);
        if (clusterTaxa == null) {
            throw new RuntimeException("Unknown cluster hash: " + clusterHash);
        }
        
        int taxaCount = clusterTaxa.length;
        
        if (taxaCount == 1) {
            int taxonId = clusterTaxa[0];
            return tree.addLeaf(geneTrees.taxa[taxonId]);
        }
        
        RangeSTBipartition choice = dpChoice.get(clusterHash);
        if (choice == null) {
            if (taxaCount == 2) {
                ArrayList<TreeNode> children = new ArrayList<>();
                for (int taxonId : clusterTaxa) {
                    children.add(tree.addLeaf(geneTrees.taxa[taxonId]));
                }
                return tree.addInternalNode(children);
            }
            throw new RuntimeException("No valid bipartition found for cluster with " + taxaCount + " taxa");
        }
        
        // Build left and right subtrees
        int[][] orderings = geneTrees.getGeneTreeOrderings();
        int[] leftCluster = choice.getTaxaIds(orderings);
        int[] rightCluster = choice.getComplementTaxaIds(orderings);
        
        Arrays.sort(leftCluster);
        Arrays.sort(rightCluster);
        String leftHash = Arrays.toString(leftCluster);
        String rightHash = Arrays.toString(rightCluster);
        
        TreeNode leftChild = buildTreeNode(leftHash, tree);
        TreeNode rightChild = buildTreeNode(rightHash, tree);
        
        ArrayList<TreeNode> children = new ArrayList<>();
        children.add(leftChild);
        children.add(rightChild);
        return tree.addInternalNode(children);
    }
    
    /**
     * Print DP table for debugging.
     */
    public void printDPTable() {
        System.out.println("Memory-Optimized DP Results:");
        System.out.println("Total clusters processed: " + dpMemo.size());
        
        // Print sample results
        int count = 0;
        for (Map.Entry<String, Double> entry : dpMemo.entrySet()) {
            if (count < 10) {  // Print first 10 for debugging
                String clusterHash = entry.getKey();
                double score = entry.getValue();
                int[] taxa = hashToTaxaArray.get(clusterHash);
                System.out.println("Cluster " + Arrays.toString(taxa) + ": Score = " + score);
                count++;
            } else {
                break;
            }
        }
        
        if (dpMemo.size() > 10) {
            System.out.println("... and " + (dpMemo.size() - 10) + " more clusters");
        }
    }
    
    /**
     * Get memory usage statistics.
     */
    public void printMemoryStats() {
        long clusterMappings = hashToTaxaArray.size();
        long candidateGroups = clusterHashToRangeBips.size();
        long dpEntries = dpMemo.size();
        
        System.out.println("\n=== Memory-Optimized Inference Statistics ===");
        System.out.println("Cluster mappings: " + clusterMappings);
        System.out.println("Candidate groups: " + candidateGroups);
        System.out.println("DP entries: " + dpEntries);
        System.out.println("Total range bipartitions: " + candidateRangeBips.size());
        System.out.println("==============================================\n");
    }
}
