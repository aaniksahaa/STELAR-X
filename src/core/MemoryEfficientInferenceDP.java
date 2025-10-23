package core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;

import preprocessing.GeneTrees;
import tree.CompactBipartition;
import tree.ClusterHash;
import tree.RangeBipartition;
import tree.MemoryEfficientBipartitionManager;
import tree.STBipartition;
import tree.Tree;
import tree.TreeNode;

/**
 * Memory-efficient DP inference that uses hash-based cluster representation
 * instead of BitSet keys, avoiding memory-intensive bitset expansions.
 * 
 * Key optimizations:
 * 1. Uses ClusterHash as keys instead of BitSet for DP memoization
 * 2. Maintains mapping from ClusterHash to CompactBipartition candidates
 * 3. Uses MemoryEfficientWeightCalculator for weight computation
 * 4. Only expands to BitSets during final tree reconstruction
 */
public class MemoryEfficientInferenceDP {
    
    private GeneTrees geneTrees;
    private List<CompactBipartition> candidateCompactBips;
    private Map<ClusterHash, List<CompactBipartition>> clusterHashToCompactBips;
    private Map<CompactBipartition, Double> compactBipWeights;
    private Map<ClusterHash, Double> dpMemo;
    private Map<ClusterHash, CompactBipartition> dpChoice;
    
    // Data structures for hash calculation and reconstruction
    private RangeBipartition.HashFunction hashFunction;
    private long[][] prefixSums;
    private long[][] prefixXORs;
    private int[][] geneTreeOrderings;
    private MemoryEfficientWeightCalculator weightCalculator;
    
    public MemoryEfficientInferenceDP(GeneTrees geneTrees, List<STBipartition> candidateSTBips) {
        System.out.println("\n==== INITIALIZING MEMORY-EFFICIENT DP ====");
        System.out.println("Constructor called with:");
        System.out.println("  GeneTrees: " + (geneTrees != null ? "Valid" : "NULL"));
        System.out.println("  Candidate STBipartitions: " + (candidateSTBips != null ? candidateSTBips.size() : "NULL"));
        
        if (geneTrees != null) {
            System.out.println("  Gene trees count: " + geneTrees.geneTrees.size());
            System.out.println("  Real taxa count: " + geneTrees.realTaxaCount);
            System.out.println("  STBipartitions: " + geneTrees.stBipartitions.size());
        }
        
        this.geneTrees = geneTrees;
        this.clusterHashToCompactBips = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
        this.hashFunction = MemoryEfficientBipartitionManager.DEFAULT_HASH_FUNCTION;
        
        System.out.println("\nStarting initialization phases...");
        
        try {
            System.out.println("Phase 1: Initializing data structures...");
            initializeDataStructures();
            System.out.println("Phase 1 completed successfully");
            
            System.out.println("Phase 2: Converting candidates to compact...");
            convertCandidatesToCompact(candidateSTBips);
            System.out.println("Phase 2 completed successfully");
            
            System.out.println("Phase 3: Preprocessing candidates...");
            preprocessCandidates();
            System.out.println("Phase 3 completed successfully");
            
            System.out.println("Phase 4: Calculating weights...");
            calculateWeights(candidateSTBips);
            System.out.println("Phase 4 completed successfully");
            
            System.out.println("\nMemory-efficient DP initialization completed successfully!");
            
        } catch (Exception e) {
            System.err.println("ERROR during MemoryEfficientInferenceDP initialization:");
            System.err.println("  Error message: " + e.getMessage());
            System.err.println("  Error type: " + e.getClass().getSimpleName());
            e.printStackTrace();
            throw e;
        }
    }
    
    /**
     * Initialize data structures needed for hash-based DP.
     */
    private void initializeDataStructures() {
        System.out.println("\n=== Initializing Memory-Efficient DP ===");
        
        int numTrees = geneTrees.geneTrees.size();
        int numTaxa = geneTrees.realTaxaCount;
        
        System.out.println("Setting up DP data structures...");
        System.out.println("  Gene trees: " + numTrees);
        System.out.println("  Taxa count: " + numTaxa);
        System.out.println("  Hash function: " + hashFunction.getName());
        
        // Initialize arrays for hash calculation
        System.out.println("\nAllocating DP-specific arrays...");
        geneTreeOrderings = new int[numTrees][];
        prefixSums = new long[numTrees][];
        prefixXORs = new long[numTrees][];
        
        // Build gene tree orderings and prefix arrays
        System.out.println("Building gene tree orderings for DP...");
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            List<Integer> ordering = collectLeavesInOrder(geneTrees.geneTrees.get(treeIdx).root);
            geneTreeOrderings[treeIdx] = ordering.stream().mapToInt(Integer::intValue).toArray();
            
            if (treeIdx < 3) {
                System.out.println("  Tree " + treeIdx + ": " + ordering.size() + " taxa ordered");
            }
        }
        
        System.out.println("Calculating prefix arrays for DP hash functions...");
        calculatePrefixArrays();
        
        // Initialize weight calculator
        System.out.println("Initializing memory-efficient weight calculator...");
        weightCalculator = new MemoryEfficientWeightCalculator(geneTrees);
        
        System.out.println("\nMemory-efficient DP initialization completed:");
        System.out.println("  Gene tree orderings: " + numTrees + " arrays");
        System.out.println("  Prefix arrays: " + (numTrees * 2) + " arrays");
        System.out.println("  Hash-based memoization: Ready");
        System.out.println("  Weight calculator: Initialized");
    }
    
    private List<Integer> collectLeavesInOrder(TreeNode node) {
        List<Integer> ordering = new ArrayList<>();
        collectLeavesInOrderRecursive(node, ordering);
        return ordering;
    }
    
    private void collectLeavesInOrderRecursive(TreeNode node, List<Integer> ordering) {
        if (node.isLeaf()) {
            ordering.add(node.taxon.id);
            return;
        }
        
        if (node.childs != null && node.childs.size() >= 2) {
            collectLeavesInOrderRecursive(node.childs.get(0), ordering);
            collectLeavesInOrderRecursive(node.childs.get(1), ordering);
        }
        
        for (int i = 2; i < (node.childs != null ? node.childs.size() : 0); i++) {
            collectLeavesInOrderRecursive(node.childs.get(i), ordering);
        }
    }
    
    private void calculatePrefixArrays() {
        for (int i = 0; i < geneTrees.geneTrees.size(); i++) {
            int[] ordering = geneTreeOrderings[i];
            prefixSums[i] = new long[ordering.length];
            prefixXORs[i] = new long[ordering.length];
            
            if (ordering.length > 0) {
                long hashedTaxon0 = hashSingleTaxon(ordering[0]);
                prefixSums[i][0] = hashedTaxon0;
                prefixXORs[i][0] = hashedTaxon0;
                
                for (int j = 1; j < ordering.length; j++) {
                    long hashedTaxon = hashSingleTaxon(ordering[j]);
                    prefixSums[i][j] = prefixSums[i][j - 1] + hashedTaxon;
                    prefixXORs[i][j] = prefixXORs[i][j - 1] ^ hashedTaxon;
                }
            }
        }
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
     * Convert STBipartition candidates to CompactBipartition representation.
     */
    private void convertCandidatesToCompact(List<STBipartition> candidateSTBips) {
        System.out.println("Converting candidate STBipartitions to compact representation...");
        
        candidateCompactBips = new ArrayList<>();
        
        // For now, create placeholder compact bipartitions
        // In a full implementation, we would extract these from the original gene tree processing
        for (STBipartition stb : candidateSTBips) {
            // Create a placeholder compact bipartition
            // This would be properly extracted from range bipartitions in a complete implementation
            CompactBipartition compactBip = new CompactBipartition(0, 0, stb.cluster1.cardinality(),
                                                                  stb.cluster1.cardinality(),
                                                                  stb.cluster1.cardinality() + stb.cluster2.cardinality());
            candidateCompactBips.add(compactBip);
        }
        
        System.out.println("Converted " + candidateCompactBips.size() + " candidates to compact representation");
    }
    
    /**
     * Preprocess candidates by grouping them by cluster hash.
     */
    private void preprocessCandidates() {
        System.out.println("\n=== Preprocessing Candidates for DP ===");
        System.out.println("Grouping candidates by cluster hash...");
        System.out.println("  Compact bipartitions to process: " + candidateCompactBips.size());
        
        int processedCount = 0;
        Map<Integer, Integer> sizeDistribution = new HashMap<>();
        
        for (CompactBipartition compactBip : candidateCompactBips) {
            // Calculate cluster hash for the full range (union of left and right)
            // This bipartition can split this full cluster
            ClusterHash clusterHash = ClusterHash.fromCompactBipartition(compactBip, hashFunction, 
                                                                        prefixSums, prefixXORs, geneTreeOrderings);
            
            clusterHashToCompactBips.computeIfAbsent(clusterHash, k -> new ArrayList<>()).add(compactBip);
            
            // Track size distribution
            int clusterSize = clusterHash.getSize();
            sizeDistribution.put(clusterSize, sizeDistribution.getOrDefault(clusterSize, 0) + 1);
            
            processedCount++;
            // Removed verbose per-item logging
        }
        
        System.out.println("\nCandidate preprocessing completed:");
        System.out.println("  Candidates processed: " + processedCount);
        System.out.println("  Unique cluster hash groups: " + clusterHashToCompactBips.size());
        System.out.println("  Average candidates per group: " + 
                          String.format("%.2f", (double) processedCount / clusterHashToCompactBips.size()));
        
        // Log size distribution
        System.out.println("  Cluster size distribution:");
        sizeDistribution.entrySet().stream()
            .sorted(Map.Entry.comparingByKey())
            .forEach(entry -> System.out.println("    Size " + entry.getKey() + ": " + entry.getValue() + " clusters"));
    }
    
    /**
     * Calculate weights using memory-efficient weight calculator.
     */
    private void calculateWeights(List<STBipartition> candidateSTBips) {
        System.out.println("Calculating weights using memory-efficient approach...");
        
        Map<STBipartition, Double> stbWeights = weightCalculator.calculateWeights(candidateSTBips);
        
        // Convert STBipartition weights to CompactBipartition weights
        compactBipWeights = new HashMap<>();
        for (int i = 0; i < candidateSTBips.size() && i < candidateCompactBips.size(); i++) {
            STBipartition stb = candidateSTBips.get(i);
            CompactBipartition compactBip = candidateCompactBips.get(i);
            Double weight = stbWeights.get(stb);
            if (weight != null) {
                compactBipWeights.put(compactBip, weight);
            }
        }
        
        System.out.println("Calculated weights for " + compactBipWeights.size() + " compact bipartitions");
    }
    
    /**
     * Solve the DP problem using hash-based cluster representation.
     */
    public double solve() {
        System.out.println("\n==== STARTING MEMORY-EFFICIENT DP SOLVE ====");
        System.out.println("Total taxa count: " + geneTrees.realTaxaCount);
        System.out.println("Available cluster hash groups: " + clusterHashToCompactBips.size());
        System.out.println("Compact bipartition weights: " + compactBipWeights.size());
        
        // Create cluster hash for all taxa - must match candidate hash calculation method
        Set<Integer> allTaxaSet = new HashSet<>();
        for (int i = 0; i < geneTrees.realTaxaCount; i++) {
            allTaxaSet.add(i);
        }
        
        // Find a candidate that spans all taxa to get the correct hash
        ClusterHash allTaxaHash = null;
        for (Map.Entry<ClusterHash, List<CompactBipartition>> entry : clusterHashToCompactBips.entrySet()) {
            if (entry.getKey().getSize() == geneTrees.realTaxaCount) {
                allTaxaHash = entry.getKey();
                System.out.println("Found matching root cluster hash from candidates: " + allTaxaHash);
                break;
            }
        }
        
        if (allTaxaHash == null) {
            System.err.println("ERROR: No candidate spans all taxa! Creating fallback hash...");
            allTaxaHash = ClusterHash.fromTaxonSet(allTaxaSet, hashFunction);
        }
        
        System.out.println("Root cluster hash: " + allTaxaHash);
        System.out.println("Root cluster size: " + allTaxaHash.getSize());
        
        double result = dp(allTaxaHash);
        
        System.out.println("DP solve completed with optimal score: " + result);
        System.out.println("DP memoization table size: " + dpMemo.size());
        
        return result;
    }
    
    /**
     * DP function using cluster hash as key instead of BitSet.
     */
    private double dp(ClusterHash clusterHash) {
        if (dpMemo.containsKey(clusterHash)) {
            return dpMemo.get(clusterHash);
        }
        
        int taxaCount = clusterHash.getSize();
        System.out.println("  DP processing cluster: size=" + taxaCount + ", hash=" + clusterHash.getHash());
        
        if (taxaCount <= 2) {
            System.out.println("    Base case: cluster size <= 2, returning score 0.0");
            dpMemo.put(clusterHash, 0.0);
            return 0.0;
        }
        
        double maxScore = Double.NEGATIVE_INFINITY;
        CompactBipartition bestChoice = null;
        
        List<CompactBipartition> candidates = clusterHashToCompactBips.get(clusterHash);
        System.out.println("    Found " + (candidates != null ? candidates.size() : 0) + " candidate bipartitions");
        
        if (candidates != null) {
            int validCandidates = 0;
            for (CompactBipartition compactBip : candidates) {
                if (isValidPartition(compactBip, clusterHash)) {
                    validCandidates++;
                    
                    // Calculate cluster hashes for left and right subclusters
                    ClusterHash leftHash = ClusterHash.fromLeftCluster(compactBip, hashFunction, prefixSums, prefixXORs, geneTreeOrderings);
                    ClusterHash rightHash = ClusterHash.fromRightCluster(compactBip, hashFunction, prefixSums, prefixXORs, geneTreeOrderings);
                    
                    double leftScore = dp(leftHash);
                    double rightScore = dp(rightHash);
                    double compactBipScore = compactBipWeights.getOrDefault(compactBip, 0.0);
                    
                    double totalScore = leftScore + rightScore + compactBipScore;
                    
                    if (totalScore > maxScore) {
                        maxScore = totalScore;
                        bestChoice = compactBip;
                    }
                }
            }
            System.out.println("    Valid candidates: " + validCandidates + "/" + candidates.size());
        }
        
        if (bestChoice == null) {
            System.out.println("    WARNING: No valid bipartition found, using score 0.0");
            maxScore = 0.0;
        } else {
            System.out.println("    Best choice found with score: " + maxScore);
        }
        
        dpMemo.put(clusterHash, maxScore);
        dpChoice.put(clusterHash, bestChoice);
        return maxScore;
    }
    
    /**
     * Check if a compact bipartition is valid for the given cluster hash.
     */
    private boolean isValidPartition(CompactBipartition compactBip, ClusterHash clusterHash) {
        // Calculate the union hash for this bipartition
        ClusterHash unionHash = ClusterHash.fromCompactBipartition(compactBip, hashFunction, 
                                                                  prefixSums, prefixXORs, geneTreeOrderings);
        return unionHash.equals(clusterHash);
    }
    
    /**
     * Reconstruct tree using the DP solution.
     * This is where we need to convert back to BitSets for compatibility with existing Tree structure.
     */
    public Tree reconstructTree() {
        System.out.println("\n=== RECONSTRUCTING TREE FROM DP SOLUTION ===");
        System.out.println("Starting tree reconstruction...");
        System.out.println("  DP memoization entries: " + dpMemo.size());
        System.out.println("  DP choice entries: " + dpChoice.size());
        System.out.println("  Cluster hash groups: " + clusterHashToCompactBips.size());
        
        // Create cluster hash for all taxa - must match candidate hash calculation method
        Set<Integer> allTaxaSet = new HashSet<>();
        for (int i = 0; i < geneTrees.realTaxaCount; i++) {
            allTaxaSet.add(i);
        }
        
        // Find a candidate that spans all taxa to get the correct hash
        ClusterHash allTaxaHash = null;
        for (Map.Entry<ClusterHash, List<CompactBipartition>> entry : clusterHashToCompactBips.entrySet()) {
            if (entry.getKey().getSize() == geneTrees.realTaxaCount) {
                allTaxaHash = entry.getKey();
                System.out.println("Using matching root cluster hash from candidates: " + allTaxaHash);
                break;
            }
        }
        
        if (allTaxaHash == null) {
            System.err.println("ERROR: No candidate spans all taxa! Creating fallback hash...");
            allTaxaHash = ClusterHash.fromTaxonSet(allTaxaSet, hashFunction);
        }
        
        System.out.println("Root cluster created:");
        System.out.println("  All taxa count: " + allTaxaSet.size());
        System.out.println("  Root cluster hash: " + allTaxaHash);
        System.out.println("  Root cluster size: " + allTaxaHash.getSize());
        
        // Check if we have a choice for the root cluster
        CompactBipartition rootChoice = dpChoice.get(allTaxaHash);
        System.out.println("  Root choice available: " + (rootChoice != null ? "YES" : "NO"));
        if (rootChoice != null) {
            System.out.println("  Root choice: " + rootChoice);
        }
        
        Tree tree = new Tree();
        tree.taxaMap = geneTrees.taxaMap;
        
        System.out.println("\nBuilding tree structure...");
        try {
            tree.root = buildTreeNode(allTaxaHash, allTaxaSet, tree);
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
            
            System.out.println("Tree reconstruction completed successfully:");
            System.out.println("  Total nodes: " + tree.nodes.size());
            System.out.println("  Leaf nodes: " + tree.leavesCount);
            System.out.println("  Internal nodes: " + (tree.nodes.size() - tree.leavesCount));
            
            return tree;
            
        } catch (Exception e) {
            System.err.println("ERROR during tree reconstruction:");
            System.err.println("  Error message: " + e.getMessage());
            System.err.println("  Error type: " + e.getClass().getSimpleName());
            
            // Print diagnostic information
            System.err.println("\nDiagnostic information:");
            System.err.println("  Root cluster hash: " + allTaxaHash);
            System.err.println("  Available cluster hashes: " + clusterHashToCompactBips.keySet().size());
            System.err.println("  DP choices available: " + dpChoice.size());
            
            throw e;
        }
    }
    
    /**
     * Build tree node recursively using cluster hash and taxon set.
     */
    private TreeNode buildTreeNode(ClusterHash clusterHash, Set<Integer> taxonSet, Tree tree) {
        int taxaCount = clusterHash.getSize();
        
        if (taxaCount == 1) {
            int taxonId = taxonSet.iterator().next();
            return tree.addLeaf(geneTrees.taxa[taxonId]);
        }
        
        CompactBipartition choice = dpChoice.get(clusterHash);
        if (choice == null) {
            if (taxaCount == 2) {
                ArrayList<TreeNode> children = new ArrayList<>();
                for (int taxonId : taxonSet) {
                    children.add(tree.addLeaf(geneTrees.taxa[taxonId]));
                }
                return tree.addInternalNode(children);
            }
            
            // Handle clusters with no valid bipartitions by creating a polytomy (star tree)
            System.out.println("No bipartition choice for cluster " + clusterHash + " (size " + taxaCount + "), creating polytomy");
            ArrayList<TreeNode> children = new ArrayList<>();
            for (int taxonId : taxonSet) {
                children.add(tree.addLeaf(geneTrees.taxa[taxonId]));
            }
            return tree.addInternalNode(children);
        }
        
        // Create cluster hashes using the same method as DP
        ClusterHash leftHash = ClusterHash.fromLeftCluster(choice, hashFunction, prefixSums, prefixXORs, geneTreeOrderings);
        ClusterHash rightHash = ClusterHash.fromRightCluster(choice, hashFunction, prefixSums, prefixXORs, geneTreeOrderings);
        
        // Reconstruct taxon sets for tree building
        Set<Integer> leftTaxonSet = reconstructLeftTaxonSet(choice);
        Set<Integer> rightTaxonSet = reconstructRightTaxonSet(choice);
        
        TreeNode leftChild = buildTreeNode(leftHash, leftTaxonSet, tree);
        TreeNode rightChild = buildTreeNode(rightHash, rightTaxonSet, tree);
        
        ArrayList<TreeNode> children = new ArrayList<>();
        children.add(leftChild);
        children.add(rightChild);
        return tree.addInternalNode(children);
    }
    
    /**
     * Reconstruct left taxon set from compact bipartition.
     */
    private Set<Integer> reconstructLeftTaxonSet(CompactBipartition compactBip) {
        Set<Integer> taxonSet = new HashSet<>();
        int[] ordering = geneTreeOrderings[compactBip.geneTreeIndex];
        
        for (int i = compactBip.leftStart; i < compactBip.leftEnd && i < ordering.length; i++) {
            taxonSet.add(ordering[i]);
        }
        
        return taxonSet;
    }
    
    /**
     * Reconstruct right taxon set from compact bipartition.
     */
    private Set<Integer> reconstructRightTaxonSet(CompactBipartition compactBip) {
        Set<Integer> taxonSet = new HashSet<>();
        int[] ordering = geneTreeOrderings[compactBip.geneTreeIndex];
        
        for (int i = compactBip.rightStart; i < compactBip.rightEnd && i < ordering.length; i++) {
            taxonSet.add(ordering[i]);
        }
        
        return taxonSet;
    }
    
    /**
     * Print DP table for debugging.
     */
    public void printDPTable() {
        System.out.println("Memory-Efficient DP Results:");
        for (Map.Entry<ClusterHash, Double> entry : dpMemo.entrySet()) {
            ClusterHash clusterHash = entry.getKey();
            double score = entry.getValue();
            CompactBipartition choice = dpChoice.get(clusterHash);
            
            System.out.print("Cluster " + clusterHash + " -> Score: " + score);
            if (choice != null) {
                System.out.println(", Choice: " + choice);
            } else {
                System.out.println(", Choice: none");
            }
        }
    }
}
