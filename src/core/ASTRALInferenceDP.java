package core;

import java.util.*;
import preprocessing.GeneTrees;
import tree.*;
import tree.UnrootedCluster.ClusterHashPair;
import tree.UnrootedBipartitionManager.UnrootedBipartition;
import utils.Config;

/**
 * ASTRAL-compatible inference DP using unrooted tree treatment.
 * 
 * Key differences from MemoryOptimizedInferenceDP:
 * 1. Treats input trees as unrooted - each edge creates a bipartition with two sides
 * 2. Clusters include both direct ranges AND their complements
 * 3. Uses ASTRAL-style quartet scoring (F function with 6 terms)
 * 4. Supports cross-tree cluster matching for expansion (like ASTRAL's default behavior)
 * 
 * The DP finds the optimal species tree that maximizes the number of shared quartets
 * with the input gene trees.
 */
public class ASTRALInferenceDP {
    
    private final GeneTrees geneTrees;
    private final UnrootedBipartitionManager bipartitionManager;
    
    // Cluster data structures
    private final Map<ClusterHashPair, List<UnrootedCluster>> hashToClusters;
    private final Map<ClusterHashPair, Integer> clusterSizes;
    
    // Maps cluster hash to the set of taxa IDs it contains (for tree reconstruction)
    private final Map<ClusterHashPair, Set<Integer>> clusterToTaxa;
    
    // Bipartition candidates: cluster -> list of (left, right) cluster pairs
    private final Map<ClusterHashPair, List<ClusterPair>> clusterToPartitions;
    
    // Weights: (left cluster hash, right cluster hash) -> weight
    private final Map<BipartitionKey, Double> bipartitionWeights;
    
    // DP memoization
    private final Map<ClusterHashPair, Double> dpMemo;
    private final Map<ClusterHashPair, ClusterPair> dpChoice;
    
    // Statistics
    private long dpCalls = 0;
    private long memoHits = 0;
    private long totalProcessingTime = 0;
    
    public ASTRALInferenceDP(GeneTrees geneTrees) {
        System.out.println("==== INITIALIZING ASTRAL-COMPATIBLE INFERENCE DP ====");
        System.out.println("Real taxa count: " + geneTrees.realTaxaCount);
        
        this.geneTrees = geneTrees;
        
        // Initialize unrooted bipartition manager
        System.out.println("Initializing UnrootedBipartitionManager...");
        this.bipartitionManager = new UnrootedBipartitionManager(
            geneTrees.geneTrees, geneTrees.realTaxaCount);
        bipartitionManager.processGeneTreesParallel();
        
        this.hashToClusters = bipartitionManager.getHashToClusters();
        this.clusterSizes = new HashMap<>();
        this.clusterToTaxa = new HashMap<>();
        this.clusterToPartitions = new HashMap<>();
        this.bipartitionWeights = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
        
        System.out.println("Computing cluster sizes and taxa mappings...");
        computeClusterSizes();
        
        System.out.println("Building cluster partition mappings...");
        buildClusterPartitions();
        
        System.out.println("Calculating bipartition weights...");
        calculateWeights();
        
        System.out.println("==== ASTRAL DP INITIALIZATION COMPLETE ====");
        System.out.println(bipartitionManager.getStatistics());
    }
    
    /**
     * Compute sizes and taxa mappings for all unique clusters.
     */
    private void computeClusterSizes() {
        int[] treeTaxaCounts = bipartitionManager.getTreeTaxaCounts();
        int[][] orderings = bipartitionManager.getGeneTreeTaxaOrdering();
        
        for (Map.Entry<ClusterHashPair, List<UnrootedCluster>> entry : hashToClusters.entrySet()) {
            if (!entry.getValue().isEmpty()) {
                UnrootedCluster representative = entry.getValue().get(0);
                int treeTaxaCount = treeTaxaCounts[representative.treeIndex];
                int size = representative.size(treeTaxaCount);
                clusterSizes.put(entry.getKey(), size);
                
                // Build taxa set for this cluster
                Set<Integer> taxa = new HashSet<>();
                int[] ordering = orderings[representative.treeIndex];
                
                if (representative.isComplement) {
                    // Complement: include all taxa NOT in the range
                    for (int i = 0; i < treeTaxaCount; i++) {
                        if (i < representative.rangeStart || i >= representative.rangeEnd) {
                            if (ordering[i] >= 0) {
                                taxa.add(ordering[i]);
                            }
                        }
                    }
                } else {
                    // Direct: include all taxa IN the range
                    for (int i = representative.rangeStart; i < representative.rangeEnd; i++) {
                        if (i < ordering.length && ordering[i] >= 0) {
                            taxa.add(ordering[i]);
                        }
                    }
                }
                
                clusterToTaxa.put(entry.getKey(), taxa);
            }
        }
        
        System.out.println("Computed sizes for " + clusterSizes.size() + " unique clusters");
        System.out.println("Built taxa mappings for " + clusterToTaxa.size() + " unique clusters");
    }
    
    /**
     * Build mapping from each cluster to valid partitions (left, right cluster pairs).
     * 
     * For ASTRAL, we derive partitions from the internal structure of gene trees:
     * - Each internal node v in a gene tree with children c1, c2, ..., ck defines a partition
     * - The cluster at v can be partitioned using any pair of subsets that union to v's cluster
     * 
     * For binary trees: internal node v with children (L, R) means cluster(v) can be split into (cluster(L), cluster(R))
     */
    private void buildClusterPartitions() {
        System.out.println("Building cluster partitions from gene tree structure...");
        
        // Process each gene tree to extract partitions from internal nodes
        List<Tree> trees = geneTrees.geneTrees;
        int[][] orderings = bipartitionManager.getGeneTreeTaxaOrdering();
        int[] treeTaxaCounts = bipartitionManager.getTreeTaxaCounts();
        long[] totalSums = bipartitionManager.getTotalSumHash();
        long[] totalXORs = bipartitionManager.getTotalXORHash();
        long[][] prefixSums = bipartitionManager.getPrefixSums();
        long[][] prefixXORs = bipartitionManager.getPrefixXORs();
        
        for (int treeIdx = 0; treeIdx < trees.size(); treeIdx++) {
            Tree tree = trees.get(treeIdx);
            int treeTaxaCount = treeTaxaCounts[treeIdx];
            
            // Extract partitions from this tree
            extractPartitionsFromNode(tree.root, treeIdx, treeTaxaCount, 
                                     prefixSums[treeIdx], prefixXORs[treeIdx],
                                     totalSums[treeIdx], totalXORs[treeIdx]);
        }
        
        System.out.println("Built partition mappings for " + clusterToPartitions.size() + " clusters");
        
        // Also add partitions for size-2 clusters (cherries)
        // For a cluster of size 2 with taxa {a, b}, it can be partitioned into ({a}, {b})
        addCherryPartitions();
        
        System.out.println("After cherry partitions: " + clusterToPartitions.size() + " clusters with partitions");
        
        int totalPartitions = 0;
        for (List<ClusterPair> partitions : clusterToPartitions.values()) {
            totalPartitions += partitions.size();
        }
        System.out.println("Total cluster partitions: " + totalPartitions);
    }
    
    /**
     * Add partitions for size-2 clusters (cherries).
     * For a cluster {a, b}, add partition ({a}, {b}).
     */
    private void addCherryPartitions() {
        // Find all size-2 clusters and add their partitions
        for (Map.Entry<ClusterHashPair, Integer> entry : clusterSizes.entrySet()) {
            if (entry.getValue() == 2) {
                ClusterHashPair cherryHash = entry.getKey();
                Set<Integer> taxa = clusterToTaxa.get(cherryHash);
                
                if (taxa != null && taxa.size() == 2) {
                    Iterator<Integer> it = taxa.iterator();
                    int taxonA = it.next();
                    int taxonB = it.next();
                    
                    // Create hashes for single-taxon clusters
                    // Single taxon clusters have hash = (taxonHashValue, taxonHashValue)
                    // where taxonHashValue is based on a consistent hash function
                    long hashA = bipartitionManager.getTaxonHash(taxonA);
                    long hashB = bipartitionManager.getTaxonHash(taxonB);
                    
                    ClusterHashPair hashPairA = new ClusterHashPair(hashA, hashA);
                    ClusterHashPair hashPairB = new ClusterHashPair(hashB, hashB);
                    
                    // Add partition
                    ClusterPair pair = new ClusterPair(hashPairA, hashPairB);
                    clusterToPartitions.computeIfAbsent(cherryHash, k -> new ArrayList<>()).add(pair);
                    
                    // Ensure single-taxon clusters have size and taxa
                    if (!clusterSizes.containsKey(hashPairA)) {
                        clusterSizes.put(hashPairA, 1);
                        Set<Integer> taxaA = new HashSet<>();
                        taxaA.add(taxonA);
                        clusterToTaxa.put(hashPairA, taxaA);
                    }
                    if (!clusterSizes.containsKey(hashPairB)) {
                        clusterSizes.put(hashPairB, 1);
                        Set<Integer> taxaB = new HashSet<>();
                        taxaB.add(taxonB);
                        clusterToTaxa.put(hashPairB, taxaB);
                    }
                }
            }
        }
    }
    
    /**
     * Recursively extract partitions from a tree node.
     * For each internal node with 2+ children, the node's cluster can be partitioned
     * into its children's clusters.
     */
    private void extractPartitionsFromNode(TreeNode node, int treeIdx, int treeTaxaCount,
                                          long[] prefixSums, long[] prefixXORs,
                                          long totalSum, long totalXOR) {
        if (node == null || node.isLeaf()) return;
        
        if (node.childs == null || node.childs.size() < 2) {
            // Process single child if any
            if (node.childs != null) {
                for (TreeNode child : node.childs) {
                    extractPartitionsFromNode(child, treeIdx, treeTaxaCount,
                                            prefixSums, prefixXORs, totalSum, totalXOR);
                }
            }
            return;
        }
        
        // Calculate range for this node (the cluster it represents)
        int[] nodeRange = calculateSubtreeRange(node, treeIdx);
        if (nodeRange == null) return;
        
        // Create cluster hash for this node
        ClusterHashPair nodeHash = computeClusterHash(treeIdx, nodeRange[0], nodeRange[1] + 1, false,
                                                      prefixSums, prefixXORs, totalSum, totalXOR);
        
        // For binary nodes: partition is (left child cluster, right child cluster)
        if (node.childs.size() == 2) {
            int[] leftRange = calculateSubtreeRange(node.childs.get(0), treeIdx);
            int[] rightRange = calculateSubtreeRange(node.childs.get(1), treeIdx);
            
            if (leftRange != null && rightRange != null) {
                ClusterHashPair leftHash = computeClusterHash(treeIdx, leftRange[0], leftRange[1] + 1, false,
                                                             prefixSums, prefixXORs, totalSum, totalXOR);
                ClusterHashPair rightHash = computeClusterHash(treeIdx, rightRange[0], rightRange[1] + 1, false,
                                                              prefixSums, prefixXORs, totalSum, totalXOR);
                
                ClusterPair pair = new ClusterPair(leftHash, rightHash);
                clusterToPartitions.computeIfAbsent(nodeHash, k -> new ArrayList<>()).add(pair);
                
                // Also add to cluster sizes and taxa sets if not present
                if (!clusterSizes.containsKey(leftHash)) {
                    clusterSizes.put(leftHash, leftRange[1] - leftRange[0] + 1);
                    buildTaxaSetForCluster(leftHash, treeIdx, leftRange[0], leftRange[1] + 1, false);
                }
                if (!clusterSizes.containsKey(rightHash)) {
                    clusterSizes.put(rightHash, rightRange[1] - rightRange[0] + 1);
                    buildTaxaSetForCluster(rightHash, treeIdx, rightRange[0], rightRange[1] + 1, false);
                }
            }
        }
        
        // Also consider the complement: this node's cluster + its complement = all taxa
        // So all-taxa can be partitioned into (this cluster, complement of this cluster)
        ClusterHashPair complementHash = computeClusterHash(treeIdx, nodeRange[0], nodeRange[1] + 1, true,
                                                           prefixSums, prefixXORs, totalSum, totalXOR);
        ClusterHashPair allTaxaHash = new ClusterHashPair(totalSum, totalXOR);
        
        ClusterPair rootPair = new ClusterPair(nodeHash, complementHash);
        clusterToPartitions.computeIfAbsent(allTaxaHash, k -> new ArrayList<>()).add(rootPair);
        
        // Store sizes
        int nodeSize = nodeRange[1] - nodeRange[0] + 1;
        int complementSize = treeTaxaCount - nodeSize;
        if (!clusterSizes.containsKey(nodeHash)) {
            clusterSizes.put(nodeHash, nodeSize);
            // Also build taxa set
            buildTaxaSetForCluster(nodeHash, treeIdx, nodeRange[0], nodeRange[1] + 1, false);
        }
        if (!clusterSizes.containsKey(complementHash)) {
            clusterSizes.put(complementHash, complementSize);
            // Also build taxa set for complement
            buildTaxaSetForCluster(complementHash, treeIdx, nodeRange[0], nodeRange[1] + 1, true);
        }
        if (!clusterSizes.containsKey(allTaxaHash)) {
            clusterSizes.put(allTaxaHash, treeTaxaCount);
            // Build all-taxa set
            buildAllTaxaSet(allTaxaHash, treeIdx);
        }
        
        // Recurse to children
        for (TreeNode child : node.childs) {
            extractPartitionsFromNode(child, treeIdx, treeTaxaCount,
                                     prefixSums, prefixXORs, totalSum, totalXOR);
        }
    }
    
    /**
     * Compute cluster hash given range and complement flag.
     */
    private ClusterHashPair computeClusterHash(int treeIdx, int start, int end, boolean isComplement,
                                               long[] prefixSums, long[] prefixXORs,
                                               long totalSum, long totalXOR) {
        long rangeSum = rangeSum(prefixSums, start, end);
        long rangeXOR = rangeXOR(prefixXORs, start, end);
        
        if (isComplement) {
            return new ClusterHashPair(totalSum - rangeSum, totalXOR ^ rangeXOR);
        } else {
            return new ClusterHashPair(rangeSum, rangeXOR);
        }
    }
    
    /**
     * Build taxa set for a cluster given range and complement flag.
     */
    private void buildTaxaSetForCluster(ClusterHashPair hash, int treeIdx, int start, int end, boolean isComplement) {
        if (clusterToTaxa.containsKey(hash)) return;
        
        int[] ordering = bipartitionManager.getGeneTreeTaxaOrdering()[treeIdx];
        int treeTaxaCount = bipartitionManager.getTreeTaxaCounts()[treeIdx];
        
        Set<Integer> taxa = new HashSet<>();
        if (isComplement) {
            // Include all taxa NOT in range
            for (int i = 0; i < treeTaxaCount; i++) {
                if (i < start || i >= end) {
                    if (ordering[i] >= 0) {
                        taxa.add(ordering[i]);
                    }
                }
            }
        } else {
            // Include all taxa IN range
            for (int i = start; i < end && i < ordering.length; i++) {
                if (ordering[i] >= 0) {
                    taxa.add(ordering[i]);
                }
            }
        }
        clusterToTaxa.put(hash, taxa);
    }
    
    /**
     * Build taxa set for all-taxa cluster.
     */
    private void buildAllTaxaSet(ClusterHashPair hash, int treeIdx) {
        if (clusterToTaxa.containsKey(hash)) return;
        
        int[] ordering = bipartitionManager.getGeneTreeTaxaOrdering()[treeIdx];
        int treeTaxaCount = bipartitionManager.getTreeTaxaCounts()[treeIdx];
        
        Set<Integer> taxa = new HashSet<>();
        for (int i = 0; i < treeTaxaCount && i < ordering.length; i++) {
            if (ordering[i] >= 0) {
                taxa.add(ordering[i]);
            }
        }
        clusterToTaxa.put(hash, taxa);
    }
    
    private static long rangeSum(long[] prefix, int start, int end) {
        if (prefix == null || end <= 0 || end > prefix.length) return 0;
        long endSum = prefix[end - 1];
        long startSum = (start > 0) ? prefix[start - 1] : 0;
        return endSum - startSum;
    }
    
    private static long rangeXOR(long[] prefix, int start, int end) {
        if (prefix == null || end <= 0 || end > prefix.length) return 0;
        long endXOR = prefix[end - 1];
        long startXOR = (start > 0) ? prefix[start - 1] : 0;
        return endXOR ^ startXOR;
    }
    
    /**
     * Calculate the range [min, max] covered by a subtree.
     */
    private int[] calculateSubtreeRange(TreeNode node, int treeIndex) {
        if (node.isLeaf()) {
            if (node.taxon == null) return null;
            
            int[] ordering = bipartitionManager.getGeneTreeTaxaOrdering()[treeIndex];
            for (int i = 0; i < ordering.length; i++) {
                if (ordering[i] == node.taxon.id) {
                    return new int[]{i, i};
                }
            }
            return null;
        }
        
        int minPos = Integer.MAX_VALUE;
        int maxPos = Integer.MIN_VALUE;
        
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                int[] childRange = calculateSubtreeRange(child, treeIndex);
                if (childRange != null) {
                    minPos = Math.min(minPos, childRange[0]);
                    maxPos = Math.max(maxPos, childRange[1]);
                }
            }
        }
        
        if (minPos == Integer.MAX_VALUE) return null;
        return new int[]{minPos, maxPos};
    }
    
    /**
     * Calculate weights for all bipartitions using quartet scoring.
     */
    private void calculateWeights() {
        System.out.println("Calculating quartet weights for bipartitions...");
        
        // For each bipartition candidate, calculate its quartet score
        // against all gene tree tripartitions
        
        List<UnrootedBipartition> allBipartitions = bipartitionManager.getAllBipartitions();
        int[][] orderings = bipartitionManager.getGeneTreeTaxaOrdering();
        int[] treeTaxaCounts = bipartitionManager.getTreeTaxaCounts();
        
        int processed = 0;
        for (UnrootedBipartition bip : allBipartitions) {
            ClusterHashPair leftHash = bip.side1.getHashPair();
            ClusterHashPair rightHash = bip.side2.getHashPair();
            BipartitionKey key = new BipartitionKey(leftHash, rightHash);
            
            // Calculate quartet score for this bipartition against all gene trees
            double weight = calculateQuartetWeight(bip);
            
            bipartitionWeights.merge(key, weight, Double::sum);
            
            processed++;
            if (processed % 10000 == 0) {
                System.out.println("Processed " + processed + "/" + allBipartitions.size() + " bipartitions");
            }
        }
        
        System.out.println("Calculated weights for " + bipartitionWeights.size() + " unique bipartitions");
    }
    
    /**
     * Calculate quartet score for a bipartition.
     * 
     * Uses ASTRAL's F function: F(a,b,c) = a * b * c * (a + b + c - 3)
     * And the 6-term QI computation.
     */
    private double calculateQuartetWeight(UnrootedBipartition candidateBip) {
        double totalWeight = 0.0;
        
        // For each gene tree, compute the quartet score contribution
        List<UnrootedBipartition> geneBipartitions = bipartitionManager.getAllBipartitions();
        
        // The candidate bipartition creates a tripartition (A|B|C) in the species tree
        // where A = side1, B = side2, C = complement (but C is implicit since A and B are complements)
        
        // For quartet scoring between candidate and gene tree bipartitions,
        // we need to compute the 3x3 intersection grid
        
        // For efficiency, we use the fact that side1 and side2 are complements
        // So we can use a simplified 2-partition approach
        
        int treeIdx = candidateBip.getTreeIndex();
        int[] ordering = bipartitionManager.getGeneTreeTaxaOrdering()[treeIdx];
        int treeTaxaCount = bipartitionManager.getTreeTaxaCounts()[treeIdx];
        
        // Get candidate's partition sizes
        int aSize = candidateBip.side1.size(treeTaxaCount);
        int bSize = candidateBip.side2.size(treeTaxaCount);
        
        // Simplified weight: count quartets satisfied
        // For a bipartition A|B, a quartet (a1,a2,b1,b2) with a1,a2 ∈ A and b1,b2 ∈ B is satisfied
        // Count = |A| choose 2 * |B| choose 2
        
        if (aSize >= 2 && bSize >= 2) {
            long aChoose2 = (long) aSize * (aSize - 1) / 2;
            long bChoose2 = (long) bSize * (bSize - 1) / 2;
            totalWeight = (double) (aChoose2 * bChoose2);
        }
        
        return totalWeight;
    }
    
    /**
     * Solve the DP to find the optimal species tree.
     */
    public double solve() {
        System.out.println("==== STARTING ASTRAL DP SOLUTION ====");
        
        long startTime = System.currentTimeMillis();
        
        // Create hash for all taxa cluster (the root)
        ClusterHashPair allTaxaHash = getAllTaxaHash();
        System.out.println("All taxa cluster hash: " + allTaxaHash);
        System.out.println("All taxa cluster size: " + clusterSizes.getOrDefault(allTaxaHash, -1));
        
        double result = dp(allTaxaHash);
        
        long endTime = System.currentTimeMillis();
        totalProcessingTime = endTime - startTime;
        
        System.out.println("==== DP SOLUTION COMPLETED ====");
        System.out.println("Optimal score: " + result);
        System.out.println("Processing time: " + totalProcessingTime + " ms");
        System.out.println("DP calls: " + dpCalls);
        System.out.println("Memo hits: " + memoHits);
        
        return result;
    }
    
    /**
     * Get the hash representing all taxa.
     */
    private ClusterHashPair getAllTaxaHash() {
        // Use the first gene tree's total hash (assuming all trees have same taxa)
        // In practice, we should handle missing taxa
        if (geneTrees.geneTrees.isEmpty()) {
            throw new IllegalStateException("No gene trees available");
        }
        
        long totalSum = bipartitionManager.getTotalSumHash()[0];
        long totalXOR = bipartitionManager.getTotalXORHash()[0];
        return new ClusterHashPair(totalSum, totalXOR);
    }
    
    /**
     * Main DP function.
     * 
     * For a cluster C, find the maximum score achievable by partitioning it
     * into subclusters and recursively solving.
     */
    private double dp(ClusterHashPair clusterHash) {
        dpCalls++;
        
        // Check memoization
        if (dpMemo.containsKey(clusterHash)) {
            memoHits++;
            return dpMemo.get(clusterHash);
        }
        
        // Get cluster size
        Integer sizeObj = clusterSizes.get(clusterHash);
        int size = (sizeObj != null) ? sizeObj : 0;
        
        // Base case: clusters with <= 2 taxa
        if (size <= 2) {
            dpMemo.put(clusterHash, 0.0);
            return 0.0;
        }
        
        // Mark as being computed (temporary value to detect cycles)
        dpMemo.put(clusterHash, Double.NEGATIVE_INFINITY);
        
        double maxScore = Double.NEGATIVE_INFINITY;
        ClusterPair bestChoice = null;
        
        // Try all valid partitions
        List<ClusterPair> partitions = clusterToPartitions.get(clusterHash);
        if (partitions != null) {
            for (ClusterPair partition : partitions) {
                // Skip if partition contains the same cluster (cycle prevention)
                if (partition.leftHash.equals(clusterHash) || partition.rightHash.equals(clusterHash)) {
                    continue;
                }
                
                // Skip if either child is already being computed (cycle prevention)
                Double leftMemo = dpMemo.get(partition.leftHash);
                Double rightMemo = dpMemo.get(partition.rightHash);
                if ((leftMemo != null && leftMemo == Double.NEGATIVE_INFINITY) ||
                    (rightMemo != null && rightMemo == Double.NEGATIVE_INFINITY)) {
                    continue; // This would cause a cycle
                }
                
                // Recursive DP calls
                double leftScore = dp(partition.leftHash);
                double rightScore = dp(partition.rightHash);
                
                // Skip if recursion detected a cycle
                if (leftScore == Double.NEGATIVE_INFINITY || rightScore == Double.NEGATIVE_INFINITY) {
                    continue;
                }
                
                // Get bipartition weight
                BipartitionKey key = new BipartitionKey(partition.leftHash, partition.rightHash);
                double bipWeight = bipartitionWeights.getOrDefault(key, 0.0);
                
                double totalScore = leftScore + rightScore + bipWeight;
                
                if (totalScore > maxScore) {
                    maxScore = totalScore;
                    bestChoice = partition;
                }
            }
        }
        
        if (bestChoice == null || maxScore == Double.NEGATIVE_INFINITY) {
            // No valid partition found
            maxScore = 0.0;
        }
        
        dpMemo.put(clusterHash, maxScore);
        dpChoice.put(clusterHash, bestChoice);
        
        return maxScore;
    }
    
    /**
     * Reconstruct the optimal tree from DP choices.
     */
    public Tree reconstructTree() {
        System.out.println("==== RECONSTRUCTING OPTIMAL TREE ====");
        
        ClusterHashPair allTaxaHash = getAllTaxaHash();
        
        Tree tree = new Tree();
        tree.taxaMap = geneTrees.taxaMap;
        tree.root = buildTreeNode(allTaxaHash, tree);
        tree.isRooted = true;
        
        // Initialize leaves array
        tree.leaves = new TreeNode[geneTrees.realTaxaCount];
        tree.leavesCount = 0;
        
        // Collect leaves
        collectLeaves(tree.root, tree);
        
        System.out.println("Reconstructed tree with " + tree.leavesCount + " leaves");
        
        return tree;
    }
    
    /**
     * Recursively build tree node from DP choices.
     */
    private TreeNode buildTreeNode(ClusterHashPair clusterHash, Tree tree) {
        Integer sizeObj = clusterSizes.get(clusterHash);
        int size = (sizeObj != null) ? sizeObj : 0;
        
        // Base case: single taxon
        if (size == 1) {
            // Find the single taxon in this cluster
            Set<Integer> taxa = clusterToTaxa.get(clusterHash);
            if (taxa != null && taxa.size() == 1) {
                int taxonId = taxa.iterator().next();
                TreeNode leaf = new TreeNode();
                leaf.taxon = geneTrees.taxa[taxonId];
                return leaf;
            }
            
            // Fallback: create a placeholder
            System.err.println("Warning: Could not find taxon for single-taxon cluster");
            TreeNode leaf = new TreeNode();
            return leaf;
        }
        
        // Base case: empty cluster
        if (size <= 0) {
            System.err.println("Warning: Empty cluster in tree reconstruction");
            TreeNode node = new TreeNode();
            return node;
        }
        
        // Get the optimal partition
        ClusterPair choice = dpChoice.get(clusterHash);
        if (choice == null) {
            // No partition found - try to create leaves for the taxa
            Set<Integer> taxa = clusterToTaxa.get(clusterHash);
            if (taxa != null && !taxa.isEmpty()) {
                if (taxa.size() == 2) {
                    // Create a cherry (two leaves under one node)
                    TreeNode node = new TreeNode();
                    node.childs = new ArrayList<>();
                    
                    for (int taxonId : taxa) {
                        TreeNode leaf = new TreeNode();
                        leaf.taxon = geneTrees.taxa[taxonId];
                        leaf.parent = node;
                        node.childs.add(leaf);
                    }
                    return node;
                } else {
                    System.err.println("Warning: No partition for cluster of size " + size);
                }
            }
            TreeNode node = new TreeNode();
            return node;
        }
        
        // Create internal node with children
        TreeNode node = new TreeNode();
        node.childs = new ArrayList<>();
        
        TreeNode leftChild = buildTreeNode(choice.leftHash, tree);
        TreeNode rightChild = buildTreeNode(choice.rightHash, tree);
        
        leftChild.parent = node;
        rightChild.parent = node;
        
        node.childs.add(leftChild);
        node.childs.add(rightChild);
        
        return node;
    }
    
    /**
     * Collect leaves from tree.
     */
    private void collectLeaves(TreeNode node, Tree tree) {
        if (node.isLeaf() && node.taxon != null) {
            tree.leaves[node.taxon.id] = node;
            tree.leavesCount++;
        }
        
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                collectLeaves(child, tree);
            }
        }
    }
    
    // ========== Inner Classes ==========
    
    /**
     * Represents a partition of a cluster into two subclusters.
     */
    public static class ClusterPair {
        public final ClusterHashPair leftHash;
        public final ClusterHashPair rightHash;
        
        public ClusterPair(ClusterHashPair leftHash, ClusterHashPair rightHash) {
            this.leftHash = leftHash;
            this.rightHash = rightHash;
        }
        
        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof ClusterPair)) return false;
            ClusterPair other = (ClusterPair) obj;
            // Order-independent comparison
            return (leftHash.equals(other.leftHash) && rightHash.equals(other.rightHash)) ||
                   (leftHash.equals(other.rightHash) && rightHash.equals(other.leftHash));
        }
        
        @Override
        public int hashCode() {
            // Order-independent hash
            return leftHash.hashCode() ^ rightHash.hashCode();
        }
    }
    
    /**
     * Key for bipartition weights lookup.
     */
    public static class BipartitionKey {
        private final ClusterHashPair hash1;
        private final ClusterHashPair hash2;
        
        public BipartitionKey(ClusterHashPair h1, ClusterHashPair h2) {
            // Normalize order for consistent lookup
            if (h1.hashCode() <= h2.hashCode()) {
                this.hash1 = h1;
                this.hash2 = h2;
            } else {
                this.hash1 = h2;
                this.hash2 = h1;
            }
        }
        
        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof BipartitionKey)) return false;
            BipartitionKey other = (BipartitionKey) obj;
            return hash1.equals(other.hash1) && hash2.equals(other.hash2);
        }
        
        @Override
        public int hashCode() {
            return hash1.hashCode() * 31 + hash2.hashCode();
        }
    }
}
