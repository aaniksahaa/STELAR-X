package tree;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import tree.UnrootedCluster.ClusterHashPair;
import utils.Threading;

/**
 * UnrootedBipartitionManager: Manages cluster extraction for ASTRAL-style unrooted tree treatment.
 * 
 * Key differences from MemoryEfficientBipartitionManager:
 * 1. Treats trees as unrooted - each edge creates a bipartition with two sides
 * 2. For each internal node, registers BOTH the subtree cluster AND its complement
 * 3. Supports cross-tree cluster matching for ASTRAL-style expansion
 * 
 * Example: For tree ((A,(B,C)),(D,E)) with edge between BC and rest:
 * - Registers cluster {B,C} (direct range)
 * - Registers cluster {A,D,E} (complement of {B,C})
 */
public class UnrootedBipartitionManager {
    
    // Gene tree data structures
    private final List<Tree> geneTrees;
    private final int realTaxaCount;
    
    // Per-tree orderings and prefix arrays
    private final int[][] geneTreeTaxaOrdering;   // [tree_index][position] = taxon_id
    private final long[][] prefixSums;            // [tree_index][position] = prefix sum
    private final long[][] prefixXORs;            // [tree_index][position] = prefix XOR
    private final long[] totalSumHash;            // [tree_index] = total sum hash for tree
    private final long[] totalXORHash;            // [tree_index] = total XOR hash for tree
    private final int[] treeTaxaCounts;           // [tree_index] = number of taxa in tree
    
    // Cluster storage
    private final Map<ClusterHashPair, List<UnrootedCluster>> hashToClusters;
    private final List<UnrootedCluster> allClusters;
    
    // Bipartition storage (pairs of clusters representing edges)
    private final List<UnrootedBipartition> allBipartitions;
    
    // Statistics
    private int totalClustersExtracted = 0;
    private int uniqueClusterHashes = 0;
    
    public UnrootedBipartitionManager(List<Tree> geneTrees, int realTaxaCount) {
        this.geneTrees = geneTrees;
        this.realTaxaCount = realTaxaCount;
        
        int numTrees = geneTrees.size();
        this.geneTreeTaxaOrdering = new int[numTrees][];
        this.prefixSums = new long[numTrees][];
        this.prefixXORs = new long[numTrees][];
        this.totalSumHash = new long[numTrees];
        this.totalXORHash = new long[numTrees];
        this.treeTaxaCounts = new int[numTrees];
        
        this.hashToClusters = new ConcurrentHashMap<>();
        this.allClusters = Collections.synchronizedList(new ArrayList<>());
        this.allBipartitions = Collections.synchronizedList(new ArrayList<>());
        
        initializeGeneTreeOrderings();
        calculatePrefixArrays();
    }
    
    /**
     * Initialize the left-to-right ordering for each gene tree.
     */
    private void initializeGeneTreeOrderings() {
        System.out.println("Initializing gene tree taxa orderings for unrooted treatment...");
        
        for (int i = 0; i < geneTrees.size(); i++) {
            Tree tree = geneTrees.get(i);
            List<Integer> ordering = new ArrayList<>();
            collectLeavesInOrder(tree.root, ordering);
            
            geneTreeTaxaOrdering[i] = ordering.stream().mapToInt(Integer::intValue).toArray();
            treeTaxaCounts[i] = ordering.size();
        }
        
        System.out.println("Initialized orderings for " + geneTrees.size() + " gene trees");
    }
    
    /**
     * Collect leaves in left-to-right order (inorder traversal).
     */
    private void collectLeavesInOrder(TreeNode node, List<Integer> ordering) {
        if (node.isLeaf()) {
            if (node.taxon != null) {
                ordering.add(node.taxon.id);
            }
            return;
        }
        
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                collectLeavesInOrder(child, ordering);
            }
        }
    }
    
    /**
     * Calculate prefix arrays and total hashes for efficient complement hash computation.
     */
    private void calculatePrefixArrays() {
        System.out.println("Calculating prefix arrays for cluster hashing...");
        
        for (int i = 0; i < geneTrees.size(); i++) {
            int[] ordering = geneTreeTaxaOrdering[i];
            int n = ordering.length;
            
            prefixSums[i] = new long[n];
            prefixXORs[i] = new long[n];
            
            if (n > 0) {
                long hashedTaxon0 = hashSingleTaxon(ordering[0]);
                prefixSums[i][0] = hashedTaxon0;
                prefixXORs[i][0] = hashedTaxon0;
                
                for (int j = 1; j < n; j++) {
                    long hashedTaxon = hashSingleTaxon(ordering[j]);
                    prefixSums[i][j] = prefixSums[i][j - 1] + hashedTaxon;
                    prefixXORs[i][j] = prefixXORs[i][j - 1] ^ hashedTaxon;
                }
                
                // Store total hashes for complement computation
                totalSumHash[i] = prefixSums[i][n - 1];
                totalXORHash[i] = prefixXORs[i][n - 1];
            }
        }
        
        System.out.println("Calculated prefix arrays for " + geneTrees.size() + " trees");
    }
    
    /**
     * Hash function for individual taxon IDs.
     */
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
     * Public method to get the hash value for a single taxon.
     * Used for creating single-taxon cluster hashes.
     */
    public long getTaxonHash(int taxonId) {
        return hashSingleTaxon(taxonId);
    }
    
    /**
     * Process all gene trees to extract clusters and bipartitions.
     * This treats trees as unrooted - each node generates both a cluster and its complement.
     */
    public void processGeneTreesParallel() {
        System.out.println("Extracting clusters from gene trees (unrooted treatment)...");
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        int chunkSize = Math.max(1, (geneTrees.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (geneTrees.size() + chunkSize - 1) / chunkSize);
        
        CountDownLatch latch = new CountDownLatch(actualThreads);
        AtomicInteger processedTrees = new AtomicInteger(0);
        AtomicInteger totalClusters = new AtomicInteger(0);
        
        System.out.println("Using " + actualThreads + " threads for parallel processing");
        
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, geneTrees.size());
            final int threadId = i;
            
            if (startIdx >= geneTrees.size()) {
                latch.countDown();
                continue;
            }
            
            Threading.execute(() -> {
                try {
                    List<UnrootedCluster> localClusters = new ArrayList<>();
                    List<UnrootedBipartition> localBipartitions = new ArrayList<>();
                    
                    for (int treeIdx = startIdx; treeIdx < endIdx; treeIdx++) {
                        Tree tree = geneTrees.get(treeIdx);
                        extractClustersFromTree(tree.root, treeIdx, localClusters, localBipartitions);
                        processedTrees.incrementAndGet();
                    }
                    
                    // Compute hashes for all local clusters
                    for (UnrootedCluster cluster : localClusters) {
                        cluster.computeHash(
                            prefixSums[cluster.treeIndex],
                            prefixXORs[cluster.treeIndex],
                            totalSumHash[cluster.treeIndex],
                            totalXORHash[cluster.treeIndex]
                        );
                    }
                    
                    // Add to global collections
                    allClusters.addAll(localClusters);
                    allBipartitions.addAll(localBipartitions);
                    totalClusters.addAndGet(localClusters.size());
                    
                    // Group by hash
                    for (UnrootedCluster cluster : localClusters) {
                        ClusterHashPair hash = cluster.getHashPair();
                        hashToClusters.computeIfAbsent(hash, k -> 
                            Collections.synchronizedList(new ArrayList<>())).add(cluster);
                    }
                    
                    System.out.println("Thread " + threadId + " processed " + 
                                      (endIdx - startIdx) + " trees, extracted " + 
                                      localClusters.size() + " clusters");
                    
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Cluster extraction interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        totalClustersExtracted = totalClusters.get();
        uniqueClusterHashes = hashToClusters.size();
        
        System.out.println("\n=== Unrooted Cluster Extraction Complete ===");
        System.out.println("Total clusters extracted: " + totalClustersExtracted);
        System.out.println("Unique cluster hashes: " + uniqueClusterHashes);
        System.out.println("Total bipartitions: " + allBipartitions.size());
    }
    
    /**
     * Extract clusters from a tree node recursively.
     * For each internal node (including root), we register:
     * 1. The subtree cluster (direct range)
     * 2. The complement cluster (complement of range)
     * 
     * For leaf nodes, we also register the single-taxon cluster and its complement.
     */
    private void extractClustersFromTree(TreeNode node, int treeIndex,
                                         List<UnrootedCluster> clusters,
                                         List<UnrootedBipartition> bipartitions) {
        // Calculate range for this subtree
        int[] range = calculateSubtreeRange(node, treeIndex);
        if (range == null) return;
        
        int start = range[0];
        int end = range[1] + 1; // Convert to exclusive end
        int treeTaxaCount = treeTaxaCounts[treeIndex];
        
        // Skip if this covers all taxa (complement would be empty)
        if (end - start >= treeTaxaCount) {
            // Still process children
            if (node.childs != null) {
                for (TreeNode child : node.childs) {
                    extractClustersFromTree(child, treeIndex, clusters, bipartitions);
                }
            }
            return;
        }
        
        // Skip if range is empty
        if (start >= end) {
            if (node.childs != null) {
                for (TreeNode child : node.childs) {
                    extractClustersFromTree(child, treeIndex, clusters, bipartitions);
                }
            }
            return;
        }
        
        // Create both the direct cluster and its complement
        UnrootedCluster directCluster = UnrootedCluster.fromRange(treeIndex, start, end);
        UnrootedCluster complementCluster = UnrootedCluster.complementOf(treeIndex, start, end);
        
        clusters.add(directCluster);
        clusters.add(complementCluster);
        
        // Create a bipartition representing this edge
        UnrootedBipartition bip = new UnrootedBipartition(directCluster, complementCluster);
        bipartitions.add(bip);
        
        // Recursively process children
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                extractClustersFromTree(child, treeIndex, clusters, bipartitions);
            }
        }
    }
    
    /**
     * Calculate the range [min, max] covered by a subtree.
     */
    private int[] calculateSubtreeRange(TreeNode node, int treeIndex) {
        if (node.isLeaf()) {
            if (node.taxon == null) return null;
            
            int[] ordering = geneTreeTaxaOrdering[treeIndex];
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
    
    // ========== Getters ==========
    
    public Map<ClusterHashPair, List<UnrootedCluster>> getHashToClusters() {
        return hashToClusters;
    }
    
    public List<UnrootedCluster> getAllClusters() {
        return allClusters;
    }
    
    public List<UnrootedBipartition> getAllBipartitions() {
        return allBipartitions;
    }
    
    public int[][] getGeneTreeTaxaOrdering() {
        return geneTreeTaxaOrdering;
    }
    
    public long[][] getPrefixSums() {
        return prefixSums;
    }
    
    public long[][] getPrefixXORs() {
        return prefixXORs;
    }
    
    public long[] getTotalSumHash() {
        return totalSumHash;
    }
    
    public long[] getTotalXORHash() {
        return totalXORHash;
    }
    
    public int[] getTreeTaxaCounts() {
        return treeTaxaCounts;
    }
    
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("UnrootedBipartitionManager Statistics:\n");
        sb.append("  Gene trees processed: ").append(geneTrees.size()).append("\n");
        sb.append("  Total clusters extracted: ").append(totalClustersExtracted).append("\n");
        sb.append("  Unique cluster hashes: ").append(uniqueClusterHashes).append("\n");
        sb.append("  Total bipartitions: ").append(allBipartitions.size()).append("\n");
        sb.append("  Average clusters per hash: ").append(
            uniqueClusterHashes > 0 ? 
            String.format("%.2f", (double)totalClustersExtracted / uniqueClusterHashes) : "N/A"
        ).append("\n");
        return sb.toString();
    }
    
    /**
     * Inner class representing a bipartition (edge) in the unrooted tree.
     * A bipartition consists of two complementary clusters.
     */
    public static class UnrootedBipartition {
        public final UnrootedCluster side1;
        public final UnrootedCluster side2;
        
        public UnrootedBipartition(UnrootedCluster side1, UnrootedCluster side2) {
            this.side1 = side1;
            this.side2 = side2;
        }
        
        public int getTreeIndex() {
            return side1.treeIndex;
        }
        
        @Override
        public String toString() {
            return "UnrootedBipartition[" + side1 + " | " + side2 + "]";
        }
    }
}
