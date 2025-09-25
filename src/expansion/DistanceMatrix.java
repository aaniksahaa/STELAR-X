package expansion;

import java.util.*;
import java.util.concurrent.CountDownLatch;
import preprocessing.GeneTrees;
import tree.STBipartition;
import tree.Tree;
import tree.TreeNode;
import utils.BitSet;
import utils.BipartitionExpansionConfig;
import utils.Threading;

/**
 * Distance matrix computation and tree inference for bipartition expansion.
 * 
 * This class implements distance-based methods similar to ASTRAL's distance matrix approach.
 * It calculates pairwise distances between taxa based on gene tree topologies and
 * uses these distances to infer species trees via UPGMA or Neighbor-Joining.
 */
public class DistanceMatrix {
    
    private GeneTrees geneTrees;
    private double[][] distanceMatrix;
    private int taxaCount;
    private String[] taxaNames;
    
    public DistanceMatrix(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
        this.taxaCount = geneTrees.realTaxaCount;
        this.taxaNames = geneTrees.taxonIdToLabel;
        this.distanceMatrix = new double[taxaCount][taxaCount];
        
        calculateDistanceMatrix();
    }
    
    /**
     * Calculate pairwise distances between taxa based on gene tree topologies.
     * Distance is based on the frequency of taxa appearing in the same clade
     * across all gene trees.
     */
    private void calculateDistanceMatrix() {
        calculateDistanceMatrixParallel();
    }
    
    /**
     * Parallel version of distance matrix calculation for improved performance.
     * Divides the gene trees among multiple threads for concurrent processing.
     */
    private void calculateDistanceMatrixParallel() {
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("Calculating distance matrix from " + geneTrees.geneTrees.size() + " gene trees in parallel...");
        }
        
        // Initialize distance matrix
        for (int i = 0; i < taxaCount; i++) {
            for (int j = 0; j < taxaCount; j++) {
                distanceMatrix[i][j] = (i == j) ? 0.0 : 1.0; // Start with maximum distance
            }
        }
        
        // Thread-safe co-occurrence counting
        int[][] coOccurrenceCount = new int[taxaCount][taxaCount];
        List<Tree> geneTrees = this.geneTrees.geneTrees;
        
        if (geneTrees.isEmpty()) {
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("No gene trees found, distance matrix remains at default values.");
            }
            return;
        }
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        // Calculate optimal number of threads to avoid invalid ranges
        int chunkSize = Math.max(1, (geneTrees.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (geneTrees.size() + chunkSize - 1) / chunkSize);
        CountDownLatch latch = new CountDownLatch(actualThreads);
        
        // Process gene trees in parallel chunks
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, geneTrees.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    // Validate range before processing
                    if (startIdx >= endIdx || startIdx >= geneTrees.size()) {
                        return;
                    }
                    
                    int[][] localCoOccurrenceCount = new int[taxaCount][taxaCount];
                    
                    if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                        System.out.println("Thread " + threadId + " processing gene trees " + startIdx + " to " + (endIdx-1));
                    }
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        Tree geneTree = geneTrees.get(j);
                        Map<Integer, BitSet> taxonToClade = new HashMap<>();
                        
                        // For each gene tree, find which clade each taxon belongs to
                        collectClades(geneTree.root, taxonToClade);
                        
                        // Count co-occurrences in this tree
                        for (int x = 0; x < taxaCount; x++) {
                            for (int y = x + 1; y < taxaCount; y++) {
                                BitSet cladeX = taxonToClade.get(x);
                                BitSet cladeY = taxonToClade.get(y);
                                
                                if (cladeX != null && cladeY != null) {
                                    // Check if taxa x and y are in the same clade
                                    if (cladeX.equals(cladeY) || isInSameClade(x, y, geneTree.root)) {
                                        localCoOccurrenceCount[x][y]++;
                                        localCoOccurrenceCount[y][x]++;
                                    }
                                }
                            }
                        }
                    }
                    
                    // Merge local counts into global counts (synchronized)
                    synchronized (coOccurrenceCount) {
                        for (int x = 0; x < taxaCount; x++) {
                            for (int y = 0; y < taxaCount; y++) {
                                coOccurrenceCount[x][y] += localCoOccurrenceCount[x][y];
                            }
                        }
                    }
                    
                    if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                        System.out.println("Thread " + threadId + " completed processing " + (endIdx - startIdx) + " gene trees");
                    }
                    
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Distance matrix calculation was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        // Convert co-occurrence counts to distances
        // Higher co-occurrence = lower distance
        for (int i = 0; i < taxaCount; i++) {
            for (int j = i + 1; j < taxaCount; j++) {
                double coOccurrenceRate = (double) coOccurrenceCount[i][j] / geneTrees.size();
                // Distance = 1 - co-occurrence rate (so high co-occurrence = low distance)
                double distance = Math.max(0.01, 1.0 - coOccurrenceRate); // Minimum distance to avoid zero
                distanceMatrix[i][j] = distance;
                distanceMatrix[j][i] = distance;
            }
        }
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("Parallel distance matrix calculation completed.");
        }
    }
    
    /**
     * Collect all clades (subtrees) for each taxon in a gene tree
     */
    private void collectClades(TreeNode node, Map<Integer, BitSet> taxonToClade) {
        if (node.isLeaf()) {
            BitSet clade = new BitSet(taxaCount);
            clade.set(node.taxon.id);
            taxonToClade.put(node.taxon.id, clade);
        } else {
            BitSet nodeClade = new BitSet(taxaCount);
            for (TreeNode child : node.childs) {
                collectClades(child, taxonToClade);
                BitSet childClade = getNodeClade(child);
                nodeClade.or(childClade);
            }
            
            // Update taxon-to-clade mapping for all taxa in this subtree
            for (int i = nodeClade.nextSetBit(0); i >= 0; i = nodeClade.nextSetBit(i + 1)) {
                taxonToClade.put(i, (BitSet) nodeClade.clone());
            }
        }
    }
    
    /**
     * Get the clade (set of taxa) represented by a tree node
     */
    private BitSet getNodeClade(TreeNode node) {
        BitSet clade = new BitSet(taxaCount);
        collectTaxaInSubtree(node, clade);
        return clade;
    }
    
    /**
     * Collect all taxa in a subtree
     */
    private void collectTaxaInSubtree(TreeNode node, BitSet clade) {
        if (node.isLeaf()) {
            clade.set(node.taxon.id);
        } else {
            for (TreeNode child : node.childs) {
                collectTaxaInSubtree(child, clade);
            }
        }
    }
    
    /**
     * Check if two taxa are in the same clade in a gene tree
     */
    private boolean isInSameClade(int taxon1, int taxon2, TreeNode root) {
        TreeNode lca = findLCA(root, taxon1, taxon2);
        if (lca == null) return false;
        
        BitSet lcaClade = getNodeClade(lca);
        return lcaClade.get(taxon1) && lcaClade.get(taxon2);
    }
    
    /**
     * Find the Lowest Common Ancestor of two taxa
     */
    private TreeNode findLCA(TreeNode root, int taxon1, int taxon2) {
        if (root == null) return null;
        
        if (root.isLeaf()) {
            return (root.taxon.id == taxon1 || root.taxon.id == taxon2) ? root : null;
        }
        
        List<TreeNode> childrenWithTargets = new ArrayList<>();
        for (TreeNode child : root.childs) {
            TreeNode result = findLCA(child, taxon1, taxon2);
            if (result != null) {
                childrenWithTargets.add(result);
            }
        }
        
        if (childrenWithTargets.size() >= 2) {
            return root; // This is the LCA
        } else if (childrenWithTargets.size() == 1) {
            return childrenWithTargets.get(0);
        } else {
            return null;
        }
    }
    
    /**
     * Build UPGMA tree from distance matrix
     */
    public Tree buildUPGMATree() {
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("Building UPGMA tree...");
        }
        
        // Create initial clusters (one per taxon)
        List<UPGMACluster> clusters = new ArrayList<>();
        for (int i = 0; i < taxaCount; i++) {
            clusters.add(new UPGMACluster(i, taxaNames[i]));
        }
        
        // UPGMA clustering algorithm
        while (clusters.size() > 1) {
            // Find closest pair of clusters
            double minDistance = Double.MAX_VALUE;
            int cluster1Index = -1, cluster2Index = -1;
            
            for (int i = 0; i < clusters.size(); i++) {
                for (int j = i + 1; j < clusters.size(); j++) {
                    double distance = calculateClusterDistance(clusters.get(i), clusters.get(j));
                    if (distance < minDistance) {
                        minDistance = distance;
                        cluster1Index = i;
                        cluster2Index = j;
                    }
                }
            }
            
            // Merge closest clusters
            UPGMACluster cluster1 = clusters.get(cluster1Index);
            UPGMACluster cluster2 = clusters.get(cluster2Index);
            UPGMACluster mergedCluster = new UPGMACluster(cluster1, cluster2, minDistance / 2.0);
            
            // Remove old clusters and add merged cluster
            clusters.remove(Math.max(cluster1Index, cluster2Index)); // Remove higher index first
            clusters.remove(Math.min(cluster1Index, cluster2Index));
            clusters.add(mergedCluster);
        }
        
        // Convert UPGMA cluster to Tree
        Tree tree = new Tree();
        tree.taxaMap = geneTrees.taxaMap;
        tree.root = convertUPGMAClusterToTreeNode(clusters.get(0), tree);
        tree.isRooted = true;
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("UPGMA tree construction completed.");
        }
        
        return tree;
    }
    
    /**
     * Calculate distance between two UPGMA clusters
     */
    private double calculateClusterDistance(UPGMACluster cluster1, UPGMACluster cluster2) {
        double totalDistance = 0.0;
        int count = 0;
        
        for (int taxon1 : cluster1.taxa) {
            for (int taxon2 : cluster2.taxa) {
                totalDistance += distanceMatrix[taxon1][taxon2];
                count++;
            }
        }
        
        return count > 0 ? totalDistance / count : Double.MAX_VALUE;
    }
    
    /**
     * Convert UPGMA cluster to TreeNode
     */
    private TreeNode convertUPGMAClusterToTreeNode(UPGMACluster cluster, Tree tree) {
        if (cluster.isLeaf()) {
            return tree.addLeaf(geneTrees.taxa[cluster.taxa.get(0)]);
        } else {
            ArrayList<TreeNode> children = new ArrayList<>();
            children.add(convertUPGMAClusterToTreeNode(cluster.leftChild, tree));
            children.add(convertUPGMAClusterToTreeNode(cluster.rightChild, tree));
            return tree.addInternalNode(children);
        }
    }
    
    /**
     * Extract bipartitions from a distance-based tree
     */
    public List<STBipartition> extractBipartitions(Tree tree) {
        List<STBipartition> bipartitions = new ArrayList<>();
        if (tree != null && tree.root != null) {
            extractBipartitionsRecursive(tree.root, bipartitions);
        }
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("  Distance tree nodes processed, extracted " + bipartitions.size() + " bipartitions");
            if (bipartitions.size() > 0) {
                System.out.println("  Sample distance bipartitions:");
                for (int i = 0; i < Math.min(3, bipartitions.size()); i++) {
                    STBipartition bip = bipartitions.get(i);
                    System.out.println("    " + bip.cluster1.cardinality() + " vs " + bip.cluster2.cardinality() + " taxa");
                }
            }
        }
        
        return bipartitions;
    }
    
    /**
     * Recursively extract bipartitions from tree
     */
    private BitSet extractBipartitionsRecursive(TreeNode node, List<STBipartition> bipartitions) {
        if (node.isLeaf()) {
            BitSet leafSet = new BitSet(taxaCount);
            leafSet.set(node.taxon.id);
            return leafSet;
        }
        
        if (node.childs.size() == 2) {
            BitSet leftCluster = extractBipartitionsRecursive(node.childs.get(0), bipartitions);
            BitSet rightCluster = extractBipartitionsRecursive(node.childs.get(1), bipartitions);
            
            // Add bipartition if both clusters are non-empty
            if (leftCluster.cardinality() > 0 && rightCluster.cardinality() > 0) {
                STBipartition bipartition = new STBipartition(leftCluster, rightCluster);
                bipartitions.add(bipartition);
            }
            
            BitSet nodeCluster = (BitSet) leftCluster.clone();
            nodeCluster.or(rightCluster);
            return nodeCluster;
        } else {
            // Handle polytomies by combining all children
            BitSet nodeCluster = new BitSet(taxaCount);
            for (TreeNode child : node.childs) {
                BitSet childCluster = extractBipartitionsRecursive(child, bipartitions);
                nodeCluster.or(childCluster);
            }
            return nodeCluster;
        }
    }
    
    /**
     * Get the distance matrix
     */
    public double[][] getDistanceMatrix() {
        return distanceMatrix;
    }
    
    /**
     * Inner class representing a cluster in UPGMA algorithm
     */
    private static class UPGMACluster {
        List<Integer> taxa;
        String name;
        UPGMACluster leftChild;
        UPGMACluster rightChild;
        double height;
        
        // Leaf cluster constructor
        public UPGMACluster(int taxonId, String taxonName) {
            this.taxa = Arrays.asList(taxonId);
            this.name = taxonName;
            this.height = 0.0;
        }
        
        // Internal cluster constructor
        public UPGMACluster(UPGMACluster left, UPGMACluster right, double height) {
            this.taxa = new ArrayList<>();
            this.taxa.addAll(left.taxa);
            this.taxa.addAll(right.taxa);
            this.name = "(" + left.name + "," + right.name + ")";
            this.leftChild = left;
            this.rightChild = right;
            this.height = height;
        }
        
        public boolean isLeaf() {
            return leftChild == null && rightChild == null;
        }
    }
}
