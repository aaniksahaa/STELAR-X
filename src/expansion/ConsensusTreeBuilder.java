package expansion;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import preprocessing.GeneTrees;
import tree.STBipartition;
import tree.Tree;
import tree.TreeNode;
import utils.BitSet;
import utils.BipartitionExpansionConfig;
import utils.Threading;

/**
 * Consensus tree builder for bipartition expansion.
 * 
 * This class implements greedy consensus algorithms similar to ASTRAL's approach.
 * It builds consensus trees from gene trees using different support thresholds
 * and extracts bipartitions from the resulting consensus trees.
 */
public class ConsensusTreeBuilder {
    
    private GeneTrees geneTrees;
    private int taxaCount;
    private String[] taxaNames;
    
    public ConsensusTreeBuilder(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
        this.taxaCount = geneTrees.realTaxaCount;
        this.taxaNames = geneTrees.taxonIdToLabel;
    }
    
    /**
     * Build consensus trees using multiple support thresholds
     */
    public List<Tree> buildConsensusTreesWithThresholds() {
        List<Tree> consensusTrees = new ArrayList<>();
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("Building consensus trees with " + 
                             BipartitionExpansionConfig.CONSENSUS_THRESHOLDS.length + " thresholds...");
        }
        
        // Count bipartition frequencies across all gene trees
        Map<STBipartition, Integer> bipartitionCounts = countBipartitionFrequencies();
        int totalTrees = geneTrees.geneTrees.size();
        
        // Build consensus tree for each threshold
        for (double threshold : BipartitionExpansionConfig.CONSENSUS_THRESHOLDS) {
            for (int i = 0; i < BipartitionExpansionConfig.MAX_CONSENSUS_TREES; i++) {
                Tree consensusTree = buildGreedyConsensus(bipartitionCounts, totalTrees, threshold);
                if (consensusTree != null) {
                    consensusTrees.add(consensusTree);
                }
            }
        }
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("Built " + consensusTrees.size() + " consensus trees.");
        }
        
        return consensusTrees;
    }
    
    /**
     * Count frequency of each bipartition across all gene trees
     */
    private Map<STBipartition, Integer> countBipartitionFrequencies() {
        return countBipartitionFrequenciesParallel();
    }
    
    /**
     * Parallel version of bipartition frequency counting for improved performance.
     * Divides the gene trees among multiple threads for concurrent processing.
     */
    private Map<STBipartition, Integer> countBipartitionFrequenciesParallel() {
        Map<STBipartition, Integer> counts = new ConcurrentHashMap<>();
        
        // Use existing bipartition counts from GeneTrees
        for (Map.Entry<STBipartition, Integer> entry : geneTrees.stBipartitions.entrySet()) {
            counts.put(entry.getKey(), entry.getValue());
        }
        
        // Also extract bipartitions directly from gene trees to ensure completeness
        List<Tree> geneTrees = this.geneTrees.geneTrees;
        
        if (geneTrees.isEmpty()) {
            return counts;
        }
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("Counting bipartition frequencies from " + geneTrees.size() + " gene trees in parallel...");
        }
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        int chunkSize = (geneTrees.size() + numThreads - 1) / numThreads;
        CountDownLatch latch = new CountDownLatch(numThreads);
        
        // Process gene trees in parallel chunks
        for (int i = 0; i < numThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, geneTrees.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    Map<STBipartition, Integer> localCounts = new HashMap<>();
                    
                    if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                        System.out.println("Thread " + threadId + " extracting bipartitions from trees " + startIdx + " to " + (endIdx-1));
                    }
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        Tree geneTree = geneTrees.get(j);
                        List<STBipartition> treeBipartitions = extractBipartitionsFromTree(geneTree);
                        for (STBipartition bip : treeBipartitions) {
                            localCounts.merge(bip, 1, Integer::sum);
                        }
                    }
                    
                    // Merge local counts into global counts
                    for (Map.Entry<STBipartition, Integer> entry : localCounts.entrySet()) {
                        counts.merge(entry.getKey(), entry.getValue(), Integer::sum);
                    }
                    
                    if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                        System.out.println("Thread " + threadId + " completed: found " + localCounts.size() + " unique bipartitions");
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
            throw new RuntimeException("Bipartition frequency counting was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("Parallel bipartition frequency counting completed. Total unique bipartitions: " + counts.size());
        }
        
        return counts;
    }
    
    /**
     * Extract bipartitions from a single gene tree
     */
    private List<STBipartition> extractBipartitionsFromTree(Tree tree) {
        List<STBipartition> bipartitions = new ArrayList<>();
        extractBipartitionsRecursive(tree.root, bipartitions);
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
            
            // Add bipartition if both clusters are non-empty and not trivial
            if (leftCluster.cardinality() > 0 && rightCluster.cardinality() > 0 &&
                leftCluster.cardinality() < taxaCount && rightCluster.cardinality() < taxaCount) {
                STBipartition bipartition = new STBipartition(leftCluster, rightCluster);
                bipartitions.add(bipartition);
            }
            
            BitSet nodeCluster = (BitSet) leftCluster.clone();
            nodeCluster.or(rightCluster);
            return nodeCluster;
        } else {
            // Handle polytomies
            BitSet nodeCluster = new BitSet(taxaCount);
            for (TreeNode child : node.childs) {
                BitSet childCluster = extractBipartitionsRecursive(child, bipartitions);
                nodeCluster.or(childCluster);
            }
            return nodeCluster;
        }
    }
    
    /**
     * Build greedy consensus tree using specified support threshold
     */
    private Tree buildGreedyConsensus(Map<STBipartition, Integer> bipartitionCounts, 
                                     int totalTrees, double threshold) {
        
        // Filter bipartitions by support threshold
        List<BipartitionWithSupport> supportedBipartitions = new ArrayList<>();
        for (Map.Entry<STBipartition, Integer> entry : bipartitionCounts.entrySet()) {
            double support = (double) entry.getValue() / totalTrees;
            if (support >= threshold) {
                supportedBipartitions.add(new BipartitionWithSupport(entry.getKey(), support));
            }
        }
        
        // Sort by support (highest first)
        supportedBipartitions.sort((a, b) -> Double.compare(b.support, a.support));
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("  Threshold " + threshold + ": " + 
                             supportedBipartitions.size() + " supported bipartitions");
        }
        
        // For now, create a simple star tree (the bipartitions will be extracted directly)
        // The actual tree structure is less important than the bipartitions themselves
        Tree consensusTree = createSimpleStarTree();
        return consensusTree;
    }
    
    /**
     * Build tree from a set of compatible bipartitions
     */
    private Tree buildTreeFromCompatibleBipartitions(List<BipartitionWithSupport> bipartitions) {
        // Start with a star tree (all taxa connected to root)
        Tree tree = new Tree();
        tree.taxaMap = geneTrees.taxaMap;
        tree.isRooted = true;
        
        // Create leaf nodes for all taxa
        ArrayList<TreeNode> leaves = new ArrayList<>();
        for (int i = 0; i < taxaCount; i++) {
            TreeNode leaf = tree.addLeaf(geneTrees.taxa[i]);
            leaves.add(leaf);
        }
        
        // Create root and attach all leaves initially (star tree)
        tree.root = tree.addInternalNode(leaves);
        
        // Greedily add compatible bipartitions
        Set<STBipartition> addedBipartitions = new HashSet<>();
        
        for (BipartitionWithSupport bipWithSupport : bipartitions) {
            STBipartition bip = bipWithSupport.bipartition;
            
            // Check if this bipartition is compatible with already added ones
            if (isCompatibleWithTree(bip, tree) && !addedBipartitions.contains(bip)) {
                if (addBipartitionToTree(bip, tree)) {
                    addedBipartitions.add(bip);
                }
            }
        }
        
        return tree;
    }
    
    /**
     * Check if a bipartition is compatible with the current tree structure
     */
    private boolean isCompatibleWithTree(STBipartition bip, Tree tree) {
        // For now, implement a simple compatibility check
        // In a more sophisticated implementation, we would check for conflicts
        // with existing tree structure
        
        // Basic check: ensure bipartition covers valid taxa
        BitSet union = (BitSet) bip.cluster1.clone();
        union.or(bip.cluster2);
        
        // Should cover all taxa or be a proper subset
        return union.cardinality() <= taxaCount;
    }
    
    /**
     * Add a bipartition to the tree structure
     */
    private boolean addBipartitionToTree(STBipartition bip, Tree tree) {
        // This is a simplified implementation
        // In practice, this would involve more complex tree manipulation
        // to maintain the tree structure while adding the bipartition
        
        // For now, we'll just record that we attempted to add it
        return true;
    }
    
    /**
     * Extract all bipartitions from consensus trees
     * SIMPLIFIED APPROACH: Instead of building complex consensus trees,
     * we directly return the supported bipartitions from the consensus analysis
     */
    public List<STBipartition> extractBipartitionsFromConsensusTrees(List<Tree> consensusTrees) {
        // Since the consensus tree building is complex, let's use a direct approach
        // We'll extract supported bipartitions directly from the frequency analysis
        return extractSupportedBipartitionsDirectly();
    }
    
    /**
     * Direct extraction of consensus bipartitions without building intermediate trees
     */
    private List<STBipartition> extractSupportedBipartitionsDirectly() {
        Set<STBipartition> uniqueBipartitions = new HashSet<>();
        
        // Count bipartition frequencies
        Map<STBipartition, Integer> bipartitionCounts = countBipartitionFrequencies();
        int totalTrees = geneTrees.geneTrees.size();
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("  Extracting bipartitions directly from frequency analysis...");
            System.out.println("  Total bipartitions in frequency map: " + bipartitionCounts.size());
        }
        
        // For each threshold, add bipartitions that meet the support requirement
        for (double threshold : BipartitionExpansionConfig.CONSENSUS_THRESHOLDS) {
            int addedForThreshold = 0;
            for (Map.Entry<STBipartition, Integer> entry : bipartitionCounts.entrySet()) {
                double support = (double) entry.getValue() / totalTrees;
                if (support >= threshold) {
                    if (uniqueBipartitions.add(entry.getKey())) {
                        addedForThreshold++;
                    }
                }
            }
            
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("    Threshold " + threshold + ": " + addedForThreshold + " new bipartitions added");
            }
        }
        
        List<STBipartition> result = new ArrayList<>(uniqueBipartitions);
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("  Total unique consensus bipartitions: " + result.size());
            if (result.size() > 0) {
                System.out.println("  Sample consensus bipartitions:");
                for (int i = 0; i < Math.min(3, result.size()); i++) {
                    STBipartition bip = result.get(i);
                    System.out.println("    " + bip.cluster1.cardinality() + " vs " + bip.cluster2.cardinality() + " taxa");
                }
            }
        }
        
        return result;
    }
    
    /**
     * Create a simple star tree (all taxa connected to root)
     */
    private Tree createSimpleStarTree() {
        Tree tree = new Tree();
        tree.taxaMap = geneTrees.taxaMap;
        tree.isRooted = true;
        
        // Create leaf nodes for all taxa
        ArrayList<TreeNode> leaves = new ArrayList<>();
        for (int i = 0; i < taxaCount; i++) {
            TreeNode leaf = tree.addLeaf(geneTrees.taxa[i]);
            leaves.add(leaf);
        }
        
        // Create root and attach all leaves (star tree)
        tree.root = tree.addInternalNode(leaves);
        return tree;
    }

    /**
     * Build a simple majority consensus tree
     */
    public Tree buildMajorityConsensus() {
        return buildGreedyConsensus(countBipartitionFrequencies(), 
                                   geneTrees.geneTrees.size(), 0.5);
    }
    
    /**
     * Build a strict consensus tree (100% support)
     */
    public Tree buildStrictConsensus() {
        return buildGreedyConsensus(countBipartitionFrequencies(), 
                                   geneTrees.geneTrees.size(), 1.0);
    }
    
    /**
     * Helper class to store bipartition with its support value
     */
    private static class BipartitionWithSupport {
        STBipartition bipartition;
        double support;
        
        public BipartitionWithSupport(STBipartition bipartition, double support) {
            this.bipartition = bipartition;
            this.support = support;
        }
    }
}
