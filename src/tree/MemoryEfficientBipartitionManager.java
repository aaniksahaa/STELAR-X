package tree;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import utils.BitSet;
import utils.Threading;

/**
 * Memory-efficient bipartition manager that reduces memory usage from O(n²k) to O(nk).
 * 
 * Instead of immediately creating BitSet representations for all bipartitions,
 * this manager:
 * 1. Represents bipartitions as ranges during initial processing
 * 2. Uses hash-based equality checking to filter unique bipartitions
 * 3. Only converts unique bipartitions to BitSets at the end
 * 
 * This optimization is particularly effective when there are many duplicate
 * bipartitions across gene trees.
 */
public class MemoryEfficientBipartitionManager {
    
    // Configuration
    public static RangeBipartition.HashFunction DEFAULT_HASH_FUNCTION = new RangeBipartition.SumHashFunction();
    public static boolean ENABLE_DOUBLE_HASHING = false;
    
    // Gene tree data structures
    private final List<Tree> geneTrees;
    private final int[][] geneTreeTaxaOrdering;  // [tree_index][taxa_position] = taxon_id
    private final long[][] prefixSums;           // [tree_index][position] = prefix sum of taxon IDs
    
    // Processing results
    private final Map<Long, List<RangeBipartition>> hashToBipartitions;
    private final Map<STBipartition, Integer> uniqueSTBipartitions;
    private final int realTaxaCount;
    
    public MemoryEfficientBipartitionManager(List<Tree> geneTrees, int realTaxaCount) {
        this.geneTrees = geneTrees;
        this.realTaxaCount = realTaxaCount;
        this.geneTreeTaxaOrdering = new int[geneTrees.size()][];
        this.prefixSums = new long[geneTrees.size()][];
        this.hashToBipartitions = new ConcurrentHashMap<>();
        this.uniqueSTBipartitions = new ConcurrentHashMap<>();
        
        initializeGeneTreeOrderings();
        calculatePrefixSums();
    }
    
    /**
     * Initialize the left-to-right ordering for each gene tree (inorder traversal of leaves).
     */
    private void initializeGeneTreeOrderings() {
        System.out.println("Initializing gene tree taxa orderings...");
        
        for (int i = 0; i < geneTrees.size(); i++) {
            Tree tree = geneTrees.get(i);
            List<Integer> ordering = new ArrayList<>();
            collectLeavesInOrder(tree.root, ordering);
            
            geneTreeTaxaOrdering[i] = ordering.stream().mapToInt(Integer::intValue).toArray();
        }
        
        System.out.println("Initialized orderings for " + geneTrees.size() + " gene trees");
    }
    
    /**
     * Collect leaves in left-to-right order (inorder traversal).
     */
    private void collectLeavesInOrder(TreeNode node, List<Integer> ordering) {
        if (node.isLeaf()) {
            ordering.add(node.taxon.id);
            return;
        }
        
        // For binary trees, visit left child, then right child
        if (node.childs != null && node.childs.size() >= 2) {
            collectLeavesInOrder(node.childs.get(0), ordering);
            collectLeavesInOrder(node.childs.get(1), ordering);
        }
        
        // Handle any additional children (though binary trees should only have 2)
        for (int i = 2; i < (node.childs != null ? node.childs.size() : 0); i++) {
            collectLeavesInOrder(node.childs.get(i), ordering);
        }
    }
    
    /**
     * Calculate prefix sums for efficient range hash computation.
     */
    private void calculatePrefixSums() {
        System.out.println("Calculating prefix sums for hash functions...");
        
        for (int i = 0; i < geneTrees.size(); i++) {
            int[] ordering = geneTreeTaxaOrdering[i];
            prefixSums[i] = new long[ordering.length];
            
            if (ordering.length > 0) {
                prefixSums[i][0] = ordering[0];
                for (int j = 1; j < ordering.length; j++) {
                    prefixSums[i][j] = prefixSums[i][j - 1] + ordering[j];
                }
            }
        }
        
        System.out.println("Calculated prefix sums for efficient range hashing");
    }
    
    /**
     * Process all gene trees to extract range bipartitions in parallel.
     */
    public Map<STBipartition, Integer> processGeneTreesParallel() {
        System.out.println("Processing gene trees with memory-efficient range bipartitions...");
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        CountDownLatch latch = new CountDownLatch(numThreads);
        AtomicInteger processedTrees = new AtomicInteger(0);
        
        int chunkSize = (geneTrees.size() + numThreads - 1) / numThreads;
        
        System.out.println("Using " + numThreads + " threads for parallel processing");
        
        // Process gene trees in parallel chunks
        for (int i = 0; i < numThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, geneTrees.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    Map<Long, List<RangeBipartition>> localHashMap = new HashMap<>();
                    
                    System.out.println("Thread " + threadId + " processing trees " + startIdx + " to " + (endIdx - 1));
                    
                    for (int treeIdx = startIdx; treeIdx < endIdx; treeIdx++) {
                        Tree tree = geneTrees.get(treeIdx);
                        if (tree.isRooted()) {
                            extractRangeBipartitions(tree.root, treeIdx, localHashMap);
                        }
                        processedTrees.incrementAndGet();
                    }
                    
                    // Merge local results into global map
                    synchronized (hashToBipartitions) {
                        for (Map.Entry<Long, List<RangeBipartition>> entry : localHashMap.entrySet()) {
                            hashToBipartitions.computeIfAbsent(entry.getKey(), k -> new ArrayList<>())
                                             .addAll(entry.getValue());
                        }
                    }
                    
                    System.out.println("Thread " + threadId + " completed processing " + (endIdx - startIdx) + " trees");
                    
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Gene tree processing was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        System.out.println("Processed " + processedTrees.get() + " gene trees");
        System.out.println("Found " + hashToBipartitions.size() + " unique hash groups");
        
        // Convert range bipartitions to STBipartitions
        convertToSTBipartitions();
        
        return uniqueSTBipartitions;
    }
    
    /**
     * Extract range bipartitions from a gene tree recursively.
     * This mirrors the original calculateSTBipartitionsUtilLocal logic.
     */
    private void extractRangeBipartitions(TreeNode node, int treeIndex, Map<Long, List<RangeBipartition>> localHashMap) {
        if (node.isLeaf()) {
            return;
        }
        
        if (node.childs != null && node.childs.size() == 2) {
            // Calculate ranges for both left and right subtrees
            int[] leftRange = calculateSubtreeRange(node.childs.get(0), treeIndex);
            int[] rightRange = calculateSubtreeRange(node.childs.get(1), treeIndex);
            
            if (leftRange != null && rightRange != null && 
                leftRange[0] <= leftRange[1] && rightRange[0] <= rightRange[1]) {
                
                // Create bipartition with both left and right ranges
                RangeBipartition rangeBip = new RangeBipartition(treeIndex, 
                    leftRange[0], leftRange[1] + 1,     // left range (exclusive end)
                    rightRange[0], rightRange[1] + 1);  // right range (exclusive end)
                
                long hash = DEFAULT_HASH_FUNCTION.calculateHash(rangeBip, prefixSums);
                
                localHashMap.computeIfAbsent(hash, k -> new ArrayList<>()).add(rangeBip);
            }
        }
        
        // Recursively process children
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                extractRangeBipartitions(child, treeIndex, localHashMap);
            }
        }
    }
    
    /**
     * Calculate the range (min and max positions) covered by a subtree.
     * Returns [min_position, max_position] or null if no valid range found.
     */
    private int[] calculateSubtreeRange(TreeNode node, int treeIndex) {
        if (node.isLeaf()) {
            // Find position of this taxon in the ordering
            int[] ordering = geneTreeTaxaOrdering[treeIndex];
            for (int i = 0; i < ordering.length; i++) {
                if (ordering[i] == node.taxon.id) {
                    return new int[]{i, i};
                }
            }
            return null; // Taxon not found
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
        
        if (minPos == Integer.MAX_VALUE || maxPos == Integer.MIN_VALUE) {
            return null;
        }
        
        return new int[]{minPos, maxPos};
    }
    
    /**
     * Convert range bipartitions to STBipartitions only for unique ones.
     */
    private void convertToSTBipartitions() {
        System.out.println("Converting unique range bipartitions to STBipartitions...");
        
        int totalRanges = 0;
        int uniqueRanges = 0;
        
        for (Map.Entry<Long, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
            List<RangeBipartition> ranges = entry.getValue();
            totalRanges += ranges.size();
            
            if (ranges.isEmpty()) continue;
            
            // Group ranges by actual equality (not just hash)
            Map<RangeBipartition, Integer> actuallyUniqueRanges = new HashMap<>();
            
            for (RangeBipartition range : ranges) {
                boolean found = false;
                
                for (RangeBipartition existing : actuallyUniqueRanges.keySet()) {
                    if (rangesAreEqual(range, existing)) {
                        actuallyUniqueRanges.put(existing, actuallyUniqueRanges.get(existing) + 1);
                        found = true;
                        break;
                    }
                }
                
                if (!found) {
                    actuallyUniqueRanges.put(range, 1);
                }
            }
            
            // Convert unique ranges to STBipartitions
            for (Map.Entry<RangeBipartition, Integer> rangeEntry : actuallyUniqueRanges.entrySet()) {
                STBipartition stb = convertRangeToSTBipartition(rangeEntry.getKey());
                if (stb != null) {
                    uniqueSTBipartitions.merge(stb, rangeEntry.getValue(), Integer::sum);
                    uniqueRanges++;
                }
            }
        }
        
        System.out.println("Memory optimization results:");
        System.out.println("  Total range bipartitions processed: " + totalRanges);
        System.out.println("  Unique STBipartitions created: " + uniqueRanges);
        System.out.println("  Memory reduction factor: " + ((double) totalRanges / Math.max(1, uniqueRanges)));
        System.out.println("  Final unique STBipartitions: " + uniqueSTBipartitions.size());
    }
    
    /**
     * Check if two range bipartitions are equal by comparing their actual taxon sets.
     */
    private boolean rangesAreEqual(RangeBipartition range1, RangeBipartition range2) {
        // Same tree - compare ranges directly
        if (range1.geneTreeIndex == range2.geneTreeIndex) {
            return range1.leftStart == range2.leftStart && range1.leftEnd == range2.leftEnd &&
                   range1.rightStart == range2.rightStart && range1.rightEnd == range2.rightEnd;
        }
        
        // Different trees - compare hash values for efficiency
        long hash1 = DEFAULT_HASH_FUNCTION.calculateHash(range1, prefixSums);
        long hash2 = DEFAULT_HASH_FUNCTION.calculateHash(range2, prefixSums);
        
        if (hash1 != hash2) {
            return false;
        }
        
        // Hash collision possible - compare actual taxon sets
        return getLeftTaxonSetForRange(range1).equals(getLeftTaxonSetForRange(range2)) &&
               getRightTaxonSetForRange(range1).equals(getRightTaxonSetForRange(range2));
    }
    
    /**
     * Get the set of taxon IDs for the left cluster of a range bipartition.
     */
    private Set<Integer> getLeftTaxonSetForRange(RangeBipartition range) {
        Set<Integer> taxonSet = new HashSet<>();
        int[] ordering = geneTreeTaxaOrdering[range.geneTreeIndex];
        
        for (int i = range.leftStart; i < range.leftEnd && i < ordering.length; i++) {
            taxonSet.add(ordering[i]);
        }
        
        return taxonSet;
    }
    
    /**
     * Get the set of taxon IDs for the right cluster of a range bipartition.
     */
    private Set<Integer> getRightTaxonSetForRange(RangeBipartition range) {
        Set<Integer> taxonSet = new HashSet<>();
        int[] ordering = geneTreeTaxaOrdering[range.geneTreeIndex];
        
        for (int i = range.rightStart; i < range.rightEnd && i < ordering.length; i++) {
            taxonSet.add(ordering[i]);
        }
        
        return taxonSet;
    }
    
    /**
     * Convert a range bipartition to an STBipartition.
     */
    private STBipartition convertRangeToSTBipartition(RangeBipartition range) {
        Set<Integer> leftSet = getLeftTaxonSetForRange(range);
        Set<Integer> rightSet = getRightTaxonSetForRange(range);
        
        if (leftSet.isEmpty() || rightSet.isEmpty()) {
            return null; // Invalid bipartition
        }
        
        // Verify that left and right sets are disjoint and cover all taxa in this subtree
        Set<Integer> union = new HashSet<>(leftSet);
        union.addAll(rightSet);
        Set<Integer> intersection = new HashSet<>(leftSet);
        intersection.retainAll(rightSet);
        
        if (!intersection.isEmpty()) {
            return null; // Overlapping clusters - invalid bipartition
        }
        
        // Create BitSets
        BitSet cluster1 = new BitSet(realTaxaCount);
        BitSet cluster2 = new BitSet(realTaxaCount);
        
        // Set cluster1 (left cluster)
        for (int taxonId : leftSet) {
            cluster1.set(taxonId);
        }
        
        // Set cluster2 (right cluster)  
        for (int taxonId : rightSet) {
            cluster2.set(taxonId);
        }
        
        if (cluster1.cardinality() == 0 || cluster2.cardinality() == 0) {
            return null; // Invalid bipartition
        }
        
        return new STBipartition(cluster1, cluster2);
    }
    
    /**
     * Get processing statistics.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Memory-Efficient Bipartition Processing Statistics:\n");
        sb.append("  Gene trees processed: ").append(geneTrees.size()).append("\n");
        sb.append("  Hash function: ").append(DEFAULT_HASH_FUNCTION.getName()).append("\n");
        sb.append("  Unique hash groups: ").append(hashToBipartitions.size()).append("\n");
        sb.append("  Final unique STBipartitions: ").append(uniqueSTBipartitions.size()).append("\n");
        
        int totalRanges = hashToBipartitions.values().stream()
                         .mapToInt(List::size)
                         .sum();
        
        sb.append("  Total range bipartitions: ").append(totalRanges).append("\n");
        sb.append("  Memory reduction factor: ").append(String.format("%.2f", 
            (double) totalRanges / Math.max(1, uniqueSTBipartitions.size()))).append("x\n");
        
        return sb.toString();
    }
}
