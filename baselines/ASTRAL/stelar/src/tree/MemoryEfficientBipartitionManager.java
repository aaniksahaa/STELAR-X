package tree;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

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
    public static boolean ENABLE_DOUBLE_HASHING = true;
    public static boolean ENABLE_EXPENSIVE_EQUALITY_CHECKS = false; // Default: trust the hash functions
    public static RangeBipartition.HashFunction DEFAULT_HASH_FUNCTION = getDefaultHashFunction();
    
    private static RangeBipartition.HashFunction getDefaultHashFunction() {
        return ENABLE_DOUBLE_HASHING ? new RangeBipartition.DoubleHashFunction() : new RangeBipartition.SumHashFunction();
    }
    
    // Gene tree data structures
    private final List<Tree> geneTrees;
    private final int[][] geneTreeTaxaOrdering;  // [tree_index][taxa_position] = taxon_id
    private final long[][] prefixSums;           // [tree_index][position] = prefix sum of taxon IDs
    private final long[][] prefixXORs;           // [tree_index][position] = prefix XOR of taxon IDs
    
    // Processing results
    private final Map<Object, List<RangeBipartition>> hashToBipartitions;
    private final Map<RangeBipartition, Integer> uniqueRangeBipartitions;
    private final int realTaxaCount;
    
    public MemoryEfficientBipartitionManager(List<Tree> geneTrees, int realTaxaCount) {
        this.geneTrees = geneTrees;
        this.realTaxaCount = realTaxaCount;
        this.geneTreeTaxaOrdering = new int[geneTrees.size()][];
        this.prefixSums = new long[geneTrees.size()][];
        this.prefixXORs = new long[geneTrees.size()][];
        this.hashToBipartitions = new ConcurrentHashMap<>();
        this.uniqueRangeBipartitions = new ConcurrentHashMap<>();
        
        initializeGeneTreeOrderings();
        calculatePrefixArrays();
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
     * Calculate prefix sums and XORs over hashed taxon IDs for efficient range hash computation.
     * This provides much better hash distribution than using raw taxon IDs.
     */
    private void calculatePrefixArrays() {
        System.out.println("Calculating prefix sums and XORs over hashed taxon IDs for hash functions...");
        
        for (int i = 0; i < geneTrees.size(); i++) {
            int[] ordering = geneTreeTaxaOrdering[i];
            prefixSums[i] = new long[ordering.length];
            prefixXORs[i] = new long[ordering.length];
            
            if (ordering.length > 0) {
                // Hash the first taxon ID and store as initial values
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
        
        System.out.println("Calculated prefix sums and XORs over hashed taxon IDs for efficient range hashing");
    }
    
    /**
     * Simple but effective hash function for individual taxon IDs.
     * Uses a combination of multiplications and XOR operations for good distribution.
     */
    private static long hashSingleTaxon(int taxonId) {
        // Use a strong hash mixing function for individual taxon IDs
        long x = taxonId;
        x ^= x >>> 16;
        x *= 0x85ebca6b;
        x ^= x >>> 13;
        x *= 0xc2b2ae35;
        x ^= x >>> 16;
        
        // Ensure we never return 0 to avoid issues with prefix operations
        return x == 0 ? 1 : x;
    }
    
    /**
     * Process all gene trees to extract range bipartitions in parallel.
     */
    public Map<RangeBipartition, Integer> processGeneTreesParallel() {
        System.out.println("Processing gene trees with memory-efficient range bipartitions...");
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        // Calculate optimal number of threads to avoid invalid ranges
        int chunkSize = Math.max(1, (geneTrees.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (geneTrees.size() + chunkSize - 1) / chunkSize);
        
        CountDownLatch latch = new CountDownLatch(actualThreads);
        AtomicInteger processedTrees = new AtomicInteger(0);
        
        System.out.println("Using " + actualThreads + " threads for parallel processing");
        
        // Process gene trees in parallel chunks
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, geneTrees.size());
            final int threadId = i;
            
            // Skip threads that would have invalid ranges
            if (startIdx >= geneTrees.size()) {
                continue;
            }
            
            Threading.execute(() -> {
                try {
                    Map<Object, List<RangeBipartition>> localHashMap = new HashMap<>();
                    
                    // Validate range before processing
                    if (startIdx >= endIdx || startIdx >= geneTrees.size()) {
                        System.out.println("Thread " + threadId + " skipped - invalid range [" + startIdx + ", " + endIdx + ")");
                        return;
                    }
                    
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
                        for (Map.Entry<Object, List<RangeBipartition>> entry : localHashMap.entrySet()) {
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
        
        // Convert range bipartitions to frequency map
        convertToFrequencyMap();
        
        return uniqueRangeBipartitions;
    }
    
    /**
     * Extract range bipartitions from a gene tree recursively.
     * This mirrors the original calculateSTBipartitionsUtilLocal logic.
     */
    private void extractRangeBipartitions(TreeNode node, int treeIndex, Map<Object, List<RangeBipartition>> localHashMap) {
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
                
                Object hash = DEFAULT_HASH_FUNCTION.calculateHash(rangeBip, prefixSums, prefixXORs);
                
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
     * Convert range bipartitions to frequency map for unique ones.
     * Uses expensive equality checks if enabled, otherwise trusts hash function uniqueness.
     */
    private void convertToFrequencyMap() {
        System.out.println("Converting unique range bipartitions to frequency map...");
        System.out.println("Expensive equality checks: " + (ENABLE_EXPENSIVE_EQUALITY_CHECKS ? "ENABLED" : "DISABLED (trusting hash)"));
        
        int totalRanges = 0;
        int uniqueRanges = 0;
        
        if (ENABLE_EXPENSIVE_EQUALITY_CHECKS) {
            System.out.println("\n\nPerforming expensive equality checks...\n\n");
            // Original expensive approach with full equality checking
            for (Map.Entry<Object, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
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
                
                // Add unique ranges to frequency map
                for (Map.Entry<RangeBipartition, Integer> rangeEntry : actuallyUniqueRanges.entrySet()) {
                    RangeBipartition range = rangeEntry.getKey();
                    if (range != null) {
                        uniqueRangeBipartitions.merge(range, rangeEntry.getValue(), Integer::sum);
                        uniqueRanges++;
                    }
                }
            }
        } else {
            // Fast approach: trust the hash function, just take the first from each hash group
            for (Map.Entry<Object, List<RangeBipartition>> entry : hashToBipartitions.entrySet()) {
                List<RangeBipartition> ranges = entry.getValue();
                totalRanges += ranges.size();
                
                if (ranges.isEmpty()) continue;
                
                // Just take the first range from each hash group and sum up all occurrences
                RangeBipartition representative = ranges.get(0);
                int totalCount = ranges.size(); // All ranges in this hash group are considered identical
                
                if (representative != null) {
                    uniqueRangeBipartitions.merge(representative, totalCount, Integer::sum);
                    uniqueRanges++;
                }
            }
        }
        
        System.out.println("Memory optimization results:");
        System.out.println("  Total range bipartitions processed: " + totalRanges);
        System.out.println("  Unique RangeBipartitions created: " + uniqueRanges);
        System.out.println("  Memory reduction factor: " + ((double) totalRanges / Math.max(1, uniqueRanges)));
        System.out.println("  Final unique RangeBipartitions: " + uniqueRangeBipartitions.size());
    }
    
    /**
     * Check if two range bipartitions are equal by comparing their actual taxon sets.
     */
    private boolean rangesAreEqual(RangeBipartition range1, RangeBipartition range2) {
        // Same tree - compare ranges directly (both orientations)
        if (range1.geneTreeIndex == range2.geneTreeIndex) {
            // Direct match: left1==left2 && right1==right2
            boolean directMatch = range1.leftStart == range2.leftStart && range1.leftEnd == range2.leftEnd &&
                                 range1.rightStart == range2.rightStart && range1.rightEnd == range2.rightEnd;
            
            // Symmetric match: left1==right2 && right1==left2
            boolean symmetricMatch = range1.leftStart == range2.rightStart && range1.leftEnd == range2.rightEnd &&
                                   range1.rightStart == range2.leftStart && range1.rightEnd == range2.leftEnd;
            
            return directMatch || symmetricMatch;
        }
        
        // Different trees - compare hash values for efficiency
        Object hash1 = DEFAULT_HASH_FUNCTION.calculateHash(range1, prefixSums, prefixXORs);
        Object hash2 = DEFAULT_HASH_FUNCTION.calculateHash(range2, prefixSums, prefixXORs);
        
        if (!hash1.equals(hash2)) {
            return false;
        }
        
        // Hash collision possible - compare actual taxon sets
        // Check both orientations: {left1, right1} == {left2, right2} OR {left1, right1} == {right2, left2}
        Set<Integer> left1 = getLeftTaxonSetForRange(range1);
        Set<Integer> right1 = getRightTaxonSetForRange(range1);
        Set<Integer> left2 = getLeftTaxonSetForRange(range2);
        Set<Integer> right2 = getRightTaxonSetForRange(range2);
        
        // Direct match: left1==left2 && right1==right2
        boolean directMatch = left1.equals(left2) && right1.equals(right2);
        
        // Symmetric match: left1==right2 && right1==left2  
        boolean symmetricMatch = left1.equals(right2) && right1.equals(left2);
        
        return directMatch || symmetricMatch;
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
     * Get processing statistics.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Memory-Efficient Bipartition Processing Statistics:\n");
        sb.append("  Gene trees processed: ").append(geneTrees.size()).append("\n");
        sb.append("  Hash function: ").append(DEFAULT_HASH_FUNCTION.getName()).append(" (using hashed taxon IDs)\n");
        sb.append("  Equality checking mode: ").append(ENABLE_EXPENSIVE_EQUALITY_CHECKS ? "EXPENSIVE" : "FAST (trust hash)").append("\n");
        sb.append("  Unique hash groups: ").append(hashToBipartitions.size()).append("\n");
        sb.append("  Final unique RangeBipartitions: ").append(uniqueRangeBipartitions.size()).append("\n");
        
        int totalRanges = hashToBipartitions.values().stream()
                         .mapToInt(List::size)
                         .sum();
        
        sb.append("  Total range bipartitions: ").append(totalRanges).append("\n");
        sb.append("  Memory reduction factor: ").append(String.format("%.2f", 
            (double) totalRanges / Math.max(1, uniqueRangeBipartitions.size()))).append("x\n");
        
        return sb.toString();
    }
    
    // Getters for accessing internal data structures (needed for memory-optimized weight calculation)
    
    /**
     * Get the hash-to-bipartitions mapping for memory-optimized processing.
     * This exposes the range bipartition data for efficient weight calculation.
     */
    public Map<Object, List<RangeBipartition>> getHashToBipartitions() {
        return hashToBipartitions;
    }
    
    /**
     * Get the gene tree taxa orderings for inverse index construction.
     */
    public int[][] getGeneTreeTaxaOrdering() {
        return geneTreeTaxaOrdering;
    }
    
    /**
     * Get the prefix sums arrays for hash computation.
     */
    public long[][] getPrefixSums() {
        return prefixSums;
    }
    
    /**
     * Get the prefix XORs arrays for hash computation.
     */
    public long[][] getPrefixXORs() {
        return prefixXORs;
    }
}

