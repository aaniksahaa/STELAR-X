package core;

import java.util.*;
import tree.Tree;
import tree.TreeNode;

/**
 * Manages inverse permutation arrays for efficient range intersection calculations.
 * 
 * This class builds and maintains:
 * 1. Gene tree orderings: [treeIndex][position] = taxonId
 * 2. Inverse index: [treeIndex][taxonId] = position
 * 
 * Enables O(min(|A|, |B|)) intersection counting instead of O(n) BitSet operations.
 * 
 * The intersection algorithm works by:
 * 1. Choosing the smaller of two ranges for iteration
 * 2. For each taxon in the smaller range, checking if its position in the other tree
 *    falls within the target range using O(1) inverse index lookup
 */
public class InverseIndexManager {
    
    private final int[][] inverseIndex; // [treeIndex][taxonId] = position
    private final int[][] geneTreeOrderings; // [treeIndex][position] = taxonId
    private final int numTrees;
    private final int numTaxa;
    
    // Statistics for logging
    private long totalIntersectionCalls = 0;
    private long totalElementsProcessed = 0;
    private long maxRangeSize = 0;
    private long minRangeSize = Long.MAX_VALUE;
    
    public InverseIndexManager(List<Tree> geneTrees, int realTaxaCount) {
        System.out.println("==== INITIALIZING INVERSE INDEX MANAGER ====");
        System.out.println("Number of gene trees: " + geneTrees.size());
        System.out.println("Number of taxa: " + realTaxaCount);
        
        this.numTrees = geneTrees.size();
        this.numTaxa = realTaxaCount;
        this.inverseIndex = new int[numTrees][numTaxa];
        this.geneTreeOrderings = new int[numTrees][];
        
        long startTime = System.currentTimeMillis();
        buildInverseIndex(geneTrees);
        long endTime = System.currentTimeMillis();
        
        System.out.println("Inverse index construction completed in " + (endTime - startTime) + " ms");
        System.out.println("Memory allocated: " + 
                         ((long) numTrees * numTaxa * 2 * 4) / (1024 * 1024) + " MB for index arrays");
        System.out.println("==== INVERSE INDEX MANAGER READY ====");
    }
    
    /**
     * Build inverse index from gene trees using left-to-right leaf ordering.
     * This mirrors the approach used in MemoryEfficientBipartitionManager.
     */
    private void buildInverseIndex(List<Tree> geneTrees) {
        System.out.println("Building inverse index mappings...");
        
        int processedTrees = 0;
        int totalLeaves = 0;
        
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            Tree tree = geneTrees.get(treeIdx);
            
            // Get left-to-right ordering (same as MemoryEfficientBipartitionManager)
            List<Integer> ordering = new ArrayList<>();
            collectLeavesInOrder(tree.root, ordering);
            
            geneTreeOrderings[treeIdx] = ordering.stream().mapToInt(Integer::intValue).toArray();
            totalLeaves += ordering.size();
            
            // Build inverse mapping: taxonId -> position
            for (int pos = 0; pos < geneTreeOrderings[treeIdx].length; pos++) {
                int taxonId = geneTreeOrderings[treeIdx][pos];
                if (taxonId >= 0 && taxonId < numTaxa) {
                    inverseIndex[treeIdx][taxonId] = pos;
                } else {
                    System.err.println("WARNING: Invalid taxon ID " + taxonId + 
                                     " in tree " + treeIdx + " at position " + pos);
                }
            }
            
            processedTrees++;
            
            // Log progress for large datasets
            if (processedTrees % 100 == 0 || processedTrees == numTrees) {
                System.out.println("Processed " + processedTrees + "/" + numTrees + 
                                 " trees, average leaves per tree: " + 
                                 (totalLeaves / (double) processedTrees));
            }
        }
        
        System.out.println("Inverse index built successfully");
        System.out.println("Total leaves processed: " + totalLeaves);
        System.out.println("Average leaves per tree: " + (totalLeaves / (double) numTrees));
        
        // Validate index consistency
        validateInverseIndex();
    }
    
    /**
     * Collect leaves in left-to-right order (inorder traversal).
     * Same implementation as MemoryEfficientBipartitionManager for consistency.
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
     * Validate that inverse index is consistent with gene tree orderings.
     */
    private void validateInverseIndex() {
        System.out.println("Validating inverse index consistency...");
        
        int validationErrors = 0;
        
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            int[] ordering = geneTreeOrderings[treeIdx];
            
            for (int pos = 0; pos < ordering.length; pos++) {
                int taxonId = ordering[pos];
                int retrievedPos = inverseIndex[treeIdx][taxonId];
                
                if (retrievedPos != pos) {
                    System.err.println("Validation error: tree " + treeIdx + 
                                     ", taxon " + taxonId + 
                                     ", expected pos " + pos + 
                                     ", got pos " + retrievedPos);
                    validationErrors++;
                    
                    if (validationErrors >= 10) {
                        System.err.println("Too many validation errors, stopping validation");
                        break;
                    }
                }
            }
            
            if (validationErrors >= 10) break;
        }
        
        if (validationErrors == 0) {
            System.out.println("Inverse index validation passed");
        } else {
            System.err.println("Inverse index validation failed with " + validationErrors + " errors");
        }
    }
    
    /**
     * Calculate intersection size between two ranges using inverse index.
     * 
     * Complexity: O(min(|range1|, |range2|)) instead of O(n) for BitSet operations.
     * 
     * @param tree1 First tree index
     * @param start1 Start position in first tree (inclusive)
     * @param end1 End position in first tree (exclusive)
     * @param tree2 Second tree index
     * @param start2 Start position in second tree (inclusive)
     * @param end2 End position in second tree (exclusive)
     * @return Number of taxa in the intersection of the two ranges
     */
    public int getRangeIntersectionSize(int tree1, int start1, int end1, 
                                       int tree2, int start2, int end2) {
        // Validate input parameters
        if (tree1 < 0 || tree1 >= numTrees || tree2 < 0 || tree2 >= numTrees) {
            System.err.println("Invalid tree indices: tree1=" + tree1 + ", tree2=" + tree2);
            return 0;
        }
        
        if (start1 < 0 || end1 < start1 || start2 < 0 || end2 < start2) {
            System.err.println("Invalid range parameters: [" + start1 + "," + end1 + ") and [" + start2 + "," + end2 + ")");
            return 0;
        }
        
        // Update statistics
        totalIntersectionCalls++;
        int size1 = end1 - start1;
        int size2 = end2 - start2;
        maxRangeSize = Math.max(maxRangeSize, Math.max(size1, size2));
        minRangeSize = Math.min(minRangeSize, Math.min(size1, size2));
        
        // Choose smaller range for iteration (complexity optimization)
        if (size1 <= size2) {
            totalElementsProcessed += size1;
            return countIntersection(tree1, start1, end1, tree2, start2, end2);
        } else {
            totalElementsProcessed += size2;
            return countIntersection(tree2, start2, end2, tree1, start1, end1);
        }
    }
    
    /**
     * Count intersection by iterating over the smaller range.
     * 
     * @param smallTree Tree index for the smaller range
     * @param smallStart Start position in smaller range
     * @param smallEnd End position in smaller range
     * @param largeTree Tree index for the larger range
     * @param largeStart Start position in larger range
     * @param largeEnd End position in larger range
     * @return Intersection count
     */
    private int countIntersection(int smallTree, int smallStart, int smallEnd,
                                 int largeTree, int largeStart, int largeEnd) {
        int count = 0;
        int[] smallOrdering = geneTreeOrderings[smallTree];
        
        // Validate array bounds
        if (smallOrdering == null) {
            System.err.println("Null ordering for tree " + smallTree);
            return 0;
        }
        
        int maxPos = Math.min(smallEnd, smallOrdering.length);
        
        for (int pos = smallStart; pos < maxPos; pos++) {
            int taxonId = smallOrdering[pos];
            
            // Validate taxon ID
            if (taxonId < 0 || taxonId >= numTaxa) {
                System.err.println("Invalid taxon ID " + taxonId + " at position " + pos + " in tree " + smallTree);
                continue;
            }
            
            int positionInLargeTree = inverseIndex[largeTree][taxonId];
            
            // Check if taxon falls within large tree's range
            if (positionInLargeTree >= largeStart && positionInLargeTree < largeEnd) {
                count++;
            }
        }
        
        return count;
    }
    
    /**
     * Get statistics about intersection calculations for performance monitoring.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Inverse Index Manager Statistics:\n");
        sb.append("  Trees: ").append(numTrees).append("\n");
        sb.append("  Taxa: ").append(numTaxa).append("\n");
        sb.append("  Total intersection calls: ").append(totalIntersectionCalls).append("\n");
        sb.append("  Total elements processed: ").append(totalElementsProcessed).append("\n");
        
        if (totalIntersectionCalls > 0) {
            sb.append("  Average elements per call: ").append(totalElementsProcessed / (double) totalIntersectionCalls).append("\n");
            sb.append("  Min range size: ").append(minRangeSize == Long.MAX_VALUE ? 0 : minRangeSize).append("\n");
            sb.append("  Max range size: ").append(maxRangeSize).append("\n");
        }
        
        return sb.toString();
    }
    
    /**
     * Reset statistics counters.
     */
    public void resetStatistics() {
        totalIntersectionCalls = 0;
        totalElementsProcessed = 0;
        maxRangeSize = 0;
        minRangeSize = Long.MAX_VALUE;
    }
    
    // Getters for external access (e.g., GPU implementation)
    public int[][] getInverseIndex() { 
        return inverseIndex; 
    }
    
    public int[][] getGeneTreeOrderings() { 
        return geneTreeOrderings; 
    }
    
    public int getNumTrees() { 
        return numTrees; 
    }
    
    public int getNumTaxa() { 
        return numTaxa; 
    }
}
