package core;

import java.util.*;
import tree.Tree;
import tree.TreeNode;

/**
 * Manages missing taxa information for efficient quartet scoring.
 * 
 * For quartet scoring (ASTRAL-style), we need to compute restricted row sums:
 *   a = |A ∩ Sg| (how many taxa of candidate side A are present in gene tree g)
 *   b = |B ∩ Sg|
 * 
 * Instead of iterating over all taxa, we leverage:
 * 1. Pre-computed missing taxa lists per gene tree (usually small)
 * 2. Membership testing via inverse index of the source tree
 * 
 * This reduces complexity from O(n) to O(|missingTaxa|) per (candidate, geneTree) pair.
 */
public class MissingTaxaManager {
    
    // Missing taxa per gene tree: missingTaxa[treeIndex] = list of missing taxon IDs
    private final List<int[]> missingTaxa;
    
    // Present taxa count per gene tree: presentCount[treeIndex] = |Sg|
    private final int[] presentCount;
    
    // Total number of taxa in the species set
    private final int totalTaxa;
    
    // Number of gene trees
    private final int numTrees;
    
    // Statistics
    private int totalMissingTaxa = 0;
    private int maxMissingPerTree = 0;
    private int treesWithMissing = 0;
    
    public MissingTaxaManager(List<Tree> geneTrees, int totalTaxaCount) {
        System.out.println("==== INITIALIZING MISSING TAXA MANAGER ====");
        System.out.println("Gene trees: " + geneTrees.size());
        System.out.println("Total taxa: " + totalTaxaCount);
        
        this.totalTaxa = totalTaxaCount;
        this.numTrees = geneTrees.size();
        this.missingTaxa = new ArrayList<>(numTrees);
        this.presentCount = new int[numTrees];
        
        buildMissingTaxaLists(geneTrees);
        
        System.out.println("Missing taxa manager initialized");
        System.out.println("Trees with missing taxa: " + treesWithMissing + "/" + numTrees);
        System.out.println("Total missing taxa: " + totalMissingTaxa);
        System.out.println("Average missing per tree: " + (totalMissingTaxa / (double) numTrees));
        System.out.println("Max missing in any tree: " + maxMissingPerTree);
        System.out.println("==== MISSING TAXA MANAGER READY ====");
    }
    
    /**
     * Build missing taxa lists for all gene trees.
     * For each gene tree, identify which taxa from the full taxon set are NOT present.
     */
    private void buildMissingTaxaLists(List<Tree> geneTrees) {
        System.out.println("Building missing taxa lists...");
        
        for (int treeIdx = 0; treeIdx < numTrees; treeIdx++) {
            Tree tree = geneTrees.get(treeIdx);
            
            // Collect taxa present in this tree
            Set<Integer> presentTaxa = new HashSet<>();
            collectTaxa(tree.root, presentTaxa);
            
            presentCount[treeIdx] = presentTaxa.size();
            
            // Find missing taxa (in full set but not in this tree)
            List<Integer> missing = new ArrayList<>();
            for (int taxonId = 0; taxonId < totalTaxa; taxonId++) {
                if (!presentTaxa.contains(taxonId)) {
                    missing.add(taxonId);
                }
            }
            
            // Convert to array for efficiency
            int[] missingArray = missing.stream().mapToInt(Integer::intValue).toArray();
            missingTaxa.add(missingArray);
            
            // Update statistics
            if (missingArray.length > 0) {
                treesWithMissing++;
                totalMissingTaxa += missingArray.length;
                maxMissingPerTree = Math.max(maxMissingPerTree, missingArray.length);
            }
            
            // Progress logging for large datasets
            if ((treeIdx + 1) % 500 == 0 || treeIdx == numTrees - 1) {
                System.out.println("Processed " + (treeIdx + 1) + "/" + numTrees + " trees");
            }
        }
    }
    
    /**
     * Recursively collect all taxa IDs in a tree.
     */
    private void collectTaxa(TreeNode node, Set<Integer> taxaSet) {
        if (node.isLeaf()) {
            taxaSet.add(node.taxon.id);
            return;
        }
        
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                collectTaxa(child, taxaSet);
            }
        }
    }
    
    /**
     * Get the number of taxa present in a gene tree.
     * 
     * @param treeIndex Gene tree index
     * @return |Sg| - number of taxa in gene tree g
     */
    public int getPresentCount(int treeIndex) {
        if (treeIndex < 0 || treeIndex >= numTrees) {
            return 0;
        }
        return presentCount[treeIndex];
    }
    
    /**
     * Get the number of missing taxa in a gene tree.
     * 
     * @param treeIndex Gene tree index
     * @return |S - Sg| - number of taxa NOT in gene tree g
     */
    public int getMissingCount(int treeIndex) {
        if (treeIndex < 0 || treeIndex >= numTrees) {
            return totalTaxa;
        }
        return missingTaxa.get(treeIndex).length;
    }
    
    /**
     * Get the list of missing taxa for a gene tree.
     * 
     * @param treeIndex Gene tree index
     * @return Array of taxon IDs that are NOT in the gene tree
     */
    public int[] getMissingTaxa(int treeIndex) {
        if (treeIndex < 0 || treeIndex >= numTrees) {
            return new int[0];
        }
        return missingTaxa.get(treeIndex);
    }
    
    /**
     * Check if a specific taxon is present in a gene tree.
     * This is a linear search through the missing list - use sparingly.
     * 
     * @param treeIndex Gene tree index
     * @param taxonId Taxon to check
     * @return true if taxon is present in the tree
     */
    public boolean isTaxonPresent(int treeIndex, int taxonId) {
        if (treeIndex < 0 || treeIndex >= numTrees) {
            return false;
        }
        
        int[] missing = missingTaxa.get(treeIndex);
        
        // Check if taxon is in the missing list
        for (int m : missing) {
            if (m == taxonId) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * Count how many taxa from a candidate cluster (represented as a range in a source tree)
     * are NOT present in a target gene tree.
     * 
     * This uses the "iterate over missing taxa" approach for efficiency when missing count is small.
     * 
     * @param inverseIndex Inverse index manager for membership testing
     * @param sourceTree Tree index where the candidate range is defined
     * @param rangeStart Start of candidate range (inclusive)
     * @param rangeEnd End of candidate range (exclusive)
     * @param targetTree Gene tree index to check against
     * @return Number of taxa in the range that are NOT in the target tree
     */
    public int countMissingInRange(InverseIndexManager inverseIndex, 
                                   int sourceTree, int rangeStart, int rangeEnd,
                                   int targetTree) {
        if (targetTree < 0 || targetTree >= numTrees) {
            return 0;
        }
        
        int[] missing = missingTaxa.get(targetTree);
        int count = 0;
        
        int[][] sourceOrderings = inverseIndex.getGeneTreeOrderings();
        int[][] targetInverseIdx = inverseIndex.getInverseIndex();
        
        if (sourceTree < 0 || sourceTree >= sourceOrderings.length ||
            sourceOrderings[sourceTree] == null) {
            return 0;
        }
        
        // Iterate over missing taxa in target tree
        // Check if each is within the range in source tree
        for (int missingTaxon : missing) {
            // Check if this taxon is in the source range
            // Use inverse index of SOURCE tree to find position
            if (sourceTree < targetInverseIdx.length && 
                missingTaxon >= 0 && missingTaxon < targetInverseIdx[sourceTree].length) {
                
                int posInSource = targetInverseIdx[sourceTree][missingTaxon];
                
                // Check if taxon exists in source tree AND is within range
                if (posInSource >= 0 && posInSource >= rangeStart && posInSource < rangeEnd) {
                    count++;
                }
            }
        }
        
        return count;
    }
    
    /**
     * Compute the restricted cluster size |A ∩ Sg| for quartet scoring.
     * 
     * Uses: |A ∩ Sg| = |A| - |A \ Sg|
     * where |A \ Sg| = number of taxa in A that are NOT in gene tree g
     * 
     * This leverages the small missing taxa list for efficiency.
     * 
     * @param inverseIndex Inverse index manager
     * @param sourceTree Tree index where cluster A is defined
     * @param rangeStart Start of A's range (inclusive)
     * @param rangeEnd End of A's range (exclusive)
     * @param targetTree Gene tree index g
     * @return |A ∩ Sg|
     */
    public int getRestrictedClusterSize(InverseIndexManager inverseIndex,
                                        int sourceTree, int rangeStart, int rangeEnd,
                                        int targetTree) {
        int clusterSize = rangeEnd - rangeStart;
        int missingInRange = countMissingInRange(inverseIndex, sourceTree, rangeStart, rangeEnd, targetTree);
        return clusterSize - missingInRange;
    }
    
    /**
     * Get statistics about missing taxa distribution.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Missing Taxa Manager Statistics:\n");
        sb.append("  Total taxa: ").append(totalTaxa).append("\n");
        sb.append("  Gene trees: ").append(numTrees).append("\n");
        sb.append("  Trees with missing taxa: ").append(treesWithMissing).append("\n");
        sb.append("  Total missing taxa: ").append(totalMissingTaxa).append("\n");
        sb.append("  Average missing per tree: ").append(String.format("%.2f", totalMissingTaxa / (double) numTrees)).append("\n");
        sb.append("  Max missing in any tree: ").append(maxMissingPerTree).append("\n");
        sb.append("  Taxa coverage: ").append(String.format("%.2f%%", 
                 100.0 * (numTrees * totalTaxa - totalMissingTaxa) / (double)(numTrees * totalTaxa))).append("\n");
        
        return sb.toString();
    }
    
    // Getters
    public int getTotalTaxa() { return totalTaxa; }
    public int getNumTrees() { return numTrees; }
    public int getTotalMissingTaxa() { return totalMissingTaxa; }
    public int getTreesWithMissing() { return treesWithMissing; }
}
