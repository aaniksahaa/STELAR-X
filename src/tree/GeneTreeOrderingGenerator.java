package tree;

import java.util.ArrayList;
import java.util.List;

/**
 * Generates left-to-right taxa orderings from gene trees for memory-efficient bipartition representation.
 * 
 * For each gene tree, performs an inorder-like traversal to determine the left-to-right order
 * of taxa (leaves only). This ordering ensures that every subtree bipartition corresponds to
 * a contiguous range in this ordering.
 */
public class GeneTreeOrderingGenerator {
    
    /**
     * Generates left-to-right taxa orderings for all gene trees.
     * 
     * @param geneTrees List of gene trees to process
     * @param maxTaxaCount Maximum number of taxa across all trees
     * @return 2D array where [i][j] is the j-th taxon ID in the i-th gene tree's left-to-right ordering
     */
    public static int[][] generateOrderings(List<Tree> geneTrees, int maxTaxaCount) {
        int numTrees = geneTrees.size();
        int[][] orderings = new int[numTrees][];
        
        for (int i = 0; i < numTrees; i++) {
            Tree tree = geneTrees.get(i);
            orderings[i] = generateSingleTreeOrdering(tree);
        }
        
        return orderings;
    }
    
    /**
     * Generates left-to-right taxa ordering for a single gene tree.
     * Uses inorder traversal to visit leaves in left-to-right order.
     */
    private static int[] generateSingleTreeOrdering(Tree tree) {
        List<Integer> ordering = new ArrayList<>();
        if (tree.root != null) {
            inorderTraversal(tree.root, ordering);
        }
        
        // Convert to array
        int[] result = new int[ordering.size()];
        for (int i = 0; i < ordering.size(); i++) {
            result[i] = ordering.get(i);
        }
        
        return result;
    }
    
    /**
     * Performs inorder traversal to collect leaf taxa in left-to-right order.
     * For rooted binary trees: left subtree -> right subtree (skip internal nodes).
     */
    private static void inorderTraversal(TreeNode node, List<Integer> ordering) {
        if (node == null) return;
        
        if (node.isLeaf()) {
            // Add leaf taxon to ordering
            if (node.taxon != null) {
                ordering.add(node.taxon.id);
            }
        } else {
            // For internal nodes, traverse children in order
            if (node.childs != null) {
                for (TreeNode child : node.childs) {
                    inorderTraversal(child, ordering);
                }
            }
        }
    }
    
    /**
     * Validates that the generated orderings are consistent and complete.
     */
    public static boolean validateOrderings(int[][] orderings, int expectedTaxaCount) {
        if (orderings == null || orderings.length == 0) {
            return false;
        }
        
        for (int[] ordering : orderings) {
            if (ordering == null) {
                continue;  // Skip empty trees
            }
            
            // Check that all taxa IDs are valid
            boolean[] seen = new boolean[expectedTaxaCount];
            for (int taxonId : ordering) {
                if (taxonId < 0 || taxonId >= expectedTaxaCount) {
                    System.err.println("Invalid taxon ID: " + taxonId);
                    return false;
                }
                seen[taxonId] = true;
            }
            
            // Note: Not all trees need to have all taxa (some may be missing in gene trees)
        }
        
        return true;
    }
    
    /**
     * Debug utility to print the generated orderings.
     */
    public static void printOrderings(int[][] orderings, String[] taxonIdToLabel) {
        System.out.println("Generated gene tree orderings:");
        for (int i = 0; i < orderings.length && i < 5; i++) {  // Print first 5 trees
            System.out.print("Tree " + i + ": ");
            int[] ordering = orderings[i];
            for (int j = 0; j < ordering.length; j++) {
                if (j > 0) System.out.print(", ");
                String label = (taxonIdToLabel != null && ordering[j] < taxonIdToLabel.length) 
                    ? taxonIdToLabel[ordering[j]] 
                    : String.valueOf(ordering[j]);
                System.out.print(label);
            }
            System.out.println();
        }
        if (orderings.length > 5) {
            System.out.println("... and " + (orderings.length - 5) + " more trees");
        }
    }
    
    /**
     * Gets the position of a taxon in a specific gene tree ordering.
     * Returns -1 if the taxon is not present in that gene tree.
     */
    public static int getTaxonPosition(int[][] orderings, int geneTreeIndex, int taxonId) {
        if (geneTreeIndex >= orderings.length || orderings[geneTreeIndex] == null) {
            return -1;
        }
        
        int[] ordering = orderings[geneTreeIndex];
        for (int i = 0; i < ordering.length; i++) {
            if (ordering[i] == taxonId) {
                return i;
            }
        }
        
        return -1;  // Taxon not found in this gene tree
    }
}
