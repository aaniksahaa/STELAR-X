package core;

import java.util.*;
import preprocessing.GeneTrees;
import tree.Tree;
import tree.TreeNode;
import utils.BitSet;

/**
 * Counts quartet frequencies for branch support calculation following ASTRAL's methodology.
 * 
 * For each internal branch in the species tree, this class counts how many quartets
 * in the gene trees support each of the three possible topologies around that branch.
 * 
 * Mathematical foundation:
 * - For an internal branch separating clusters A|B from C|D
 * - Count quartets (a,b,c,d) where a∈A, b∈B, c∈C, d∈D
 * - Three possible topologies: ((a,b),(c,d)), ((a,c),(b,d)), ((a,d),(b,c))
 * - f1 = frequency of main topology (species tree topology)
 * - f2, f3 = frequencies of two alternative topologies
 */
public class QuartetFrequencyCounter {
    
    private GeneTrees geneTrees;
    private Tree speciesTree;
    
    /**
     * Data structure to hold quartet frequencies for a branch
     */
    public static class QuartetFrequencies {
        public double f1;  // Main topology frequency (supporting species tree)
        public double f2;  // Alternative topology 1 frequency
        public double f3;  // Alternative topology 2 frequency
        public double n;   // Total effective number of quartets
        public long quartCount; // Total number of possible quartets for this branch
        
        public QuartetFrequencies(double f1, double f2, double f3, double n, long quartCount) {
            this.f1 = f1;
            this.f2 = f2;
            this.f3 = f3;
            this.n = n;
            this.quartCount = quartCount;
        }
        
        @Override
        public String toString() {
            return String.format("f1=%.2f, f2=%.2f, f3=%.2f, n=%.2f, quartets=%d", 
                               f1, f2, f3, n, quartCount);
        }
    }
    
    /**
     * Represents a quartet topology around a branch
     */
    private static class Quadripartition {
        BitSet cluster1, cluster2, cluster3, cluster4;
        
        public Quadripartition(BitSet c1, BitSet c2, BitSet c3, BitSet c4) {
            this.cluster1 = c1;
            this.cluster2 = c2;
            this.cluster3 = c3;
            this.cluster4 = c4;
        }
    }
    
    public QuartetFrequencyCounter(GeneTrees geneTrees, Tree speciesTree) {
        this.geneTrees = geneTrees;
        this.speciesTree = speciesTree;
    }
    
    /**
     * Count quartet frequencies for a specific branch in the species tree.
     * 
     * @param branchNode The node defining the branch (the branch is above this node)
     * @return QuartetFrequencies containing f1, f2, f3, and effective n
     */
    public QuartetFrequencies countQuartetFrequencies(TreeNode branchNode) {
        if (branchNode.isLeaf()) {
            return new QuartetFrequencies(0, 0, 0, 0, 0);
        }
        
        if (branchNode.parent == null) {
            return new QuartetFrequencies(0, 0, 0, 0, 0);
        }
        
        // Get the four clusters around this branch
        Quadripartition quad = getQuadripartitionForBranch(branchNode);
        if (quad == null) {
            return new QuartetFrequencies(0, 0, 0, 0, 0);
        }
        
        // Count total possible quartets for this branch
        long totalQuartets = (long) quad.cluster1.cardinality() * 
                            quad.cluster2.cardinality() * 
                            quad.cluster3.cardinality() * 
                            quad.cluster4.cardinality();
        
        if (totalQuartets == 0) {
            return new QuartetFrequencies(0, 0, 0, 0, 0);
        }
        
        // Count frequencies in gene trees
        double f1 = 0, f2 = 0, f3 = 0;
        double effectiveN = 0;
        
        for (Tree geneTree : geneTrees.geneTrees) {
            QuartetFrequencies geneFreqs = countQuartetsInGeneTree(geneTree, quad);
            f1 += geneFreqs.f1;
            f2 += geneFreqs.f2;
            f3 += geneFreqs.f3;
            effectiveN += geneFreqs.n;
        }
        
        
        return new QuartetFrequencies(f1, f2, f3, effectiveN, totalQuartets);
    }
    
    /**
     * Get the quadripartition (four clusters) around a branch in the species tree
     */
    private Quadripartition getQuadripartitionForBranch(TreeNode branchNode) {
        // For a binary tree, we need to identify the four clusters around the branch
        
        // The branch separates the subtree rooted at branchNode from the rest
        BitSet subtreeCluster = getSubtreeCluster(branchNode);
        BitSet remainingCluster = (BitSet) subtreeCluster.clone();
        remainingCluster.flip(0, geneTrees.realTaxaCount);
        
        // For the subtree, we need two child clusters
        if (branchNode.childs == null) {
            return null;
        }
        
        if (branchNode.childs.size() != 2) {
            return null;
        }
        
        BitSet cluster1 = getSubtreeCluster(branchNode.childs.get(0));
        BitSet cluster2 = getSubtreeCluster(branchNode.childs.get(1));
        
        // For the remaining part, we need to split it appropriately
        // This is more complex and depends on the tree structure
        TreeNode sibling = getSibling(branchNode);
        if (sibling == null) {
            return null;
        }
        
        BitSet cluster3, cluster4;
        if (sibling.childs != null && sibling.childs.size() == 2) {
            cluster3 = getSubtreeCluster(sibling.childs.get(0));
            cluster4 = getSubtreeCluster(sibling.childs.get(1));
        } else {
            // If sibling is a leaf or has different structure
            cluster3 = getSubtreeCluster(sibling);
            cluster4 = new BitSet(geneTrees.realTaxaCount);
            // Get the remaining taxa not in cluster1, cluster2, cluster3
            for (int i = 0; i < geneTrees.realTaxaCount; i++) {
                if (!cluster1.get(i) && !cluster2.get(i) && !cluster3.get(i)) {
                    cluster4.set(i);
                }
            }
        }
        
        // Ensure all clusters are non-empty
        if (cluster1.cardinality() == 0 || cluster2.cardinality() == 0 || 
            cluster3.cardinality() == 0 || cluster4.cardinality() == 0) {
            return null;
        }
        
        return new Quadripartition(cluster1, cluster2, cluster3, cluster4);
    }
    
    /**
     * Get the sibling node of a given node
     */
    private TreeNode getSibling(TreeNode node) {
        if (node.parent == null) {
            return null;
        }
        
        if (node.parent.childs == null) {
            return null;
        }
        
        for (TreeNode child : node.parent.childs) {
            if (child != node) {
                return child;
            }
        }
        
        return null;
    }
    
    /**
     * Get all taxa in the subtree rooted at a node
     */
    private BitSet getSubtreeCluster(TreeNode node) {
        BitSet cluster = new BitSet(geneTrees.realTaxaCount);
        collectSubtreeTaxa(node, cluster);
        return cluster;
    }
    
    /**
     * Recursively collect all taxa in a subtree
     */
    private void collectSubtreeTaxa(TreeNode node, BitSet cluster) {
        if (node.isLeaf()) {
            cluster.set(node.taxon.id);
        } else {
            for (TreeNode child : node.childs) {
                collectSubtreeTaxa(child, cluster);
            }
        }
    }
    
    /**
     * Count quartet frequencies in a single gene tree for a given quadripartition
     */
    private QuartetFrequencies countQuartetsInGeneTree(Tree geneTree, Quadripartition speciesQuad) {
        double f1 = 0, f2 = 0, f3 = 0;
        double n = 0;
        
        // For each possible quartet (one taxon from each cluster)
        for (int a = speciesQuad.cluster1.nextSetBit(0); a >= 0; a = speciesQuad.cluster1.nextSetBit(a + 1)) {
            for (int b = speciesQuad.cluster2.nextSetBit(0); b >= 0; b = speciesQuad.cluster2.nextSetBit(b + 1)) {
                for (int c = speciesQuad.cluster3.nextSetBit(0); c >= 0; c = speciesQuad.cluster3.nextSetBit(c + 1)) {
                    for (int d = speciesQuad.cluster4.nextSetBit(0); d >= 0; d = speciesQuad.cluster4.nextSetBit(d + 1)) {
                        
                        // Check if all four taxa are present in this gene tree
                        if (!geneTree.isTaxonPresent(a) || !geneTree.isTaxonPresent(b) || 
                            !geneTree.isTaxonPresent(c) || !geneTree.isTaxonPresent(d)) {
                            continue;
                        }
                        
                        // Determine the topology of this quartet in the gene tree
                        int topology = getQuartetTopology(geneTree, a, b, c, d);
                        
                        if (topology == 1) {
                            f1++; // ((a,b),(c,d)) - matches species tree
                        } else if (topology == 2) {
                            f2++; // ((a,c),(b,d)) - alternative 1
                        } else if (topology == 3) {
                            f3++; // ((a,d),(b,c)) - alternative 2
                        } else {
                            // topology == 0 means unresolved - should we count this?
                            // For now, let's still count it as effective N but not assign to any topology
                        }
                        
                        n++; // Count as an effective quartet
                    }
                }
            }
        }
        
        return new QuartetFrequencies(f1, f2, f3, n, 0);
    }
    
    /**
     * Determine the topology of a quartet (a,b,c,d) in a gene tree.
     * Returns:
     * 1 for ((a,b),(c,d))
     * 2 for ((a,c),(b,d)) 
     * 3 for ((a,d),(b,c))
     * 0 for unresolved/error
     */
    private int getQuartetTopology(Tree geneTree, int a, int b, int c, int d) {
        // Find the LCA (Lowest Common Ancestor) relationships
        TreeNode lcaAB = findLCA(geneTree, a, b);
        TreeNode lcaCD = findLCA(geneTree, c, d);
        TreeNode lcaAC = findLCA(geneTree, a, c);
        TreeNode lcaBD = findLCA(geneTree, b, d);
        TreeNode lcaAD = findLCA(geneTree, a, d);
        TreeNode lcaBC = findLCA(geneTree, b, c);
        
        // Find the deepest LCA among all pairs
        TreeNode lcaAll = findLCA(geneTree, Arrays.asList(a, b, c, d));
        
        // Determine topology based on LCA relationships
        if (lcaAB != lcaAll && lcaCD != lcaAll && 
            isAncestor(lcaAll, lcaAB) && isAncestor(lcaAll, lcaCD)) {
            return 1; // ((a,b),(c,d))
        } else if (lcaAC != lcaAll && lcaBD != lcaAll && 
                   isAncestor(lcaAll, lcaAC) && isAncestor(lcaAll, lcaBD)) {
            return 2; // ((a,c),(b,d))
        } else if (lcaAD != lcaAll && lcaBC != lcaAll && 
                   isAncestor(lcaAll, lcaAD) && isAncestor(lcaAll, lcaBC)) {
            return 3; // ((a,d),(b,c))
        }
        
        return 0; // Unresolved
    }
    
    /**
     * Find the Lowest Common Ancestor of two taxa in a tree
     */
    private TreeNode findLCA(Tree tree, int taxon1, int taxon2) {
        TreeNode leaf1 = tree.leaves[taxon1];
        TreeNode leaf2 = tree.leaves[taxon2];
        
        if (leaf1 == null || leaf2 == null) {
            return null;
        }
        
        // Get paths to root
        List<TreeNode> path1 = getPathToRoot(leaf1);
        List<TreeNode> path2 = getPathToRoot(leaf2);
        
        // Find LCA
        TreeNode lca = null;
        int minLen = Math.min(path1.size(), path2.size());
        
        for (int i = 0; i < minLen; i++) {
            TreeNode node1 = path1.get(path1.size() - 1 - i);
            TreeNode node2 = path2.get(path2.size() - 1 - i);
            
            if (node1 == node2) {
                lca = node1;
            } else {
                break;
            }
        }
        
        return lca;
    }
    
    /**
     * Find LCA of multiple taxa
     */
    private TreeNode findLCA(Tree tree, List<Integer> taxa) {
        if (taxa.isEmpty()) return null;
        
        TreeNode lca = tree.leaves[taxa.get(0)];
        for (int i = 1; i < taxa.size(); i++) {
            lca = findLCA(tree, lca, tree.leaves[taxa.get(i)]);
        }
        return lca;
    }
    
    /**
     * Find LCA of two nodes
     */
    private TreeNode findLCA(Tree tree, TreeNode node1, TreeNode node2) {
        if (node1 == null || node2 == null) return null;
        
        List<TreeNode> path1 = getPathToRoot(node1);
        List<TreeNode> path2 = getPathToRoot(node2);
        
        TreeNode lca = null;
        int minLen = Math.min(path1.size(), path2.size());
        
        for (int i = 0; i < minLen; i++) {
            TreeNode n1 = path1.get(path1.size() - 1 - i);
            TreeNode n2 = path2.get(path2.size() - 1 - i);
            
            if (n1 == n2) {
                lca = n1;
            } else {
                break;
            }
        }
        
        return lca;
    }
    
    /**
     * Get path from a node to the root
     */
    private List<TreeNode> getPathToRoot(TreeNode node) {
        List<TreeNode> path = new ArrayList<>();
        TreeNode current = node;
        
        while (current != null) {
            path.add(current);
            current = current.parent;
        }
        
        return path;
    }
    
    /**
     * Check if ancestor is an ancestor of descendant
     */
    private boolean isAncestor(TreeNode ancestor, TreeNode descendant) {
        if (ancestor == null || descendant == null) return false;
        
        TreeNode current = descendant;
        while (current != null) {
            if (current == ancestor) {
                return true;
            }
            current = current.parent;
        }
        
        return false;
    }
    
    /**
     * Count quartet frequencies for all internal branches in the species tree
     */
    public Map<TreeNode, QuartetFrequencies> countAllBranchFrequencies() {
        Map<TreeNode, QuartetFrequencies> frequencies = new HashMap<>();
        
        // Post-order traversal to process all internal nodes
        countAllBranchFrequenciesUtil(speciesTree.root, frequencies);
        
        return frequencies;
    }
    
    /**
     * Utility method for recursive branch frequency counting
     */
    private void countAllBranchFrequenciesUtil(TreeNode node, Map<TreeNode, QuartetFrequencies> frequencies) {
        if (node.isLeaf()) {
            return;
        }
        
        // Process children first
        for (TreeNode child : node.childs) {
            countAllBranchFrequenciesUtil(child, frequencies);
        }
        
        // Process this node (represents the branch above it)
        if (!node.isRoot()) {
            QuartetFrequencies freqs = countQuartetFrequencies(node);
            frequencies.put(node, freqs);
        }
    }
}
