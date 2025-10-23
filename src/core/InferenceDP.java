package core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import preprocessing.GeneTrees;
import tree.STBipartition;
import tree.Tree;
import tree.TreeNode;
import utils.BitSet;

public class InferenceDP {
    
    private GeneTrees geneTrees;
    private List<STBipartition> candidateSTBips;
    private Map<BitSet, List<STBipartition>> clusterToSTBips;
    private Map<STBipartition, Double> stbipWeights;
    private Map<BitSet, Double> dpMemo;
    private Map<BitSet, STBipartition> dpChoice;
    
    // Memory-optimized implementation
    private MemoryOptimizedInferenceDP memoryOptimizedDP;
    
    public InferenceDP(GeneTrees geneTrees, List<STBipartition> candidateSTBips) {
        System.out.println("==== INITIALIZING INFERENCE DP (MEMORY-OPTIMIZED) ====");
        System.out.println("Delegating to memory-optimized implementation...");
        
        // Delegate to memory-optimized implementation
        this.memoryOptimizedDP = new MemoryOptimizedInferenceDP(geneTrees, candidateSTBips);
        
        // Keep these for compatibility but they won't be used
        this.geneTrees = geneTrees;
        this.candidateSTBips = candidateSTBips;
        this.clusterToSTBips = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
    }
    
    private void preprocessCandidates() {
        for (STBipartition stbip : candidateSTBips) {
            BitSet union = (BitSet) stbip.cluster1.clone();
            union.or(stbip.cluster2);
            
            clusterToSTBips.computeIfAbsent(union, k -> new ArrayList<>()).add(stbip);
        }
    }
    
    private void calculateWeights() {
        WeightCalculator calculator = new WeightCalculator(geneTrees);
        stbipWeights = calculator.calculateWeights(candidateSTBips);
    }
    
    public double solve() {
        System.out.println("Solving using memory-optimized DP...");
        return memoryOptimizedDP.solve();
    }
    
    private double dp(BitSet cluster) {
        if (dpMemo.containsKey(cluster)) {
            return dpMemo.get(cluster);
        }
        
        int taxaCount = cluster.cardinality();
        if (taxaCount <= 2) {
            dpMemo.put(cluster, 0.0);
            return 0.0;
        }
        
        double maxScore = Double.NEGATIVE_INFINITY;
        STBipartition bestChoice = null;
        
        List<STBipartition> candidates = clusterToSTBips.get(cluster);
        if (candidates != null) {
            for (STBipartition stbip : candidates) {
                if (isValidPartition(stbip, cluster)) {
                    double leftScore = dp(stbip.cluster1);
                    double rightScore = dp(stbip.cluster2);
                    double stbipScore = stbipWeights.getOrDefault(stbip, 0.0);
                    
                    double totalScore = leftScore + rightScore + stbipScore;
                    
                    if (totalScore > maxScore) {
                        maxScore = totalScore;
                        bestChoice = stbip;
                    }
                }
            }
        }
        
        if (bestChoice == null) {
            maxScore = 0.0;
        }
        
        dpMemo.put(cluster, maxScore);
        dpChoice.put(cluster, bestChoice);
        return maxScore;
    }
    
    private boolean isValidPartition(STBipartition stbip, BitSet cluster) {
        BitSet union = (BitSet) stbip.cluster1.clone();
        union.or(stbip.cluster2);
        return union.equals(cluster);
    }
    
    public Tree reconstructTree() {
        System.out.println("Reconstructing tree using memory-optimized DP...");
        return memoryOptimizedDP.reconstructTree();
    }
    
    private TreeNode buildTreeNode(BitSet cluster, Tree tree) {
        int taxaCount = cluster.cardinality();
        
        if (taxaCount == 1) {
            int taxonId = cluster.nextSetBit(0);
            return tree.addLeaf(geneTrees.taxa[taxonId]);
        }
        
        STBipartition choice = dpChoice.get(cluster);
        if (choice == null) {
            if (taxaCount == 2) {
                ArrayList<TreeNode> children = new ArrayList<>();
                for (int i = cluster.nextSetBit(0); i >= 0; i = cluster.nextSetBit(i + 1)) {
                    children.add(tree.addLeaf(geneTrees.taxa[i]));
                }
                return tree.addInternalNode(children);
            }
            throw new RuntimeException("No valid bipartition found for cluster: " + cluster + ", taxaCount: " + taxaCount);
        }
        TreeNode leftChild = buildTreeNode(choice.cluster1, tree);
        TreeNode rightChild = buildTreeNode(choice.cluster2, tree);
        ArrayList<TreeNode> children = new ArrayList<>();
        children.add(leftChild);
        children.add(rightChild);
        return tree.addInternalNode(children);
    }
    
    public void printDPTable() {
        System.out.println("DP Results:");
        for (Map.Entry<BitSet, Double> entry : dpMemo.entrySet()) {
            BitSet cluster = entry.getKey();
            double score = entry.getValue();
            STBipartition choice = dpChoice.get(cluster);
            
            System.out.print("Cluster " + cluster + " -> Score: " + score);
            if (choice != null) {
                System.out.println(", Choice: " + choice.print(geneTrees.taxonIdToLabel));
            } else {
                System.out.println(", Choice: none");
            }
        }
    }
} 