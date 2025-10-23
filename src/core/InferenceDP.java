package core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import preprocessing.GeneTrees;
import tree.RangeBipartition;
import tree.Tree;
import tree.TreeNode;
import utils.BitSet;

public class InferenceDP {
    
    private GeneTrees geneTrees;
    private List<RangeBipartition> candidateRangeBips;
    private Map<BitSet, List<RangeBipartition>> clusterToRangeBips;
    private Map<RangeBipartition, Double> rangeBipWeights;
    private Map<BitSet, Double> dpMemo;
    private Map<BitSet, RangeBipartition> dpChoice;
    
    // Memory-optimized implementation
    private MemoryOptimizedInferenceDP memoryOptimizedDP;
    
    public InferenceDP(GeneTrees geneTrees, List<RangeBipartition> candidateRangeBips) {
        System.out.println("==== INITIALIZING INFERENCE DP (MEMORY-OPTIMIZED) ====");
        System.out.println("Delegating to memory-optimized implementation...");
        
        // Delegate to memory-optimized implementation
        this.memoryOptimizedDP = new MemoryOptimizedInferenceDP(geneTrees, candidateRangeBips);
        
        // Keep these for compatibility but they won't be used
        this.geneTrees = geneTrees;
        this.candidateRangeBips = candidateRangeBips;
        this.clusterToRangeBips = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
    }
    
    private void preprocessCandidates() {
        // No-op: preprocessing is handled by MemoryOptimizedInferenceDP
    }
    
    private void calculateWeights() {
        // No-op: weight calculation is handled by MemoryOptimizedInferenceDP
    }
    
    public double solve() {
        System.out.println("Solving using memory-optimized DP...");
        return memoryOptimizedDP.solve();
    }
    
    private double dp(BitSet cluster) {
        // No-op: DP is handled by MemoryOptimizedInferenceDP
        return 0.0;
    }
    
    private boolean isValidPartition(RangeBipartition rangeBip, BitSet cluster) {
        // No-op: validation is handled by MemoryOptimizedInferenceDP
        return true;
    }
    
    public Tree reconstructTree() {
        System.out.println("Reconstructing tree using memory-optimized DP...");
        return memoryOptimizedDP.reconstructTree();
    }
    
    private TreeNode buildTreeNode(BitSet cluster, Tree tree) {
        // No-op: tree building is handled by MemoryOptimizedInferenceDP
        return null;
    }
    
    public void printDPTable() {
        System.out.println("Delegating DP table printing to memory-optimized implementation...");
        memoryOptimizedDP.printDPTable();
    }
} 