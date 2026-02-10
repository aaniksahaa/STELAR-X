package core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import preprocessing.GeneTrees;
import tree.RangeBipartition;
import tree.MixedBipartition;
import tree.Tree;
import tree.TreeNode;
import utils.BitSet;
import utils.Config;

public class InferenceDP {
    
    private GeneTrees geneTrees;
    private List<RangeBipartition> candidateRangeBips;
    private Map<BitSet, List<RangeBipartition>> clusterToRangeBips;
    private Map<RangeBipartition, Double> rangeBipWeights;
    private Map<BitSet, Double> dpMemo;
    private Map<BitSet, RangeBipartition> dpChoice;
    
    // Memory-optimized implementation (for TRIPLET mode)
    private MemoryOptimizedInferenceDP memoryOptimizedDP;
    
    // ASTRAL-compatible implementation (for QUARTET mode)
    private ASTRALInferenceDP astralDP;
    
    // Flag to track which mode we're using
    private boolean useASTRALMode;
    
    public InferenceDP(GeneTrees geneTrees, List<RangeBipartition> candidateRangeBips) {
        this.geneTrees = geneTrees;
        this.candidateRangeBips = candidateRangeBips;
        this.clusterToRangeBips = new HashMap<>();
        this.dpMemo = new HashMap<>();
        this.dpChoice = new HashMap<>();
        
        // Decide which implementation to use based on scoring mode
        this.useASTRALMode = (Config.SCORING_MODE == Config.ScoringMode.QUARTET);
        
        if (useASTRALMode) {
            System.out.println("==== INITIALIZING INFERENCE DP (ASTRAL-COMPATIBLE) ====");
            System.out.println("Using ASTRAL-style quartet scoring with unrooted tree treatment...");
            this.astralDP = new ASTRALInferenceDP(geneTrees);
        } else {
            System.out.println("==== INITIALIZING INFERENCE DP (MEMORY-OPTIMIZED) ====");
            System.out.println("Using STELAR-X triplet scoring...");
            this.memoryOptimizedDP = new MemoryOptimizedInferenceDP(geneTrees, candidateRangeBips);
        }
    }
    
    /**
     * Enable mixed bipartitions (cross-tree recombination) for extended candidate set.
     * This should be called before solve() if mixed bipartitions are to be used.
     * 
     * @param mixedBips List of mixed bipartitions from CandidateExtender
     */
    public void enableMixedBipartitions(List<MixedBipartition> mixedBips) {
        if (useASTRALMode) {
            // ASTRAL mode handles expansion internally
            System.out.println("Note: ASTRAL mode uses internal expansion, mixed bipartitions ignored");
        } else {
            memoryOptimizedDP.enableMixedBipartitions(mixedBips);
        }
    }
    
    private void preprocessCandidates() {
        // No-op: preprocessing is handled by the respective implementation
    }
    
    private void calculateWeights() {
        // No-op: weight calculation is handled by the respective implementation
    }
    
    public double solve() {
        if (useASTRALMode) {
            System.out.println("Solving using ASTRAL-compatible DP...");
            return astralDP.solve();
        } else {
            System.out.println("Solving using memory-optimized DP...");
            return memoryOptimizedDP.solve();
        }
    }
    
    private double dp(BitSet cluster) {
        // No-op: DP is handled by the respective implementation
        return 0.0;
    }
    
    private boolean isValidPartition(RangeBipartition rangeBip, BitSet cluster) {
        // No-op: validation is handled by the respective implementation
        return true;
    }
    
    public Tree reconstructTree() {
        if (useASTRALMode) {
            System.out.println("Reconstructing tree using ASTRAL-compatible DP...");
            return astralDP.reconstructTree();
        } else {
            System.out.println("Reconstructing tree using memory-optimized DP...");
            return memoryOptimizedDP.reconstructTree();
        }
    }
    
    private TreeNode buildTreeNode(BitSet cluster, Tree tree) {
        // No-op: tree building is handled by the respective implementation
        return null;
    }
    
    public void printDPTable() {
        if (useASTRALMode) {
            System.out.println("ASTRAL DP table printing not yet implemented");
        } else {
            System.out.println("Delegating DP table printing to memory-optimized implementation...");
            memoryOptimizedDP.printDPTable();
        }
    }
} 