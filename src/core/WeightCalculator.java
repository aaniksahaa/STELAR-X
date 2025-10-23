package core;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;

import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
import com.sun.jna.Structure;
import com.sun.jna.Memory;

import preprocessing.GeneTrees;
import tree.STBipartition;
import utils.BitSet;
import utils.Config;
import utils.Threading;

// The parallelism is mainly happening here
public class WeightCalculator {
    
    private GeneTrees geneTrees;
    
    // Updated CUDA interface for memory-optimized approach
    public interface WeightCalcLib extends Library {
        WeightCalcLib INSTANCE = Native.load("weight_calc", WeightCalcLib.class);
        
        // Legacy BitSet-based structure (kept for backward compatibility)
        @Structure.FieldOrder({"cluster1", "cluster2", "bitsetSize"})
        public static class Bipartition extends Structure {
            public Pointer cluster1;
            public Pointer cluster2;
            public int bitsetSize;
            
            public Bipartition() {
                super();
            }
            
            public Bipartition(Pointer p) {
                super(p);
            }
        }
        
        // New compact bipartition structure for memory-optimized approach
        @Structure.FieldOrder({"geneTreeIndex", "leftStart", "leftEnd", "rightStart", "rightEnd"})
        public static class CompactBipartition extends Structure {
            public int geneTreeIndex;
            public int leftStart;
            public int leftEnd;
            public int rightStart;
            public int rightEnd;
            
            public CompactBipartition() {
                super();
            }
            
            public CompactBipartition(Pointer p) {
                super(p);
            }
        }
        
        // Legacy BitSet-based kernel (for fallback)
        void launchWeightCalculation(
            Bipartition[] candidates,
            Bipartition[] geneTreeBips,
            int[] frequencies,
            double[] weights,
            int numCandidates,
            int numGeneTreeBips,
            int bitsetSize
        );
        
        // New memory-optimized kernel using compact ranges and inverse indices
        void launchCompactWeightCalculation(
            CompactBipartition[] candidates,
            CompactBipartition[] geneTreeBips,
            int[] frequencies,
            double[] weights,
            Pointer inverseIndexPtr,    // Flattened [tree*numTaxa + taxon] = position
            Pointer orderingPtr,        // Flattened [tree*numTaxa + position] = taxon
            int numCandidates,
            int numGeneTreeBips,
            int numTrees,
            int numTaxa
        );
    }
    
    public WeightCalculator(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
    }
    
    private double calculateScore(STBipartition stb1, STBipartition stb2) {
        // First configuration: (A|B) with (X|Y)
        BitSet aIntersectX = (BitSet) stb1.cluster1.clone();
        aIntersectX.and(stb2.cluster1);
        int p1 = aIntersectX.cardinality();
        
        BitSet bIntersectY = (BitSet) stb1.cluster2.clone();
        bIntersectY.and(stb2.cluster2);
        int p2 = bIntersectY.cardinality();
        
        double score1 = 0;
        if (p1 + p2 >= 2) {
            score1 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        // Second configuration: (A|B) with (Y|X) - cross configuration
        BitSet aIntersectY = (BitSet) stb1.cluster1.clone();
        aIntersectY.and(stb2.cluster2);
        p1 = aIntersectY.cardinality();
        
        BitSet bIntersectX = (BitSet) stb1.cluster2.clone();
        bIntersectX.and(stb2.cluster1);
        p2 = bIntersectX.cardinality();
        
        double score2 = 0;
        if (p1 + p2 >= 2) {
            score2 = p1 * p2 * (p1 + p2 - 2) / 2.0;
        }
        
        return score1 + score2;
    }
    
    public Map<STBipartition, Double> calculateWeights(List<STBipartition> candidates) {
        System.out.println("==== WEIGHT CALCULATION STARTED (MEMORY-OPTIMIZED) ====");
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Number of candidates: " + (candidates != null ? candidates.size() : "NULL"));
        System.out.println("Gene trees available: " + (geneTrees != null ? "YES" : "NO"));
        if (geneTrees != null) {
            System.out.println("Gene trees stBipartitions count: " + 
                             (geneTrees.stBipartitions != null ? geneTrees.stBipartitions.size() : "NULL"));
            System.out.println("Gene trees realTaxaCount: " + geneTrees.realTaxaCount);
        }
        
        // Use memory-optimized weight calculator for all modes
        System.out.println("Using memory-optimized weight calculation approach");
        MemoryOptimizedWeightCalculator memOptCalculator = new MemoryOptimizedWeightCalculator(geneTrees);
        return memOptCalculator.calculateWeights(candidates);
    }
    
    // Legacy methods removed - now using memory-optimized approach exclusively
}