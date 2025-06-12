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
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.Memory;

import preprocessing.GeneTrees;
import tree.STBipartition;
import utils.BitSet;
import utils.Config;
import utils.Threading;

public class WeightCalculator {
    
    private GeneTrees geneTrees;
    
    // JNA interface for CUDA library
    public interface WeightCalcLib extends Library {
        WeightCalcLib INSTANCE = Native.load("weight_calc", WeightCalcLib.class);
        
        // Structure to match CUDA Bipartition struct
        @Structure.FieldOrder({"cluster1", "cluster2", "bitsetSize"})
        public static class Bipartition extends Structure {
            public Pointer cluster1;
            public Pointer cluster2;
            public int bitsetSize;
        }
        
        void launchWeightCalculation(
            Bipartition[] candidates,
            Bipartition[] geneTreeBips,
            int[] frequencies,
            double[] weights,
            int numCandidates,
            int numGeneTreeBips,
            int bitsetSize
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
        switch (Config.COMPUTATION_MODE) {
            case CPU_SINGLE:
                return calculateWeightsSingleThread(candidates);
            case CPU_PARALLEL:
                return calculateWeightsMultiThread(candidates);
            case GPU_PARALLEL:
                return calculateWeightsGPU(candidates);
            default:
                throw new IllegalStateException("Unknown computation mode: " + Config.COMPUTATION_MODE);
        }
    }
    
    private Map<STBipartition, Double> calculateWeightsSingleThread(List<STBipartition> candidates) {
        Map<STBipartition, Double> weights = new HashMap<>();
        
        for (STBipartition candidate : candidates) {
            double totalScore = 0.0;
            
            for (Map.Entry<STBipartition, Integer> entry : geneTrees.stBipartitions.entrySet()) {
                STBipartition geneTreeSTB = entry.getKey();
                int frequency = entry.getValue();
                
                double score = calculateScore(candidate, geneTreeSTB);
                totalScore += score * frequency;
            }
            
            weights.put(candidate, totalScore);
        }
        
        return weights;
    }
    
    private Map<STBipartition, Double> calculateWeightsMultiThread(List<STBipartition> candidates) {
        Map<STBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        Threading.startThreading(numThreads);
        
        int chunkSize = (candidates.size() + numThreads - 1) / numThreads;
        CountDownLatch latch = new CountDownLatch(numThreads);
        
        for (int i = 0; i < numThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, candidates.size());
            
            Threading.execute(() -> {
                try {
                    for (int j = startIdx; j < endIdx; j++) {
                        STBipartition candidate = candidates.get(j);
                        double totalScore = 0.0;
                        
                        Map<STBipartition, Integer> localBipartitions = new HashMap<>(geneTrees.stBipartitions);
                        
                        for (Map.Entry<STBipartition, Integer> entry : localBipartitions.entrySet()) {
                            STBipartition geneTreeSTB = entry.getKey();
                            int frequency = entry.getValue();
                            
                            double score = calculateScore(candidate, geneTreeSTB);
                            totalScore += score * frequency;
                        }
                        
                        weights.put(candidate, totalScore);
                    }
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Weight calculation was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        return weights;
    }
    
    private Map<STBipartition, Double> calculateWeightsGPU(List<STBipartition> candidates) {
        try {
            // Convert Java data structures to CUDA-compatible format
            int numCandidates = candidates.size();
            int numGeneTreeBips = geneTrees.stBipartitions.size();
            int bitsetSize = (geneTrees.realTaxaCount + 63) / 64; // Size in 64-bit words
            
            // Allocate arrays for CUDA
            WeightCalcLib.Bipartition[] cudaCandidates = new WeightCalcLib.Bipartition[numCandidates];
            WeightCalcLib.Bipartition[] cudaGeneTreeBips = new WeightCalcLib.Bipartition[numGeneTreeBips];
            int[] frequencies = new int[numGeneTreeBips];
            double[] weights = new double[numCandidates];
            
            // Convert candidates
            for (int i = 0; i < numCandidates; i++) {
                STBipartition stb = candidates.get(i);
                WeightCalcLib.Bipartition bip = new WeightCalcLib.Bipartition();
                bip.cluster1 = convertBitSetToPointer(stb.cluster1);
                bip.cluster2 = convertBitSetToPointer(stb.cluster2);
                bip.bitsetSize = bitsetSize;
                cudaCandidates[i] = bip;
            }
            
            // Convert gene tree bipartitions
            int idx = 0;
            for (Map.Entry<STBipartition, Integer> entry : geneTrees.stBipartitions.entrySet()) {
                STBipartition stb = entry.getKey();
                WeightCalcLib.Bipartition bip = new WeightCalcLib.Bipartition();
                bip.cluster1 = convertBitSetToPointer(stb.cluster1);
                bip.cluster2 = convertBitSetToPointer(stb.cluster2);
                bip.bitsetSize = bitsetSize;
                cudaGeneTreeBips[idx] = bip;
                frequencies[idx] = entry.getValue();
                idx++;
            }
            
            // Launch CUDA kernel
            WeightCalcLib.INSTANCE.launchWeightCalculation(
                cudaCandidates,
                cudaGeneTreeBips,
                frequencies,
                weights,
                numCandidates,
                numGeneTreeBips,
                bitsetSize
            );
            
            // Convert results back to Java Map
            Map<STBipartition, Double> result = new HashMap<>();
            for (int i = 0; i < numCandidates; i++) {
                result.put(candidates.get(i), weights[i]);
            }
            
            return result;
            
        } catch (Exception e) {
            System.err.println("GPU computation failed, falling back to CPU: " + e.getMessage());
            return calculateWeightsMultiThread(candidates);
        }
    }
    
    private Pointer convertBitSetToPointer(BitSet bitSet) {
        // Convert Java BitSet to native memory
        long[] words = bitSet.toLongArray();
        Pointer ptr = new Memory(words.length * 8); // 8 bytes per long
        ptr.write(0, words, 0, words.length);
        return ptr;
    }
} 