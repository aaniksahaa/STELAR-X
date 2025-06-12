package core;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import preprocessing.GeneTrees;
import tree.STBipartition;
import utils.BitSet;
import utils.Threading;

public class WeightCalculator {
    
    private GeneTrees geneTrees;
    
    public WeightCalculator(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
    }
    
    private double calculateScore(STBipartition stb1, STBipartition stb2) {
        // First configuration: (A|B) with (X|Y)
        // p1 = |A ∩ X|, p2 = |B ∩ Y|
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
        // p1 = |A ∩ Y|, p2 = |B ∩ X|
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
        // Use ConcurrentHashMap for thread-safe operations
        Map<STBipartition, Double> weights = new ConcurrentHashMap<>();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        // Initialize thread pool with available processors
        Threading.startThreading(numThreads);
        
        int chunkSize = (candidates.size() + numThreads - 1) / numThreads;
        CountDownLatch latch = new CountDownLatch(numThreads);
        
        // Create and submit tasks for each thread
        for (int i = 0; i < numThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, candidates.size());
            
            Threading.execute(() -> {
                try {
                    for (int j = startIdx; j < endIdx; j++) {
                        STBipartition candidate = candidates.get(j);
                        double totalScore = 0.0;
                        
                        // Create a local copy of stBipartitions to avoid concurrent modification
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
            // Wait for all threads to complete
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Weight calculation was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        return weights;
    }
} 